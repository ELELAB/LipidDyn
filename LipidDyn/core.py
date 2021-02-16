#!/usr/bin/env python


# Copyright (C) 2019, Simone Scrima <simonescrima@gmail.com>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program. If not, see <http://www.gnu.org/licenses/>.


from __future__ import absolute_import
from __future__ import print_function
import MDAnalysis as mda
import numpy as np
import math
import os, sys
from collections import OrderedDict
import errno
import subprocess
import argparse
import shutil
import multiprocessing as mp


#/******************************************************************
#*    Title: calcOrderParameters.py
#*    Author: J. Melcr with contribution of H.Antila
#*    Date: 2018/03/26 
#*    Availability: https://github.com/NMRLipids/MATCH/scripts
#*


# In Angstroms, max distance between atoms for reasonable OP 
# calculation (PBC and sanity check)
bond_len_max=1.5  
bond_len_max_sq=bond_len_max**2

class OrderParameter:
    
    # Class for storing and manipulating
    # order parameter (OP) related metadata (definition, name, ...),
    # trajectories and methods to evaluate OPs.

    # OP definition consist of:
    #    - name of the OP
    #    - residue name
    #    - involved atoms (exactly 2)
    #    + extra: mean, std.dev. & err. estimate 
    #             (using standard error of the means from individual residues)
    #             of the OP (when reading-in an already calculated result)
    

    def __init__(self, name, resname, atom_A_name, atom_B_name, *args):
        
        """Initialization of an instance of this class.
        Parameters
        ----------
        name : str 
                Name of the order parameter, a label
        resname : str
                Name of residue atoms are in
        atom_A_name : str
                XXX
        atom_B_name : str
                XXX

        """

        self.name = name             # Name of the order parameter, a label
        self.resname = resname       # Name of residue atoms are in
        self.atAname = atom_A_name
        self.atBname = atom_B_name

        # variables for error estimate -- standard error of the mean (STEM)
        self.avg   = None   # average/mean value from all residues
        self.means = None   # list of mean values from each individual residue
        self.std   = None   # standard deviation (sqrt(variance))
        self.stem  = None   # STandard Error of the Mean
        
        # Trajectory as list
        self.traj = []  # for storing OPs
        
        for field in self.__dict__:
            if not isinstance(field, str):
                raise UserWarning("provided name >> {} << is not a string! \n \
                Unexpected behaviour might occur.".format(field))
            else:
                if not field.strip():
                    raise RuntimeError("provided name >> {} << is empty! \n \
                    Cannot use empty names for atoms and OP definitions.".format(field))
        
        # extra optional arguments allow setting avg,std values -- suitable for reading-in results of this script
        if len(args) == 2:
            self.avg = args[0]
            self.std = args[1]
        elif len(args) == 3:
            self.avg  = args[0]
            self.std  = args[1]
            self.stem = args[2]
        else:
            if len(args) != 0:
                raise UserWarning("Number of optional positional arguments is {len}, \
                    not 3, 2 or 0. Args: {args}\n Wrong file format?".format(len=len(args), args=args))


    def calc_OP(self, atoms):
        
        # Calculates Order Parameter according to equation
        # S = 1/2 * (3*cos(theta)^2 -1)

        
        vec = atoms[1].position - atoms[0].position
        d2 = np.square(vec).sum()
        
        if d2>bond_len_max_sq:
            raise UserWarning("Atomic distance for atoms \
                {at1} and {at2} in residue no. {resnr} is suspiciously \
                long: {d}!\nPBC removed???".format(at1=atoms[0].name, \
                at2=atoms[1].name, resnr=atoms[0].resid, d=math.sqrt(d2)
                )
                )             
        
        cos2 = vec[2]**2/d2
        S = 0.5*(3.0*cos2-1.0)
        
        return(S)


    def calc_angle(self, atoms, z_dim = 45.0):
        
        # Calculates the angle between the vector and z-axis in degrees
        # no PBC check!
        
        vec = atoms[1].position - atoms[0].position
        d = math.sqrt(np.square(vec).sum())
        cos = vec[2]/d

        
        # values for the bottom leaflet are inverted so that 
        # they have the same nomenclature as the top leaflet

        cos *= math.copysign(1.0, atoms[0].position[2]-z_dim*0.5)
        try:
            angle = math.degrees(math.acos(cos))
        except ValueError:
            if abs(cos)>=1.0:
                print("Cosine is too large = {} --> truncating it to +/-1.0".format(cos))
                cos = math.copysign(1.0, cos)
                angle = math.degrees(math.acos(cos))
        return angle

    @property
    def get_avg_std_OP(self):
        
        # Provides average and stddev of all OPs in self.traj
        # This method becomes deprecated after the introduction of 
        # error estimation
        
        # convert to numpy array
        return(np.mean(self.traj), np.std(self.traj))


    @property
    def get_avg_std_stem_OP(self):
        
        # Provides average, stddev and standard error of mean for 
        # all OPs in self.traj
        
        self.means = np.mean(self.traj, axis=0)
        return( np.mean(self.traj), 
                 np.std(self.means), 
                 np.std(self.means)/np.sqrt(len(self.means)) )  


def read_trajs_calc_OPs(ordPars, top, trajs):
    
    """Procedure that creates MDAnalysis (mda) Universe instance 
    with topology top,reads in trajectories trajs and then goes 
    through every frame and evaluates each Order Parameter "S" 
    from the list of OPs ordPars.
    Parameters
    ----------
    ordPars : list of OrderParameter class instances
       each item in this list describes an Order parameter to be calculated in the trajectory
    top : str
        filename of a top file (e.g. conf.gro)
    trajs : list of strings
        filenames of trajectories
    """

    # read-in topology and trajectory
    mol = mda.Universe(top, trajs)

    # make atom selections for each OP and store it as its attribute for later use with trajectory
    for op in ordPars.values():
        # selection = pairs of atoms, split-by residues
        # this selection format preserves the order of the atoms (atA, atB) independent of their order in the topology
        selection = mol.select_atoms("resname {rnm} and name {atA}".format(
                                        rnm=op.resname, atA=op.atAname),
                                     "resname {rnm} and name {atB}".format(
                                        rnm=op.resname, atB=op.atBname)
                                    ).atoms.split("residue")
        
        for res in selection:
            
            # check if we have only 2 atoms (A & B) selected
            if res.n_atoms != 2:
                print(res.resnames, res.resids)
                for atom in res.atoms:
                    print(atom.name, atom.id)
                raise UserWarning("Selection >> name {atA} {atB} << \
                contains {nat} atoms, but should contain exactly 2!".format(
                atA=op.atAname, atB=op.atBname, nat=res.n_atoms)
                )
        op.selection = selection
      

    # Go through trajectory frame-by-frame
    # and calculate each OP from the list of OPs
    # for each residue separately
    for frame in mol.trajectory:
        for op in ordPars.values():
            
            # temporary list of order parameters for 
            # each individual residue for the given frame
            temp_S = []
            
            for residue in op.selection:
                if "vec" in op.name:
                    S = op.calc_angle(residue, z_dim=frame.dimensions[2])
                else:
                    S = op.calc_OP(residue)
                temp_S.append(S)

            # resulting S-trajectory will be a list of lists
            # so that individual residues can be easily distinquished
            op.traj.append(temp_S)


def parse_op_input(def_file):
    
    """
    parses input file with Order Parameter definitions
    file format is as follows:
    OP_name    resname    atom1    atom2  +extra: OP_mean  OP_std
    (flexible cols)
    fname : string
        string 
    returns : dictionary 
        with OrderParameters class instances
    """
    # Using ordered dict since it preserves the read-in order. Might come in handy when comparing to experiments.
    ordPars = OrderedDict()
    for line in def_file:
        if not line.startswith("#"):
            items = line.split()
            ordPars[items[0]] = OrderParameter(*items)
    return ordPars

#*******************************************************************


class Density:

    def __init__ (self,
                  universe,
                  ):

        """Initialize the class 

        Parameters
        -------------
        universe : object 
            MDAnalysis universe object  
        """

        self.universe = universe
        self.bins = 0.02
        self.ts = self.universe.trajectory[0] 
        # divide by ten the coordinates to convert in nm 
        self.n1 =  int(round((self.ts.dimensions[0]*0.1)/self.bins)) # X coordinate of box
        self.n2 =  int(round((self.ts.dimensions[1]*0.1)/self.bins)) # Y coordinate of box
       

    def run_density(self,
                    selection):

        """Run the density map calculations for the
        chosen selection of the universe.
    
        Parameters
        ----------
        selection : object
            MDAnalysis AtomGroups selections 
        """

        # if the user wants also the raw data

        box1 = 0
        box2 = 0 
        grid = np.zeros((self.n1,self.n2)) # grid with shape of n1 and n2 
        traj = self.universe.trajectory
        for ts in traj: # begin cycle through traj
            # divide by ten the coordinates to convert in nm
            invcellvol = self.n1*self.n2
            invcellvol /= np.linalg.det(ts.triclinic_dimensions*0.1)
            box1 += (ts.dimensions[0]*0.1)
            box2 += (ts.dimensions[1]*0.1)
            for atom in selection: 
                # x coordinate of the atom divided by the x dimension of box
                # find the fraction of box the atom is in on X
                # divide by ten the coordinates to convert in nm
                m1 = (atom.position[0]*0.1)/(ts.dimensions[0]*0.1) 
                if  m1 >= 1 : 
                    m1 -= 1 # pbc maybe subtracting 1 
                if m1 < 0 : 
                    m1 +=1
                m2 = (atom.position[1]*0.1)/(ts.dimensions[1]*0.1)
                if  m2 >= 1 : 
                    m2 -= 1 # pbc maybe subtracting 1 
                if m2 < 0 :
                    m2 +=1
                grid[int(m1*self.n1)][int(m2*self.n2)] +=  invcellvol  
 
        # normalize grid point by the number of frames        
        grid = np.true_divide(grid, len(self.universe.trajectory))  
        
        #normalize box1 and tick_x; tick_x == tick_y 
        box1 = box1/len(self.universe.trajectory)
        box2 = box2/len(self.universe.trajectory)
        tick_x = [(i*box1)/self.n1 for i in range(self.n1)]
        tick_y = [(i*box2)/self.n2 for i in range(self.n1)]  
        grid = np.insert(grid, 0, tick_x, 0)   # add a row
        tick_y = np.append(0, tick_y)  # add one 0 to make the shape the same
        grid = np.insert(grid, 0, tick_y, 1)# add a column
        grid = np.around(grid,decimals=5)
        
        return(grid)
    
    def run_enrichment(self,
                       selection):

        """Run the lipid enrichment calculations 
        for the 
    
        Parameters
        ----------
        selection : object
            MDAnalysis AtomGroups selections 
        """

        d = {}
        # It can be the entire membrane or the upper/lower leaflets
        d["membrane"] = self.run_density(selection.residues.atoms)
        membrane = np.unique(selection.residues.resnames)
        
        # add protein to dictionary of arrays form the universe 
        protein = self.universe.select_atoms("protein")
        d["protein"]= self.run_density(protein)

        # add all the different lipid residues from the previous selection
        for lipid in membrane:
            # select all the lipid residues in that selection
            sel = self.selection.select_atoms("resname " + lipid)
            d[lipid] = self.run_density(sel)

        
        sum_densmap = []
        d1 = {}
        for key in d:
            # ignore all the membrane and the protein
            if key =="membrane" or key == "protein": 
                continue
            # normalize dividing single lipids by the membrane
            array = np.divide(d[key],
                              d["membrane"],
                              out=np.zeros_like(d[key]),
                              where=d["membrane"]!=0)   
            d1[key] = array
            sum_densmap.append(array) # store in densmap
        
        
        # sum all the normalized density maps to have sum_denmaps 
        # and substract the protein array to the sum_densmap
        sum_densmap = np.sum(sum_densmap,axis=0) - d["protein"]
     
        
        d2 = {}
        for key in d1:
            array = np.divide(d1[key],
                          sum_densmap,
                          out=np.zeros_like(d1[key]),
                          where=sum_densmap!=0) # normalize  ignoring 0s  
            # substitute the first column and rows with their original values
            array[:,0] = d[key][:,0] 
            array[0,:] = d[key][0,:]
            d2[key] = array
        return(d2)


class Fatslim:
     
    def __init__ (self,
                  trajectory,
                  gro,
                  headgroups_ndx_file,
                  ncore,
                  apl_cutoff,
                  thk_cutoff,
                  raw
                  ):

        """Initialize the class 

        Parameters
        -------------
        trajectory : file 
            Name of the .xtc file.
        gro : file 
            Name of the gro file 
        headgroups_ndx_file =  file
            Name of the ndx file containing the lipids headgroup
        ncore : str
            Number of core to parallelize the future metods
        apl_cutoff/thk_cutoff : float
            Float for the cutoff of apl or thickness command         
        """

        self.trajectory = os.path.abspath(trajectory)
        self.gro = os.path.abspath(gro)
        self.headgroups_ndx_file = os.path.abspath(headgroups_ndx_file)
        self.ncore = ncore
        self.apl_cutoff = float(apl_cutoff)
        self.thk_cutoff = float(thk_cutoff)
        self.raw = raw

    def run_thickness(self,
                      out_file):

        """Execute the thickness command of fatslim,
        which compute the average thickness of lower,
        upper leaflet and the entire membrane in a
        ''.xvg file along the trajectory.
    
        Parameters
        ----------
        out_file : str 
                Name of the output file.
        """

        # if the user wants also the raw data
        if self.raw:
            cmd = ['fatslim', 'thickness',
                   '-c',self.gro,
                   '-n',self.headgroups_ndx_file,
                   '-t',self.trajectory,
                   '--nthread',str(self.ncore),
                   '--plot-thickness', out_file,
                   '--export-thickness-raw',out_file,
                   '--thickness-cutoff', str(self.thk_cutoff)]                    
        else:
            cmd = ['fatslim', 'thickness',
                   '-c',self.gro,
                   '-n',self.headgroups_ndx_file,
                   '-t',self.trajectory,
                   '--nthread',str(self.ncore),
                   '--plot-thickness', out_file,
                   '--thickness-cutoff', str(self.thk_cutoff)]        

        a = subprocess.call(cmd,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)

    def run_apl(self,
                out_file):

        """Execute the APL command of fatslim,
        which compute the average area per lipid 
        of lower,upper leaflet and the entire membrane 
        in a ''.xvg file along the trajectory.
    
        Parameters
        ----------
        out_file : str 
                Name of the output file.
        """

        if self.raw:
            cmd = ['fatslim', 'apl',
                   '-c',self.gro,
                   '-n',self.headgroups_ndx_file,
                   '-t',self.trajectory,
                   '--nthread',str(self.ncore),
                   '--plot-apl',out_file,
                   '--export-apl-raw',out_file,
                   '--apl-cutoff',str(self.apl_cutoff)]
        else:
            cmd = ['fatslim', 'apl',
                   '-c',self.gro,
                   '-n',self.headgroups_ndx_file,
                   '-t',self.trajectory,
                   '--nthread',str(self.ncore),
                   '--plot-apl',out_file,
                   '--apl-cutoff',str(self.apl_cutoff)]
        
        a = subprocess.call(cmd,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)




