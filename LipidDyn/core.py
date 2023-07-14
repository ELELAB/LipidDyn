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
import pandas as pd
import logging
from lipyphilic.lib.order_parameter import SCC
import multiprocessing as mp
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.leaflet import LeafletFinder
import re
mp.set_start_method("fork")
 

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
    

    def __init__(self, u, name, resname, atAname, atBname):
        
        """
        Initialization of an instance of this class.
        Parameters
        ----------
        u : mda.Universe instance
                universe to perform selection of the atoms below
        name : str 
                Name of the order parameter, a label
        resname : str
                Name of residue atoms are in
        atAname : str
                Name of the first atom from the definition file
        atBname : str
                Name of the second atom from the definition file

        """

        self.name = name             # Name of the order parameter, a label
        self.resname = resname       # Name of residue atoms are in
        self.atAname = atAname
        self.atBname = atBname
        
        # Trajectory list for storing OPs
        self.traj = []

        # Select atoms from the mda.universe
        self.selection = self._select_atoms(u, resname, atAname, atBname)


    def _select_atoms(self, u, resname, atAname, atBname):

        '''
        Internal function for the constructor. Returns mda.universe selection with
        the provided atoms.
        ----------
        Returns : pairs of atoms selection, split-by residues
        '''

        # This selection format preserves the order of the atoms (atA, atB) independent of their order in the topology
        selection = u.select_atoms("resname {rnm} and name {atA}".format(
                                        rnm=resname, atA=atAname),
                                     "resname {rnm} and name {atB}".format(
                                        rnm=resname, atB=atBname)
                                    ).atoms.split("residue")

        # Sanity check for having only 2 atoms (A & B) selected
        for res in selection:
            if res.n_atoms != 2:
                logging.warning("Selection >> name {atA} {atB} << \
                contains {n_at} atoms, but should contain exactly 2!".format(
                atA=atAname, atB=atBname, n_at=res.n_atoms))

        return selection


    def calc_OP(self, atoms):
        
        # Calculates Order Parameter according to equation
        # S = 1/2 * (3*cos(theta)^2 -1)

        vec = atoms[1].position - atoms[0].position
        d2 = np.square(vec).sum()
        
        if d2 > bond_len_max_sq:
            logging.warning("Atomic distance for atoms \
                {at1} and {at2} in residue no. {resnr} is suspiciously \
                long: {d}!\nPBC removed???".format(at1=atoms[0].name, \
                at2=atoms[1].name, resnr=atoms[0].resid, d=math.sqrt(d2)))
        
        cos2 = vec[2]**2/d2
        S = 0.5*(3.0*cos2-1.0)
        
        return S


    def calc_angle(self, atoms, z_dim=45.0):
        
        # Calculates the angle between the vector and z-axis in degrees
        # no PBC check!
        
        vec = atoms[1].position - atoms[0].position
        d = math.sqrt(np.square(vec).sum())
        cos = vec[2]/d

        # values for the bottom leaflet are inverted so that 
        # they have the same nomenclature as the top leaflet

        cos *= math.copysign(1.0, atoms[0].position[2] - z_dim*0.5)
        try:
            angle = math.degrees(math.acos(cos))
        except ValueError:
            if abs(cos) >= 1.0:
                logging.debug("Cosine is too large = {} --> truncating it to +/-1.0".format(cos))
                cos = math.copysign(1.0, cos)
                angle = math.degrees(math.acos(cos))

        return angle


    # def get_avg_std_OP(self):
        
    #     # Provides average and stddev of all OPs in self.traj
    #     # This method becomes deprecated after the introduction of 
    #     # error estimation
        
    #     # convert to numpy array
    #     return(np.mean(self.traj), np.std(self.traj))


    def calc_avg_std_stem_OP(self):
        
        # Provides average, stddev and standard error of mean for all OPs in self.traj
        self.means = np.mean(self.traj, axis=0)  # mean over frames
        self.avg = np.mean(self.traj)
        self.std = np.std(self.means)
        self.stem = np.std(self.means)/np.sqrt(len(self.means))


def get_OP_cg(u,lipid_dict,lipid_name, sn):
    
    """
    Compute the SCC Order Parameter 
    for CG systems using the Lipyphilic module SCC()
    ---------
    u: MDAnalysis Universe object

    lipid_dict: dict
                parsed dict from yaml file with the sn 
                and C-beads definition for each lipid type
    lipid_name: str
                name of the lipid residue
    """
    
      
    op_dict={} # dict to store the op values

    # for each CC-atoms bond defined
    for bead in lipid_dict.keys():

        # compute the SCC op 
        sn_cc = SCC(u,tail_sel = "resname"+ " "\
                    + lipid_name + " "\
                    + "and name" + " " \
                    + lipid_dict[bead][0]+ " " \
                    + lipid_dict[bead][1])
        
        sn_cc.run()

        # define the array (n_res, n_time) storing the OP values
        array = sn_cc.SCC
        # compute the op_mean over n_res and n_time    
        sn_OP_avg = np.mean(array)
        # get means only over time frames
        array_means = np.mean(array, axis = 1) 
        # compute the std
        sn_OP_std = np.std(array_means)  
        # compute the stem where N = n_res
        sn_OP_stem = sn_OP_std / math.sqrt(len(array_means))
   
        
        # store mean,std,stem into dict
        op_dict[sn + "_"+ bead] = sn_OP_avg, sn_OP_std, sn_OP_stem
        
    return op_dict

def read_trajs_calc_OPs(u, ordPars):
    
    """
    Procedure that reads in trajectories trajs and then goes 
    through every frame and evaluates each Order Parameter "S" 
    from the list of OPs ordPars. Finally providesaverage, stddev and
    standard error of mean for all OPs in self.traj.
    Parameters
    ----------
    u : MDAnalysis.core.universe.Universe object 
        MDAnalysis universe object
    ordPars : list of OrderParameter class instances
        each item in this list describes an Order parameter to be calculated in the trajectory
    """

    # Go through trajectory frame-by-frame and calculate each OP
    # from the list of OPs for each residue separately
    for op in ordPars.values():
        for frame in u.trajectory:
            
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
        
        # calculate average, std dev and std error for the trajectory
        # it is stored as attribute in each OrderParameter instance
        op.calc_avg_std_stem_OP()


def parse_op_input(u, def_file_path):
    
    """
    Reads order parameter definition file to create OrderParameter instances
    File format is as follows (ignores blank lines and lines starting by "#"):
    OP_name    resname    atom1    atom2
    
    Parameters
    ----------
    u : mda.universe
        universe instance of the system to make atom selection
    def_file_path : string
        path to the definition file
    returns : dictionary 
        with OrderParameters class instances, one for each arguments line in input file
    """

    # Using ordered dict since it preserves the read-in order. Might come in handy when comparing experiments.
    ordPars = OrderedDict()
    with open(def_file_path, 'r') as def_file:
        for line in def_file:
            if not line.startswith("#") and line.strip():
                args = line.strip().split()
                ordPars[args[0]] = OrderParameter(u, *args)

    return ordPars



class Density:

    def __init__ (self,
                  universe,
                  ncore,
                  bin_size = 0.02
                  ):

        """Initialize the class 

        Parameters
        -------------
        universe : MDAnalysis.core.universe.Universe object 
            MDAnalysis universe object
        bin_size : float
            Bin size for the grid spacing
        ncore : int 
            Number of cores employed for the calculation    
        """
        self.universe = universe
        self.bin_size = bin_size
        self.ncore = ncore 

        ts = self.universe.trajectory[0] 
        # divide by ten the coordinates to convert in nm 
        self.n1 =  int(round((ts.dimensions[0]*0.1)/self.bin_size)) # X coordinate of box
        self.n2 =  int(round((ts.dimensions[1]*0.1)/self.bin_size)) # Y coordinate of box
       

    def multiprocessing(self,
                        selection):
        
        """Function to parallelize and speed up the 
           calculations of the 2D density maps.


            Parameters
            -------------
            selection : MDAnalysis.core.groups.AtomGroup object
                selection of atoms on which density will be calculated

            Returns:
            ------------
            grid: numpy.ndarray object
                Numpy array containig the normalized density 
                values
        """

        u = self.universe
        ncore = self.ncore

        # divide the frames in chunk      
        chunks= np.linspace(0, 
                            len(u.trajectory), 
                            ncore+1, 
                            dtype=int)

        traj_sliced=[]
        for ix, i in enumerate(chunks[:-1]):
            traj_sliced.append(u.trajectory[i:chunks[ix+1]])
             
 
        # can be shared between processes         
        final_grid = mp.Manager().list()  
        processes = []
        for piece in traj_sliced:
            # Passing the list
            p = mp.Process(target=self.calculate_density,
                           args=(selection,
                                 piece,
                                 final_grid))  
            p.start()
            processes.append(p)
        for p in processes:
            p.join()
            
        # sum all the grids generated in multiple processes   
        final_grid = sum(list(final_grid))
        grid = self.average(final_grid) # normalize
        grid[:,0] /= 10.0
        grid[0,:] /= 10.0
        grid[1:,1:] *= 1000.0
        return(grid)

    def calculate_density(self,
                          selection,
                          trajectory,
                          l_grids):

        """Calculate the density following Gromacs algorithm.
           The computed density is appended to a Manager() 
           type list object from Multiprocessing library to
           speed up the calculation.

        Parameters
        ----------
        selection : MDAnalysis.core.groups.AtomGroup object
                selection of atoms on which density will be calculated
        trajectory : MDAnalysis.coordinates.XTC.XTCReader object
              MDAnalysis trajectory 
        l_grids : Object/list
              Manager() type list object from Multiprocessing library
        
        Returns:
        ----------
        l_grids : Object/list
              Manager() type list object from Multiprocessing library.
        """
        
        grid = np.zeros((self.n1,self.n2)) # grid with shape of n1 and n2 
        n1n2 = self.n1 * self.n2 
        for ts in trajectory: # begin cycle through traj
            
            invcellvol = n1n2 / np.linalg.det(ts.triclinic_dimensions)

            # x coordinate of the atom divided by the x dimension of box
            # find the fraction of box the atom is in on X 
            m1 = selection.positions[:,0]/ts.dimensions[0]
            m1[m1 >= 1.0] -= 1.0
            m1[m1 < 0.0 ] += 1.0
            m1 *= self.n1

            m2 = selection.positions[:,1]/ts.dimensions[1]
            m2[m2 >= 1.0 ] -= 1.0
            m2[m2 <  0.0 ] += 1.0
            m2 *= self.n2

            grid_coords = np.array([m1, m2], dtype=np.int64)

            np.add.at(grid, tuple(grid_coords), invcellvol)
        l_grids.append(grid) 
        return(l_grids)  # return list object for multiprocessing library
             


    def average(self,
                final_grid):
        
        """Time averaging of the density maps array

        Parameters
        ----------
        final_grid : numpy.ndarray object
            Numpy array containing the grid not averaged

        Returns:
        ----------
        grid : numpy.ndarray object 
            Numpy array containing the final values of the 2D density maps
        """

        box1 = 0
        box2 = 0 
        for ts in self.universe.trajectory:
            box1 += (ts.dimensions[0]) # X coordinate
            box2 += (ts.dimensions[1]) # Y coordinate
 
        # normalize grid point by the number of frames        
        grid = np.true_divide(final_grid, len(self.universe.trajectory))  

        #normalize box1 and tick_x; tick_x == tick_y 
        box1 = box1/len(self.universe.trajectory)
        box2 = box2/len(self.universe.trajectory)
        tick_x = [(i*box1)/self.n1 for i in range(self.n1)]
        tick_y = [(i*box2)/self.n2 for i in range(self.n1)]  
        grid = np.insert(grid, 0, tick_x, 0)   # add a row
        tick_y = np.append(0, tick_y)  # add one 0 to make the shape the same
        grid = np.insert(grid, 0, tick_y, 1)# add a column
        grid = np.around(grid, decimals=8)
        return(grid)





class FatslimCommands:
     
    def __init__ (self,
                  trajectory,
                  gro,
                  headgroups_ndx_file,
                  ncore,
                  apl_cutoff,
                  thk_cutoff
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
        self.apl_cutoff = apl_cutoff
        self.thk_cutoff = thk_cutoff
    
    def membranes(self,
                  out_file):

        """Execute the membranes command of fatslim,
        which identifies the upper and lower leaflet of
        a membrane simulation.

        Parameters
        ----------
        out_file : str 
                Name of the output file.
        """

        a = subprocess.call(['fatslim', 'membranes',
                             '-c',self.gro,
                             '-t',self.gro,
                             '-n',self.headgroups_ndx_file,
                             '--nthread', str(self.ncore),
                             '--output-index','bilayer.ndx'],
                             stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL
                            ) 

    def raw_thickness(self,
                      out_file):

        """Execute the thickness command of fatslim,
        which compute the average thickness of lower,
        upper leaflet and the entire membrane in a
        ''.xvg file.
        In this case we compute the raw values for 
        each frame of the trajectory, creating multiple
        files.
    
        Parameters
        ----------
        out_file : str 
                Name of the output file.
        """
        
             
        # the user modified the cut-off
        if self.thk_cutoff != 6.0:
            a = subprocess.call(['fatslim', 'thickness',
                                 '-c',self.gro,
                                 '-n',self.headgroups_ndx_file,
                                 '-t',self.trajectory,
                                 '--nthread',str(self.ncore),
                                 '--export-thickness-raw', out_file,
                                 '--thickness-cutoff',str(self.thk_cutoff)],
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL
                                )
        # default cut-off
        else: 
            a = subprocess.call(['fatslim', 'thickness',
                                 '-c',self.gro,
                                 '-n',self.headgroups_ndx_file,
                                 '-t',self.trajectory,
                                 '--nthread',str(self.ncore),
                                 '--export-thickness-raw', out_file],
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL
                                )


    def thickness(self,
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
        
        if self.thk_cutoff != 6.0: 
            a = subprocess.call(['fatslim', 'thickness',
                                 '-c',self.gro,
                                 '-n',self.headgroups_ndx_file,
                                 '-t',self.trajectory,
                                 '--nthread', str(self.ncore),
                                 '--plot-thickness', out_file,
                                 '--thickness-cutoff',str(self.thk_cutoff)],
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL
                                )
        else:
            a = subprocess.call(['fatslim',
                                 'thickness',
                                 '-c',self.gro,
                                 '-n',self.headgroups_ndx_file,
                                 '-t',self.trajectory,
                                 '--nthread',str(self.ncore),
                                 '--plot-thickness', out_file],
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL
                                )


    def raw_AreaPerLipid(self,
                         out_file):

        """Execute the APL command of fatslim,
        which compute the area per lipid thickness 
        of lower,upper leaflet and the entire membrane
        in a''.xvg file.
        In this case we compute the raw values for 
        each frame of the trajectory, creating multiple
        files.
    
        Parameters
        ----------
        out_file : str 
                Name of the output file.
        """

        if self.apl_cutoff != 3.0:
            a = subprocess.call(['fatslim', 'apl',
                                 '-c',self.gro,
                                 '-n',self.headgroups_ndx_file,
                                 '-t',self.trajectory,
                                 '--nthread',str(self.ncore),
                                 '--export-apl-raw',out_file,
                                 '--apl-cutoff',str(self.apl_cutoff)],
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL
                                )
        else:
            a = subprocess.call(['fatslim', 'apl',
                                 '-c',self.gro,
                                 '-n',self.headgroups_ndx_file,
                                 '-t',self.trajectory,
                                 '--nthread',str(self.ncore),
                                 '--export-apl-raw',out_file],
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL
                                )


    def AreaPerLipid(self,
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

        if self.apl_cutoff != 3.0:
            a = subprocess.call(['fatslim', 'apl',
                                 '-c',self.gro,
                                 '-n',self.headgroups_ndx_file,
                                 '-t',self.trajectory,
                                 '--nthread',str(self.ncore),
                                 '--plot-apl',out_file,
                                 '--apl-cutoff', str(self.apl_cutoff)],
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL
                                )
        else:
            a = subprocess.call(['fatslim', 'apl',
                                 '-c',self.gro,
                                 '-n',self.headgroups_ndx_file,
                                 '-t',self.trajectory,
                                 '--nthread',str(self.ncore),
                                 '--plot-apl',out_file],
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL
                                )

    def SpeciesXVG(self,
                   RawDir,
                   XvgDir,
                   RS_convert,
                   out_file):

        """Execute the APL command of fatslim,
        which compute the average area per lipid
        of lower,upper leaflet and the entire membrane
        in a ''.xvg file along the trajectory.

        Parameters
        ----------
        RawDir : str
                Path to directory with raw files

        XvgDir : str
                Path to directory with xvg files

        RS_convert : df
                Dataframe correlating residnumber and lipid species
        out_file : str
                Name of output file
        """

        # open Raw data
        RawDict = {}
        for filename in os.listdir(RawDir):
            FilePath = os.path.join(RawDir, filename)
            raw_df = pd.read_csv(FilePath)

            # convert residue numbers to lipid species and rename the column
            raw_df['resid'] = raw_df['resid'].map(RS_convert.set_index('resid')['Species'])
            raw_df.rename(columns={'resid':'species'}, inplace=True)

            # add to dictonary with a name indicating frame number
            nFrame = int(re.findall("\d+", filename)[0].lstrip("_"))
            RawDict[nFrame] = raw_df
        RawDict = dict(sorted(RawDict.items()))

        # open average xvg file
        with open(XvgDir) as f:
            xvgFile = f.read()
            xvgFile = xvgFile.split('\n')

        # define either APL or Thickness
        analysis = RawDir.split('/')[-1]

        # create an xvg file for each lipid class. For each class the following is done:
        # iterating trough each raw file to calulate the average property for each frame.
        # a copy of the average xvg file is copied and the content changed with the calulated.
        Species = list(dict.fromkeys(RS_convert['Species']))
        for nSpecies in Species:
            i = -1
            SpeContent = []
            for frame in RawDict:
                # indexes of species in csv (membrane, upper and lower leaflet)
                mAtoms = (RawDict[frame]['species'] == nSpecies)
                lAtoms = (RawDict[frame]['species'] == nSpecies) & (RawDict[frame]['leaflet'] == 'lower leaflet')
                uAtoms = (RawDict[frame]['species'] == nSpecies) & (RawDict[frame]['leaflet'] == 'upper leaflet')

                # average properties in csv (membrane, upper and lower leaflet)
                Membrane = format(RawDict[frame].iloc[: , -1][mAtoms].mean(),'.3f')
                LowLeaf = format(RawDict[frame].iloc[: , -1][lAtoms].mean(),'.3f')
                UppLeaf = format(RawDict[frame].iloc[: , -1][uAtoms].mean(),'.3f')

                # frame time is taken from orignal average xvg
                i = i + 1
                time = xvgFile[15:][i][0:8]

                # combine and edit to fit xvg format and append to List.
                Row = [time, Membrane, LowLeaf, UppLeaf]
                Row = str("    ".join(str(e) for e in Row) + ' ')
                SpeContent.append(Row)

            # The xvg files are saved in species specific folders.
            if not os.path.exists(out_file + '/' + nSpecies):
                try:
                    os.mkdir((out_file + '/' + nSpecies))
                except OSError as exc:
                    if exc.errno != errno.EEXIST:
                        raise

            # avereage xvg content is replaced
            xvgFile[15:] = SpeContent
            with open(out_file + '/' + nSpecies + '/' + analysis + '_' + nSpecies + '.xvg', 'w') as f:
                for line in xvgFile:
                    f.write(line)
                    f.write('\n')


#--------------- MEMBRANE CURVATURE ---------------
"""
MembraneCurvature
=======================================

:Author: Estefania Barreto-Ojeda
:Year: 2021
:Copyright: GNU Public License v3

MembraneCurvature calculates the mean and Gaussian curvature of
surfaces derived from a selection of reference.

Mean curvature is calculated in units of Å :sup:`-1` and Gaussian curvature
in units of Å :sup:`-2`.
"""

# curvature calculations 
def gaussian_curvature(Z):
    """
    Calculate gaussian curvature from Z cloud points.

    Parameters
    ---------
      Z: np.ndarray.
      Multidimensional array of shape (n,n).

    Returns
    -------
    K : np.ndarray.
      The result of gaussian curvature of Z. Returns multidimensional
      array object with values of gaussian curvature of shape `(n, n)`. 
    """

    Zy, Zx = np.gradient(Z)
    Zxy, Zxx = np.gradient(Zx)
    Zyy, _ = np.gradient(Zy)

    K = (Zxx * Zyy - (Zxy ** 2)) / (1 + (Zx ** 2) + (Zy ** 2)) ** 2

    return(K)


def mean_curvature(Z):
        
    """
    Calculates mean curvature from Z cloud points.

    Parameters
    ----------
    Z: np.ndarray.
      Multidimensional array of shape (n,n).

    Returns
    -------
    H : np.ndarray.
      The result of gaussian curvature of Z. Returns multidimensional
      array object with values of gaussian curvature of shape `(n, n)`.
    """

    Zy, Zx = np.gradient(Z)
    Zxy, Zxx = np.gradient(Zx)
    Zyy, _ = np.gradient(Zy)

    H = (1 + Zx**2) * Zyy + (1 + Zy**2) * Zxx - 2 * Zx * Zy * Zxy
    H = H / (2 * (1 + Zx**2 + Zy**2)**(1.5))

    return(H)


# surface calculations   
def derive_surface(atoms, 
                   n_cells_x, 
                   n_cells_y, 
                   max_width_x, 
                   max_width_y
                   ): 
    
    """
    Derive surface from `atom` positions in `AtomGroup`.
    
    Parameters
    ----------
    atoms: AtomGroup.
      AtomGroup of reference selection to define the surface of the membrane.
   
    n_cells_x: int.
      number of cells in the grid of size `max_width_x`, `x` axis.
   
    n_cells_y: int.
      number of cells in the grid of size `max_width_y`, `y` axis.
   
    max_width_x: float.
      Maximum width of simulation box in x axis. (Determined by simulation box dimensions)
   
    max_width_y: float.
      Maximum width of simulation box in y axis. (Determined by simulation box dimensions)
    
    Returns
    -------
    z_coordinates: numpy.ndarray
      Average z-coordinate values. Return Numpy array of floats of
      shape `(n_cells_x, n_cells_y)`.
    """
        
    coordinates = atoms.positions

    return(get_z_surface(coordinates, 
                         n_x_bins = n_cells_x, 
                         n_y_bins = n_cells_y,
                         x_range = (0, max_width_x),
                         y_range = (0, max_width_y)))


def get_z_surface(coordinates, 
                  n_x_bins = 10, 
                  n_y_bins = 10, 
                  x_range = (0, 100),
                  y_range = (0, 100)
                  ):

    """
    Derive surface from distribution of z coordinates in grid.
    
    Parameters
    ----------
    coordinates : numpy.ndarray 
      Coordinates of AtomGroup. Numpy array of shape=(n_atoms, 3).
   
    n_x_bins : int.
      Number of bins in grid in the `x` dimension. 
   
    n_y_bins : int.
      Number of bins in grid in the `y` dimension. 
   
    x_range : tuple of (float, float)
      Range of coordinates (min, max) in the `x` dimension with shape=(2,).
   
    y_range : tuple of (float, float)
      Range of coordinates (min, max) in the `y` dimension with shape=(2,). 

    Returns
    -------
    z_surface: np.ndarray
      Surface derived from set of coordinates in grid of `x_range, y_range` dimensions.
      Returns Numpy array of floats of shape (`n_x_bins`, `n_y_bins`)
    """

    grid_z_coordinates = np.zeros((n_x_bins, n_y_bins))
    grid_norm_unit = np.zeros((n_x_bins, n_y_bins))

    x_factor = n_x_bins / (x_range[1] - x_range[0])
    y_factor = n_y_bins / (y_range[1] - y_range[0])

    x_coords, y_coords, z_coords = coordinates.T

    cell_x_floor = np.floor(x_coords * x_factor).astype(int)
    cell_y_floor = np.floor(y_coords * y_factor).astype(int)

    for l, m, z in zip(cell_x_floor, cell_y_floor, z_coords):
        try: 
            # negative coordinates
            if l < 0 or m < 0:
                logging.warning(("Atom with negative coordinates falls " \
                                "outside grid boundaries. Element " \
                                "({},{}) in grid can't be assigned." \
                                " Skipping atom.").format(l, m))
            else:
                grid_z_coordinates[l, m] += z
                grid_norm_unit[l, m] += 1

        # too large positive coordinates
        except IndexError:
            logging.warning(("Atom coordinates exceed size of grid " \
                            "and element ({},{}) can't be assigned. " \
                            "Maximum (x,y) coordinates must be < ({}, {}). " \
                            "Skipping atom.").format(l, m, x_range[1], y_range[1]))

    z_surface = normalized_grid(grid_z_coordinates, grid_norm_unit)
    return z_surface


def normalized_grid(grid_z_coordinates, grid_norm_unit):
    
    """
    Calculates average `z` coordinates in unit cell.

    Parameters
    ----------
    z_ref: np.array
      Empty array of `(l,m)`
   
    grid_z_coordinates: np.array
      Array of size `(l,m)` with `z` coordinates stored in unit cell.
   
    grid_norm_unit: np.array
      Array of size `(l,m)` with number of atoms in unit cell.
    
    Returns
    -------
    z_surface: np.ndarray
      Normalized `z` coordinates in grid.        
      Returns Numpy array of floats of shape (`n_x_bins`, `n_y_bins`)
    """
    
    grid_norm_unit = np.where(grid_norm_unit > 0, grid_norm_unit, np.nan)
    z_normalized = grid_z_coordinates / grid_norm_unit
    
    return(z_normalized)


#------------------ Own function --------------------
def curvature_data_extraction(mc):
    # All frames = AF, upper and lower leaflets
    AF_surface = mc.results.z_surface
    AF_mean_curvature = mc.results.mean
    AF_gaussian_curvature = mc.results.gaussian

    # Average over frames, upper and lower leaflets
    Avg_surface = mc.results.average_z_surface
    Avg_mean_curvature = mc.results.average_mean
    Avg_gaussian_curvature = mc.results.average_gaussian

    return([AF_surface, 
            AF_mean_curvature, 
            AF_gaussian_curvature, 
            Avg_surface, 
            Avg_mean_curvature, 
            Avg_gaussian_curvature])


class MembraneCurvature(AnalysisBase):

    """
    MembraneCurvature is a tool to calculate membrane curvature.

    Parameters
    ----------
    universe : Universe or AtomGroup
      An MDAnalysis Universe object.
    
    select : str or iterable of str, optional.
      The selection string of an atom selection to use as a reference to derive a surface.
    
    wrap : bool, optional
      Apply coordinate wrapping to pack atoms into the primary unit cell.
    
    n_x_bins : int, optional, default: '100'
      Number of bins in grid in the x dimension.
    
    n_y_bins : int, optional, default: '100'
      Number of bins in grid in the y dimension.
    
    x_range : tuple of (float, float), optional, default: (0, `universe.dimensions[0]`)
      Range of coordinates (min, max) in the x dimension.
    
    y_range : tuple of (float, float), optional, default: (0, `universe.dimensions[1]`)
      Range of coordinates (min, max) in the y dimension.

    Attributes
    ----------
    results.z_surface : ndarray
      Surface derived from atom selection in every frame
      Array of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    
    results.mean_curvature : ndarray
      Mean curvature associated to the surface.
      Array of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    
    results.gaussian_curvature : ndarray
      Gaussian curvature associated to the surface.
      Arrays of shape (`n_frames`, `n_x_bins`, `n_y_bins`
    
    results.average_z_surface : ndarray
      Average of the array elements in `z_surface`.
      Each array has shape (`n_x_bins`, `n_y_bins`)
    
    results.average_mean_curvature : ndarray
      Average of the array elements in `mean_curvature`.
      Each array has shape (`n_x_bins`, `n_y_bins`)
    
    results.average_gaussian_curvature: ndarray
      Average of the array elements in `gaussian_curvature`.
      Each array has shape (`n_x_bins`, `n_y_bins`)
    """

    def __init__(self,
                 universe,
                 select = 'all',
                 n_x_bins = 100,
                 n_y_bins = 100,
                 x_range = None,
                 y_range = None,
                 wrap = True, 
                 **kwargs
                 ):

        super().__init__(universe.universe.trajectory, **kwargs)
        self.ag = universe.select_atoms(select)
        self.wrap = wrap
        self.n_x_bins = n_x_bins
        self.n_y_bins = n_y_bins
        self.x_range = x_range if x_range else (0, universe.dimensions[0])
        self.y_range = y_range if y_range else (0, universe.dimensions[1])

        # Raise if selection doesn't exist
        if len(self.ag) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup.")

        # Only checks the first frame. NPT simulations not properly covered here.
        # Warning message if range doesn't cover entire dimensions of simulation box
        for dim_string, dim_range, num in [('x', self.x_range, 0), 
                                           ('y', self.y_range, 1)]:
            if (dim_range[1] < universe.dimensions[num]):
                logging.warning(f"Grid range in {dim_string} does not cover " \
                                "entire dimensions of simulation box.\n" \
                                "Minimum dimensions must be equal to " \
                                "simulation box.")

        # Apply wrapping coordinates
        if not self.wrap:
            logging.warning("`wrap == False` may result in inaccurate calculation " \
                            "of membrane curvature. Surfaces will be derived from " \
                            "a reduced number of atoms. \n " \
                            " Ignore this warning if your trajectory has " \
                            " rotational/translational fit rotations! ")

    def _prepare(self):
        # Initialize empty np.array with results
        self.results.z_surface = np.full((self.n_frames,
                                          self.n_x_bins,
                                          self.n_y_bins), np.nan)
        self.results.mean = np.full((self.n_frames,
                                     self.n_x_bins,
                                     self.n_y_bins), np.nan)
        self.results.gaussian = np.full((self.n_frames,
                                         self.n_x_bins,
                                         self.n_y_bins), np.nan)

    def _single_frame(self):
        # Apply wrapping coordinates
        if self.wrap:
            self.ag.wrap()
        # Populate a slice with np.arrays of surface, mean, and gaussian per frame
        self.results.z_surface[self._frame_index] = get_z_surface(self.ag.positions,
                                                                  n_x_bins = self.n_x_bins,
                                                                  n_y_bins = self.n_y_bins,
                                                                  x_range = self.x_range,
                                                                  y_range = self.y_range)
        self.results.mean[self._frame_index] = mean_curvature(self.results.z_surface[self._frame_index])
        self.results.gaussian[self._frame_index] = gaussian_curvature(self.results.z_surface[self._frame_index])

    def _conclude(self):
        self.results.average_z_surface = np.nanmean(self.results.z_surface, axis=0)
        self.results.average_mean = np.nanmean(self.results.mean, axis=0)
        self.results.average_gaussian = np.nanmean(self.results.gaussian, axis=0)


class prot_lip_interaction:
    """
    Protein-lipid interaction analysis class
    Class for storing the protein and lipid information to perform the
    interaction analysis

    Parameters
    ----------
    u: MDAnalysis universe class
        Universe containing the system for the analysis

    p: MDAnalysis AtomGroup
        Selection from u containing the protein to analyse

    lipid_resnames: list of str
        Iterable containing the name of the lipid residues to analyse

    int_doms: list of str
        Iterable containing the sections of the protein to analyse. The
        selections may be MDAnalysis selection commands or ranges (e.g.,
        35-48) that will select the range of residue indices
        automatically. If the int_doms is None or empty, the whole
        protein group will be analysed.

    Attributes
    ----------

    compute_DE -> pandas.DataFrame
        calculates the depletion-enrichment factor for each frame for each
        lipid residue

    compte_DE_byres -> pandas.DataFrame
        calculates the depletion-enrichment factor for each frame for each
        lipid residue and for each individual protein residue

    """

    def __init__(self, u, p, lipid_resnames, int_doms=None):
        self.u = u

        # Parse interaction domains
        if int_doms:
            self.prot = mda.AtomGroup([], u)
            for dom in int_doms:
                try:
                    # Default range of residue indices
                    if re.search('^\d+\-\d+$', dom):
                        dom_group = p.select_atoms(f'resindex {dom}')
                    
                    # General MDAnalysis selection command
                    else:
                        dom_group = p.select_atoms(dom)

                    if dom_group:
                        self.prot = self.prot.union(dom_group)
                    else:
                        logging.warning(f'Selection command {dom} resulted in empty group')

                except BaseException as error:
                    logging.warning(f'Selection command {dom} produced error: {error}')

        # Select all the protein if no domains are provided      
        else:
            self.prot = p

        if not self.prot:
            logging.error("Error: couldn't find protein. Plase review protein" \
                          "selection and interaction domains in the configuration file")
        
        self.prot_idxs = np.unique(self.prot.resindices)

        # Store lipid residues
        self.lipid_resnames = lipid_resnames
        self.lips = mda.AtomGroup([], u)
        for lipid in lipid_resnames:
            self.lips = self.lips.union(u.select_atoms(f'resname {lipid}'))


    def compute_DE(self, lipid_resnames, cutoff=6):
        frames = self.u.trajectory
        prot = self.prot
        lips = self.lips

        # Get data dimensions
        L = len(lipid_resnames)
        F = frames.n_frames
        
        # Calculate global specific lipid count
        anywhere = np.empty((1, L), dtype=np.int64)
        for l,lipid in enumerate(lipid_resnames):
            anywhere[0, l] = lips.select_atoms(f'resname {lipid}').n_atoms
        
        # Calculate total lipid count
        total_anywhere = np.sum(anywhere)

        # Calculate specific lipid count around selected protein residues
        around = np.zeros((F, L), dtype=np.int64)
        for f,ts in enumerate(frames):
            # Selection is performed from the lipid group, but taking the protein group as reference
            around_by_resnames = lips.select_atoms(f'around {cutoff} global group prot', prot=prot).groupby('resnames')
            for l,lipid in enumerate(lipid_resnames):
                if lipid in around_by_resnames:
                    around[f, l] = around_by_resnames[lipid].n_atoms

        # Calculate total lipid count around selected protein residues
        total_around = np.sum(around, axis=1, keepdims=True)

        # Apply adjusted depletion-enrichment formula
        around_DE = np.zeros((F, L), dtype=np.float64)
        where = np.array(total_around, dtype=bool)
        np.divide(around, total_around, out=around_DE, where=where)
        total_DE = np.divide((anywhere - around), (total_anywhere - total_around))
        DE = around_DE / total_DE
        
        # Create dataframe
        columns = pd.Index(lipid_resnames, name='lipid')
        indices = pd.Index([ts.time for ts in frames], name='frame')
        DE_df = pd.DataFrame(DE, index=indices, columns=columns)

        return DE_df


    def compute_DE_byres(self, lipid_resnames, cutoff=6):
        frames = self.u.trajectory
        lips = self.lips
        
        # Get data dimensions
        L = len(lipid_resnames)
        P = len(self.prot_idxs)
        F = frames.n_frames
        
        # Calculate global specific lipid count
        anywhere = np.empty((1, 1, L), dtype=np.int64)
        for l,lipid in enumerate(lipid_resnames):
            anywhere[0, 0, l] = lips.select_atoms(f'resname {lipid}').n_atoms

        # Calculate total lipid count
        total_anywhere = np.sum(anywhere)

        # Calculate specific lipid count around selected protein residues
        around = np.zeros((F, P, L), dtype=np.int64)
        for f,ts in enumerate(frames):
            for p,res in enumerate(self.prot_idxs):
                # Selection is performed from the lipid group, but taking the protein residue numbers as reference
                around_by_resnames = lips.select_atoms(f'around {cutoff} global resnum {res}').groupby('resnames')
                for l,lipid in enumerate(lipid_resnames):
                    if lipid in around_by_resnames:
                        around[f, p, l] = around_by_resnames[lipid].n_atoms

        # Calculate total lipid count around selected protein residues
        total_around = np.sum(around, axis=2, keepdims=True)

        # Apply adjusted depletion-enrichment formula
        around_DE = np.zeros((F, P, L), dtype=np.float64)
        where = np.array(total_around, dtype=bool)
        np.divide(around, total_around, out=around_DE, where=where)
        total_DE = np.divide((anywhere - around), (total_anywhere - total_around))
        DE = around_DE / total_DE
        
        # Create dataframe
        columns = pd.MultiIndex.from_product((self.prot_idxs, lipid_resnames), names=('position', 'lipid'))
        indices = pd.Index([ts.time for ts in frames], name='frame')
        DE_df = pd.DataFrame(DE.reshape((F, P*L)), index=indices, columns=columns)

        return DE_df



class lipid_scrambling:
    """
    Lipid scrambling analysis class
    Class for storing the lipid and information to perform the
    scrambling analysis based on the lipid headgroups z-coordinate

    Parameters
    ----------
    u: MDAnalysis universe class
        Universe containing the system for the analysis
    
    g: MDanalysis AtomGroup
        Selection from u with the lipid headgroups of the system. It is used to
        estimate the surface of the membrane.

    lipid_resnames: list of str or None
        Iterable containing the name of the lipid residues to analyse. If None
        or empty list is provided, all resnames in g will be used.

    p: MDAnalysis AtomGroup or None
        Selection from u containing the protein to analyse. If None is provided,
        contact analysis will not be available.

    int_doms: list of str
        Iterable containing the sections of the protein to analyse. The
        selections may be MDAnalysis selection commands or ranges (X-Y)
        that will select the range of residue indices automatically. If
        the int_doms is None or empty, the whole protein group will be
        analysed.

    Attributes
    ----------

    compute_z_data -> np.ndarray
        Estimates membrane surface by linear regression and calculates the
        leaflet z-positions for each target lipid
    
    filter_scrambling -> list
        From the z_data array obtained in compute_z_data, designates as scrambling
        lipids those that fulfill certain requirements depending on input parameters.
        Returns list with information about scrambling lipids.

    """

    def __init__(self, u, g, lipid_resnames=[], p=None, int_doms=None):
        self.u = u

        # Store lipid residues
        self.g = g
        if not lipid_resnames:
            self.heads = g
        else:
            self.heads = mda.AtomGroup([], u)
            for lipid in lipid_resnames:
                self.heads = self.heads.union(g.select_atoms(f'resname {lipid}'))
        
        # Parse interaction domains
        if p and int_doms:
            self.prot = mda.AtomGroup([], u)
            for dom in int_doms:
                try:
                    # Default range of residue indices
                    if re.search('^\d+\-\d+$', dom):
                        dom_group = p.select_atoms(f'resindex {dom}')
                    
                    # General MDAnalysis selection command
                    else:
                        dom_group = p.select_atoms(dom)

                    if dom_group:
                        self.prot = self.prot.union(dom_group)
                    else:
                        logging.warning(f'Selection command {dom} resulted in empty group')

                except BaseException as error:
                    logging.error(f'Selection command {dom} produced error: {error}')
        
        # Select all the protein if no domains are provided
        else:
            self.prot = p

        if self.prot:
            self.prot_idxs = np.unique(self.prot.resindices)
        else:
            self.prot_idxs = np.empty(0)

    
    def compute_z_data(self, reg_params=None):
        """
        Estimates membrane surface through linear regression, and calculates
        leaflet position for each lipid in self.heads. Regression is polynomial
        regression with degree <reg_params>, default 4 (fairly smooth surface).
        This regression is not perfect, but produces satisfactory results if the
        shape is not too irregular, like those with simple grooves. Each leaflet
        is estimated independently.
        The data for regression is taken from the leaflets obtained by MDAnalysis'
        LeafletFinder. This can fail if the leaflets are connected. The function
        will retry with different cutoffs until two big groups are differentiated.
        If the groups become too small before, the program is unable to get the
        groups in that frame and it will be skipped.

        Parameters
        ----------
        reg_params: int
            Polynomial degree for linear regression. E.g., reg_params=2 will use
            [a + x + x^2 + y + y^2 = z] to estimate the surface, a being the
            intercept and x, y and z the atom coordinates.

        """
        if reg_params is None:
            reg_params = 4

        n_g = len(self.g)
        n_lips = len(self.heads)
        ones = np.ones((n_lips, 1))
        powers = np.arange(1, reg_params+1)
        z_data = np.empty((0, 1 + n_lips * 3))
        for i,ts in enumerate(self.u.trajectory):
            repeat = False
            cutoff = 15.0
            min_cutoff = 7.0
            L = LeafletFinder(self.u, self.g, cutoff=cutoff)
            if len(L.groups()) > 1:
                idx1, idx2 = sorted(range(len(L.groups())), key=L.sizes().get, reverse=True)[:2]
                if L.sizes().get(idx1) > 0.8*n_g:
                    repeat = True
            else:
                repeat = True

            while repeat and cutoff >= min_cutoff:
                cutoff -= 2
                L.update(cutoff=cutoff)
                if len(L.groups()) > 1:
                    idx1, idx2 = sorted(range(len(L.groups())), key=L.sizes().get, reverse=True)[:2]
                    if L.sizes().get(idx1) <= 0.8*n_g:
                        repeat = False

            if cutoff < min_cutoff or L.sizes().get(idx1) + L.sizes().get(idx2) < 0.8*n_g:
                logging.warning(f'Not able to differentiate leaflets, skipping frame {i} at {ts.time} ps')
                continue

            leaf1 = L.group(idx1)
            leaf2 = L.group(idx2)
            coords1 = leaf1.positions
            coords2 = leaf2.positions

            ones1 = np.ones((len(leaf1), 1))
            ones2 = np.ones((len(leaf2), 1))
            x1 = np.c_[ones1, np.power(coords1[:, 0:1], powers), np.power(coords1[:, 1:2], powers)]
            x2 = np.c_[ones2, np.power(coords2[:, 0:1], powers), np.power(coords2[:, 1:2], powers)]
            y1 = coords1[:, 2]
            y2 = coords2[:, 2]

            sol1, res1, _, _ = np.linalg.lstsq(x1, y1, rcond=None)
            sol2, res2, _, _ = np.linalg.lstsq(x2, y2, rcond=None)

            coords = self.heads.positions
            x = np.c_[ones, np.power(coords[:,0:1], powers), np.power(coords[:,1:2], powers)]
            z1 = np.matmul(x, sol1)
            z2 = np.matmul(x, sol2)
            z_up = np.max((z1, z2), axis=0)
            z_low = np.min((z1, z2), axis=0)

            z_data = np.r_[z_data, np.r_[ts.time, np.c_[z_up, z_low, coords[:, 2]].flatten()][None, :]]

        return z_data


    def filter_scrambling(self, z_data, upper_z=None, lower_z=None,
                          upper_rel_z=None, lower_rel_z=None,
                          upper_ratio=None, lower_ratio=None,
                          cutoff=None):
        n_frames = z_data.shape[0]
        n_lips = int((z_data.shape[1] - 1) / 3)
        z_data = z_data[:, 1:].reshape((n_frames, n_lips, 3))

        # Define parameters and defaults
        if upper_ratio is None:
            min_frames_up = 1
        else:
            min_frames_up = upper_ratio * n_frames
        if lower_ratio is None:
            min_frames_low = 1
        else:
            min_frames_low = lower_ratio * n_frames
        if upper_z is None:
            upper_z=-np.inf
        if lower_z is None:
            lower_z=np.inf
        if upper_rel_z is None:
            upper_rel_z=-np.inf
        if lower_rel_z is None:
            lower_rel_z=np.inf
        if cutoff is None:
            cutoff = 6

        rel_z = (z_data[:,:,2] - z_data[:,:,1]) / (z_data[:,:,0] - z_data[:,:,1])
        
        # Determine per-frame for each lipid if it is up/low
        is_up = np.logical_and(z_data[:,:,2] > upper_z, rel_z > upper_rel_z)
        is_low = np.logical_and(z_data[:,:,2] < lower_z, rel_z < lower_rel_z)

        # Determine per-lipid if it has been up or low for the minimum time required
        is_up_filt = np.count_nonzero(is_up, axis=0) >= min_frames_up
        is_low_filt = np.count_nonzero(is_low, axis=0) >= min_frames_low

        # Declare as scrambled lipids those that have been both up and low for the minimum time required
        scrmb = np.logical_and(is_up_filt, is_low_filt)
        scrmb_idxs = np.nonzero(scrmb)[0]

        scrmb_list = []
        for i in scrmb_idxs:
            resname = self.heads.resnames[i]
            resid = self.heads.resids[i]
            z_data_head = z_data[:, i]

            # If protein is defined, compute per-frame contacts with headgroup
            # Protein residue indices (probably starting by 0) are stored!!!
            # This allows differentiating between chains
            if self.prot:
                contacts = []
                head = self.heads.select_atoms(f'resid {resid}')
                for f,ts in enumerate(self.u.trajectory):
                    prot_res = np.unique(self.prot.select_atoms(f'around {cutoff} global group head', head=head).resindices)
                    contacts.append(','.join([str(res) for res in prot_res]))
            else:
                contacts = None

            scrmb_list.append([i, resname, resid, z_data_head, contacts])

        return scrmb_list