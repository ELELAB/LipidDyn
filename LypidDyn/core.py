#!/usr/bin/env python

"""
-------------------------------------------------------------------------------------------------
__Author__ : Simone Scrima
Group : Computational Biology Laboratory ; Danish Cancer Society Research Center ; CPH ; Denmark 
PLATFORM : Linux-64
Year : 2019 
Version : 0.1
-------------------------------------------------------------------------------------------------

The pipeline exploits a free-available softwares for the automatic computation
of different biophysical parameters for lipids membranes (i.e. GROMACS, FATSLiM).
To be used in the same directory of all the files coming from the raw productive 
simulations (i.e *.xtc,*.gro, *.tpr etc..).
It allows to compute : 
1) 2D density maps 
2) Thickness and Area Per lipid 
3) Order Parameter 
4) "Movements" ( .txt file containing the X Y Z coordinates for the invidual lipids molecule)
 

The method to compute 3) was originally developed by J. Melcr. with the contribution 
from  H. Antila for NMRlipids project and readapted for the purpose of this work
"""



# coding: utf-8

from __future__ import absolute_import
from __future__ import print_function
import MDAnalysis as mda
import numpy as np
import math
import os, sys
from optparse import OptionParser
from collections import OrderedDict
import errno
import subprocess
import argparse
import shutil
from multiprocessing import Pool

# Import GROMACS-related modules
import gromacs
import gromacs.tools as tools
import gromacs.setup as setup
gromacs.environment.flags['capture_output'] = False

# set-up GROMACS
gromacs.config.setup()


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
        input file name
    returns : dictionary 
        with OrderParameters class instances
    """
   # Using ordered dict since it preserves the read-in order. Might come in handy when comparing to experiments.
    ordPars = OrderedDict()
    with open(def_file,"r") as f:
        for line in f.readlines():
            if not line.startswith("#"):
                items = line.split()
                ordPars[items[0]] = OrderParameter(*items)
    return ordPars

#*******************************************************************





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
        self.apl_cutoff = float(apl_cutoff)
        self.thk_cutoff = float(thk_cutoff)
    
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
        if self.thk_cutoff != 2.0:
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
        
        if self.thk_cutoff != 2.0: 
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


###############################################################################################################################################################    
   
    
# class TrajCommandModified:

#     def __init__(self,
#                  trajectory,
#                  topology,
#                  headgroups_ndx_file,
#                  ncore
#                  ):
    
#         self.trajectory = os.path.abspath(trajectory)
#         self.topology = os.path.abspath(topology)
#         self.headgroups_ndx_file = os.path.abspath(headgroups_ndx_file)
#         self.ncore = int(ncore)

#     def traj(self,
#              list_ndx):
        
#         """Execute the traj command of gromacs,
#         which extract x,y,z coordinates of a selection
#         given by the index; in our case we pass lipids
#         residues in loop. 
    
#         Parameters
#         ----------
#         list_ndx : list 
#                 List of index files of lipid residues
#         """

#         for ndx_file in list_ndx:


#             traj = tools.Traj(f = self.trajectory,
#                               s = self.topology,
#                               ox = ndx_file.split('.')[0] + '.xvg',
#                               com = True,
#                               n = ndx_file,
#                               stdout=False,
#                               stderr=False
#                               )
#             traj.run()


#     def leaflets(self):
          
#         """Creates two index files from the true_bilayer.ndx :
#            index_lower_leaflet_res.ndx : index with lipid residues of lower l.
#            index_upper_leaflet_res.ndx : index with lipid residues of upper l.
#         """        

        
       
#         # Lower Leaflet index
#         make_ndx = tools.Make_ndx(f = self.topology, \
#                                   n = 'true_bilayer.ndx',
#                                   o = 'index_lower_leaflet_res.ndx', \
#                                   input = ('del 1',
#                                            'name 0 lower_leaflet', \
#                                            'splitres 0', \
#                                            'del 0',
#                                            'q'),
#                                   stdout=False,
#                                   stderr=False
#                                   )
#         make_ndx.run()


#         # Upper leaflet index

#         make_ndx = tools.Make_ndx(f = self.topology, \
#                                   n = 'true_bilayer.ndx',
#                                   o = 'index_upper_leaflet_res.ndx', \
#                                   input = ('del 0',
#                                            'name 0 upper_leaflet', \
#                                            'splitres 0', \
#                                            'del 0',
#                                            'q'),
#                                   stdout=False,
#                                   stderr=False
#                                       )
#         make_ndx.run()




#     def get_lipids_indexes(self,
#                            leaflet_ndx):

#         """ Method to create, for each single lipid molecule
#         constituting upper and lower leaflet, a corresponding
#         ''.ndx file , from the main ''.ndx file.

#         Parameters
#         ----------
#         leaflet_ndx : ''.ndx file
#                 Upper or lower leaflet ndx
#         """ 
        

#         # Create different txt file as many as the different lipid / res in the ndx file
        
#         data = {}
#         with open(leaflet_ndx, "r") as f:
#             keep_parsing = True
#             for line in f:
#                 if "[" in line:
#                     header = line.strip("\n").strip("[").strip("]").strip(" ").split("_")[2] \
#                     + "_" +  line.strip("\n").strip("[").strip("]").strip(" ").split("_")[3]
#                     keep_parsing = False
#                     header = header.replace('*', '')
#                 else:
#                     if not keep_parsing:
#                         data[header] = line
#                     else:
#                         data[header] += line
#                     keep_parsing = True


#         # Creates all the different files

#         for header in data.keys():
#             with open(header +".ndx", "w") as out:
#                     out.write("[ %s ]\n%s" % (header, data[header]))


#     def get_xvg_lipids(self):
        
#         """Method to produce the different xvg
#         file using the ndxs produced by the 
#         get_lipids_indexes().
#         One xvg for each index, to extract the
#         X,Y,Z coordinates of each single lipid residue
#         Parameter
#         --------- 
#         dir_name : str
#             Name of the folder in which the method operate 
#         """

#         #starting_directory = os.getcwd() 

#         #dirct = starting_directory +'/'+ dir_name # get the path
#         #os.mkdir(dirct)
#         #os.chdir(dirct)
        
#         l_ndx = [i for i in os.listdir() if i.endswith('.ndx') ]
        
#         #for i in os.listdir():
            
#             #if i.endswith('.ndx'): #Execute the traj command for each file in the folder
                
#                 #l_ndx.append(i)
        
#         chunks = [l_ndx[i::self.ncore] for i in range(self.ncore)] # divide in chunks to parallelize
        
#         # Execute the parallelization 
#         pool = Pool(processes = self.ncore)
#         result = pool.map(self.traj, chunks)
#         pool.close()
#         pool.join()

#         # Modifies the xvg to remove the first  25 lines
#         # and substitute the first line with the name of
#         # the lipid residue and its number
#         os.system('for i in *.xvg ;\
#                    do var=`echo $i | sed "s/.xvg//"` ; \
#                    sed -i "1,24d" ${var}.xvg ; \
#                    sed -i "1 s/.*/$var/g" ${var}.xvg;\
#                    done'
#                    )
        

