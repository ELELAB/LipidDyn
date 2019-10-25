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

# Import GROMACS-related modules
import gromacs
import gromacs.tools as tools
import gromacs.setup as setup

# set-up GROMACS
gromacs.config.setup()

###############################################################################################################################################################
d = { "DOPC" :("sn1_C2a DOPC C22 H2R\n\
sn1_C2b DOPC C22 H2S\n\
sn1_C3a DOPC C23 H3R\n\
sn1_C3b DOPC C23 H3S\n\
sn1_C4a DOPC C24 H4R\n\
sn1_C4b DOPC C24 H4S\n\
sn1_C5a DOPC C25 H5R\n\
sn1_C5b DOPC C25 H5S\n\
sn1_C6a DOPC C26 H6R\n\
sn1_C6b DOPC C26 H6S\n\
sn1_C7a DOPC C27 H7R\n\
sn1_C7b DOPC C27 H7S\n\
sn1_C8a DOPC C28 H8R\n\
sn1_C8b DOPC C28 H8S\n\
sn1_C9a DOPC C29 H9R\n\
sn1_C10a DOPC C210 H10R\n\
sn1_C11a DOPC C211 H11R\n\
sn1_C11b DOPC C211 H11S\n\
sn1_C12a DOPC C212 H12R\n\
sn1_C12b DOPC C212 H12S\n\
sn1_C13a DOPC C213 H13R\n\
sn1_C13b DOPC C213 H13S\n\
sn1_C14a DOPC C214 H14R\n\
sn1_C14b DOPC C214 H14S\n\
sn1_C15a DOPC C215 H15R\n\
sn1_C15b DOPC C215 H15S\n\
sn1_C16a DOPC C216 H16R\n\
sn1_C16b DOPC C216 H16S\n\
sn1_C17a DOPC C217 H17R\n\
sn1_C17b DOPC C217 H17S\n\
sn1_C18a DOPC C218 H18R\n\
sn1_C18b DOPC C218 H18S\n\
sn2_C2a DOPC C32 H2X\n\
sn2_C2b DOPC C32 H2Y\n\
sn2_C3a DOPC C33 H3X\n\
sn2_C3b DOPC C33 H3Y\n\
sn2_C4a DOPC C34 H4X\n\
sn2_C4b DOPC C34 H4Y\n\
sn2_C5a DOPC C35 H5X\n\
sn2_C5b DOPC C35 H5Y\n\
sn2_C6a DOPC C36 H6X\n\
sn2_C6b DOPC C36 H6Y\n\
sn2_C7a DOPC C37 H7X\n\
sn2_C7b DOPC C37 H7Y\n\
sn2_C8a DOPC C38 H8X\n\
sn2_C8b DOPC C38 H8Y\n\
sn2_C9a DOPC C39 H9X\n\
sn2_C10a DOPC C310 H10X\n\
sn2_C11a DOPC C311 H11X\n\
sn2_C11b DOPC C311 H11Y\n\
sn2_C12a DOPC C312 H12X\n\
sn2_C12b DOPC C312 H12Y\n\
sn2_C13a DOPC C313 H13X\n\
sn2_C13b DOPC C313 H13Y\n\
sn2_C14a DOPC C314 H14X\n\
sn2_C14b DOPC C314 H14Y\n\
sn2_C15a DOPC C315 H15X\n\
sn2_C15b DOPC C315 H15Y\n\
sn2_C16a DOPC C316 H16X\n\
sn2_C16b DOPC C316 H16Y\n\
sn2_C17a DOPC C317 H17X\n\
sn2_C17b DOPC C317 H17Y\n\
sn2_C18a DOPC C318 H18X\n\
sn2_C18b DOPC C318 H18X"),
"POPC":("sn1_C2a POPC C32 H2X\n\
sn1_C2b POPC C32 H2Y\n\
sn1_C3a POPC C33 H3X\n\
sn1_C3b POPC C33 H3Y\n\
sn1_C4a POPC C34 H4X\n\
sn1_C4b POPC C34 H4Y\n\
sn1_C5a POPC C35 H5X\n\
sn1_C5b POPC C35 H5Y\n\
sn1_C6a POPC C36 H6X\n\
sn1_C6b POPC C36 H6Y\n\
sn1_C7a POPC C37 H7X\n\
sn1_C7b POPC C37 H7Y\n\
sn1_C8a POPC C38 H8X\n\
sn1_C8b POPC C38 H8Y\n\
sn1_C9a POPC C39 H9X\n\
sn1_C9b POPC C39 H9Y\n\
sn1_C10a POPC C310 H10X\n\
sn1_C10b POPC C310 H10Y\n\
sn1_C11a POPC C311 H11X\n\
sn1_C11b POPC C311 H11Y\n\
sn1_C12a POPC C312 H12X\n\
sn1_C12b POPC C312 H12Y\n\
sn1_C13a POPC C313 H13X\n\
sn1_C13b POPC C313 H13Y\n\
sn1_C14a POPC C314 H14X\n\
sn1_C14b POPC C314 H14Y\n\
sn1_C15a POPC C315 H15X\n\
sn1_C15b POPC C315 H15Y\n\
sn1_C16a POPC C316 H16X\n\
sn1_C16b POPC C316 H16Y\n\
sn2_C2a POPC C22 H2R\n\
sn2_C2b POPC C22 H2S\n\
sn2_C3a POPC C23 H3R\n\
sn2_C3b POPC C23 H3S\n\
sn2_C4a POPC C24 H4R\n\
sn2_C4b POPC C24 H4S\n\
sn2_C5a POPC C25 H5R\n\
sn2_C5b POPC C25 H5S\n\
sn2_C6a POPC C26 H6R\n\
sn2_C6b POPC C26 H6S\n\
sn2_C7a POPC C27 H7R\n\
sn2_C7b POPC C27 H7S\n\
sn2_C8a POPC C28 H8R\n\
sn2_C8b POPC C28 H8S\n\
sn2_C9a POPC C29 H91\n\
sn2_C10a POPC C210 H101\n\
sn2_C11a POPC C211 H11R\n\
sn2_C11b POPC C211 H11S\n\
sn2_C12a POPC C212 H12R\n\
sn2_C12b POPC C212 H12S\n\
sn2_C13a POPC C213 H13R\n\
sn2_C13b POPC C213 H13S\n\
sn2_C14a POPC C214 H14R\n\
sn2_C14b POPC C214 H14S\n\
sn2_C15a POPC C215 H15R\n\
sn2_C15b POPC C215 H15S\n\
sn2_C16a POPC C216 H16R\n\
sn2_C16b POPC C216 H16S\n\
sn2_C17a POPC C217 H17R\n\
sn2_C17b POPC C217 H17S\n\
sn2_C18a POPC C218 H18R\n\
sn2_C18b POPC C218 H18S"),
"LSM":("sn2_C2a LSM C2F H2F\n\
sn2_C2b LSM C2F H2G\n\
sn2_C3a LSM C3F H3F\n\
sn2_C3b LSM C3F H3G\n\
sn2_C4a LSM C4F H4F\n\
sn2_C4b LSM C4F H4G\n\
sn2_C5a LSM C5F H5F\n\
sn2_C5b LSM C5F H5G\n\
sn2_C6a LSM C6F H6F\n\
sn2_C6b LSM C6F H6G\n\
sn2_C7a LSM C7F H7F\n\
sn2_C7b LSM C7F H7G\n\
sn2_C8a LSM C8F H8F\n\
sn2_C8b LSM C8F H8G\n\
sn2_C9a LSM C9F H9F\n\
sn2_C9b LSM C9F H9G\n\
sn2_C10a LSM C10F H10F\n\
sn2_C10b LSM C10F H10G\n\
sn2_C11a LSM C11F H11F\n\
sn2_C11b LSM C11F H11G\n\
sn2_C12a LSM C12F H12F\n\
sn2_C12b LSM C12F H12G\n\
sn2_C13a LSM C13F H13F\n\
sn2_C13b LSM C13F H13G\n\
sn2_C14a LSM C14F H14F\n\
sn2_C14b LSM C14F H14G\n\
sn2_C15a LSM C15F H15F\n\
sn2_C15b LSM C15F H15G\n\
sn2_C16a LSM C16F H16F\n\
sn2_C16b LSM C16F H16G\n\
sn2_C17a LSM C17F H17F\n\
sn2_C17b LSM C17F H17G\n\
sn2_C18a LSM C18F H18F\n\
sn2_C18b LSM C18F H18G\n\
sn2_C19a LSM C19F H19F\n\
sn2_C19b LSM C19F H19G\n\
sn2_C20a LSM C20F H20F\n\
sn2_C20b LSM C20F H20G\n\
sn2_C21a LSM C21F H21F\n\
sn2_C21b LSM C21F H21G\n\
sn2_C22a LSM C22F H22F\n\
sn2_C22b LSM C22F H22G\n\
sn2_C23a LSM C23F H23F\n\
sn2_C23b LSM C23F H23G\n\
sn2_C24a LSM C24F H24F\n\
sn2_C24b LSM C24F H24G\n\
sn1_C4a LSM C4S H4S\n\
sn1_C5a LSM C5S H5S\n\
sn1_C6a LSM C6S H6S\n\
sn1_C6b LSM C6S H6T\n\
sn1_C7a LSM C7S H7S\n\
sn1_C7b LSM C7S H7T\n\
sn1_C8a LSM C8S H8S\n\
sn1_C8b LSM C8S H8T\n\
sn1_C9a LSM C9S H9S\n\
sn1_C9b LSM C9S H9T\n\
sn1_C10a LSM C10S H10S\n\
sn1_C10b LSM C10S H10T\n\
sn1_C11a LSM C11S H11S\n\
sn1_C11b LSM C11S H11T\n\
sn1_C12a LSM C12S H12S\n\
sn1_C12b LSM C12S H12T\n\
sn1_C13a LSM C13S H13S\n\
sn1_C13b LSM C13S H13T\n\
sn1_C14a LSM C14S H14S\n\
sn1_C14b LSM C14S H14T\n\
sn1_C15a LSM C15S H15S\n\
sn1_C15b LSM C15S H15T\n\
sn1_C16a LSM C16S H16S\n\
sn1_C16b LSM C16S H16T\n\
sn1_C17a LSM C17S H17S\n\
sn1_C17b LSM C17S H17T\n\
sn1_C18a LSM C18S H18S\n\
sn1_C18b LSM C18S H18T"),
"PSM":("sn2_C2a PSM C2F H2F\n\
sn2_C2b PSM C2F H2G\n\
sn2_C3a PSM C3F H3F\n\
sn2_C3b PSM C3F H3G\n\
sn2_C4a PSM C4F H4F\n\
sn2_C4b PSM C4F H4G\n\
sn2_C5a PSM C5F H5F\n\
sn2_C5b PSM C5F H5G\n\
sn2_C6a PSM C6F H6F\n\
sn2_C6b PSM C6F H6G\n\
sn2_C7a PSM C7F H7F\n\
sn2_C7b PSM C7F H7G\n\
sn2_C8a PSM C8F H8F\n\
sn2_C8b PSM C8F H8G\n\
sn2_C9a PSM C9F H9F\n\
sn2_C9b PSM C9F H9G\n\
sn2_C10a PSM C10F H10F\n\
sn2_C10b PSM C10F H10G\n\
sn2_C11a PSM C11F H11F\n\
sn2_C11b PSM C11F H11G\n\
sn2_C12a PSM C12F H12F\n\
sn2_C12b PSM C12F H12G\n\
sn2_C13a PSM C13F H13F\n\
sn2_C13b PSM C13F H13G\n\
sn2_C14a PSM C14F H14F\n\
sn2_C14b PSM C14F H14G\n\
sn2_C15a PSM C15F H15F\n\
sn2_C15b PSM C15F H15G\n\
sn2_C16a PSM C16F H16F\n\
sn2_C16b PSM C16F H16G\n\
sn1_C4a PSM C4S H4S\n\
sn1_C5a PSM C5S H5S\n\
sn1_C6a PSM C6S H6S\n\
sn1_C6b PSM C6S H6T\n\
sn1_C7a PSM C7S H7S\n\
sn1_C7b PSM C7S H7T\n\
sn1_C8a PSM C8S H8S\n\
sn1_C8b PSM C8S H8T\n\
sn1_C9a PSM C9S H9S\n\
sn1_C9b PSM C9S H9T\n\
sn1_C10a PSM C10S H10S\n\
sn1_C10b PSM C10S H10T\n\
sn1_C11a PSM C11S H11S\n\
sn1_C11b PSM C11S H11T\n\
sn1_C12a PSM C12S H12S\n\
sn1_C12b PSM C12S H12T\n\
sn1_C13a PSM C13S H13S\n\
sn1_C13b PSM C13S H13T\n\
sn1_C14a PSM C14S H14S\n\
sn1_C14b PSM C14S H14T\n\
sn1_C15a PSM C15S H15S\n\
sn1_C15b PSM C15S H15T\n\
sn1_C16a PSM C16S H16S\n\
sn1_C16b PSM C16S H16T\n\
sn1_C17a PSM C17S H17S\n\
sn1_C17b PSM C17S H17T\n\
sn1_C18a PSM C18S H18S\n\
sn1_C18b PSM C18S H18T"),
"NSM":("sn1_C2a NSM C2F H2F\n\
sn2_C2b NSM C2F H2G\n\
sn2_C3a NSM C3F H3F\n\
sn2_C3b NSM C3F H3G\n\
sn2_C4a NSM C4F H4F\n\
sn2_C4b NSM C4F H4G\n\
sn2_C5a NSM C5F H5F\n\
sn2_C5b NSM C5F H5G\n\
sn2_C6a NSM C6F H6F\n\
sn2_C6b NSM C6F H6G\n\
sn2_C7a NSM C7F H7F\n\
sn2_C7b NSM C7F H7G\n\
sn2_C8a NSM C8F H8F\n\
sn2_C8b NSM C8F H8G\n\
sn2_C9a NSM C9F H9F\n\
sn2_C9b NSM C9F H9G\n\
sn2_C10a NSM C10F H10F\n\
sn2_C10b NSM C10F H10G\n\
sn2_C11a NSM C11F H11F\n\
sn2_C11b NSM C11F H11G\n\
sn2_C12a NSM C12F H12F\n\
sn2_C12b NSM C12F H12G\n\
sn2_C13a NSM C13F H13F\n\
sn2_C13b NSM C13F H13G\n\
sn2_C14a NSM C14F H14F\n\
sn2_C14b NSM C14F H14G\n\
sn2_C15a NSM C15F H15F\n\
sn2_C16a NSM C16F H16F\n\
sn2_C17a NSM C17F H17F\n\
sn2_C17b NSM C17F H17G\n\
sn2_C18a NSM C18F H18F\n\
sn2_C18b NSM C18F H18G\n\
sn2_C19a NSM C19F H19F\n\
sn2_C19b NSM C19F H19G\n\
sn2_C20a NSM C20F H20F\n\
sn2_C20b NSM C20F H20G\n\
sn2_C21a NSM C21F H21F\n\
sn2_C21b NSM C21F H21G\n\
sn2_C22a NSM C22F H22F\n\
sn2_C22b NSM C22F H22G\n\
sn2_C23a NSM C23F H23F\n\
sn2_C23b NSM C23F H23G\n\
sn2_C24a NSM C24F H24F\n\
sn2_C24b NSM C24F H24G\n\
sn1_C4a NSM C4S H4S\n\
sn1_C5a NSM C5S H5S\n\
sn1_C6a NSM C6S H6S\n\
sn1_C6b NSM C6S H6T\n\
sn1_C7a NSM C7S H7S\n\
sn1_C7b NSM C7S H7T\n\
sn1_C8a NSM C8S H8S\n\
sn1_C8b NSM C8S H8T\n\
sn1_C9a NSM C9S H9S\n\
sn1_C9b NSM C9S H9T\n\
sn1_C10a NSM C10S H10S\n\
sn1_C10b NSM C10S H10T\n\
sn1_C11a NSM C11S H11S\n\
sn1_C11b NSM C11S H11T\n\
sn1_C12a NSM C12S H12S\n\
sn1_C12b NSM C12S H12T\n\
sn1_C13a NSM C13S H13S\n\
sn1_C13b NSM C13S H13T\n\
sn1_C14a NSM C14S H14S\n\
sn1_C14b NSM C14S H14T\n\
sn1_C15a NSM C15S H15S\n\
sn1_C15b NSM C15S H15T\n\
sn1_C16a NSM C16S H16S\n\
sn1_C16b NSM C16S H16T\n\
sn1_C17a NSM C17S H17S\n\
sn1_C17b NSM C17S H17T\n\
sn1_C18a NSM C18S H18S\n\
sn1_C18b NSM C18S H18T")
        }

bond_len_max=1.5  # in Angstroms, max distance between atoms for reasonable OP calculation (PBC and sanity check)
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
        return (np.mean(self.traj), np.std(self.traj))


    @property
    def get_avg_std_stem_OP(self):
        
        # Provides average, stddev and standard error of mean for 
        # all OPs in self.traj
        
        self.means = np.mean(self.traj, axis=0)
        return ( np.mean(self.traj), 
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


def parse_op_input(resname):
    
    """Parses input file with Order Parameter definitions
    file format is as follows:
    OP_name    resname    atom1    atom2  +extra: OP_mean  OP_std
    (flexible cols)
    Parameters
    ----------
    fname : string
        input file name

    returns : dictionary 
        with OrderParameters class instances
    """
    
    # Using ordered dict since it preserves the read-in order. Might come in handy when comparing to experiments.
    ordPars = OrderedDict()
    l = resname.split("\n")
    o = []

    for items in l:
        items = items.split(" ")
        ordPars[items[0]] = OrderParameter(*items)
    return(ordPars)

###############################################################################################################################################################


class FatslimCommands:
     
    def __init__ (self,
                  gro,
                  headgroups_ndx_file,
                  thread,
                  ):

        self.gro = gro
        self.thread = thread
        self.headgroups_ndx_file = headgroups_ndx_file

    
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
                             '--nthreads',self.thread,
                             '--output-index','bilayer.ndx'],
                            ) 

    # def raw_thickness(self,
    #                   trajectory,
    #                   cutoff, 
    #                   out_file):

    #     """Execute the thickness command of fatslim,
    #     which compute the average thickness of lower,
    #     upper leaflet and the entire membrane in a
    #     ''.xvg file.
    #     In this case we compute the raw values for 
    #     each frame of the trajectory, creating multiple
    #     files.
    
    #     Parameters
    #     ----------
    #     trajectory : str 
    #         Name of the .xtc file.
    #     out_file : str 
    #             Name of the output file.
    #     """
    #     cutoff = float(cutoff)
    #     cutoff =str(cutoff)
             
    #     # the user modified the cut-off
    #     if cutoff != 2.0:
    #         a = subprocess.call(['fatslim', 'thickness',
    #                              '-c',self.gro,
    #                              '-n',self.headgroups_ndx_file,
    #                              '-t',trajectory,
    #                              '--nthreads',self.thread,
    #                              '--export-thickness-raw', out_file,
    #                              '--thickness-cutoff',cutoff]
    #                             )
    #     # default cut-off
    #     else: 
    #         a = subprocess.call(['fatslim', 'thickness',
    #                              '-c',self.gro,
    #                              '-n',self.headgroups_ndx_file,
    #                              '-t',trajectory,
    #                              '--nthreads',self.thread,
    #                              '--export-thickness-raw', out_file],
    #                             )


    def thickness(self,
                  trajectory,
                  cutoff, 
                  out_file):

        """Execute the thickness command of fatslim,
        which compute the average thickness of lower,
        upper leaflet and the entire membrane in a
        ''.xvg file along the trajectory.
    
        Parameters
        ----------
        trajectory : str 
            Name of the .xtc file.
        out_file : str 
                Name of the output file.
        """
        cutoff = float(cutoff)
        cutoff =str(cutoff)
        if cutoff != 2.0:
            a = subprocess.call(['fatslim', 'thickness',
                                 '-c',self.gro,
                                 '-n',self.headgroups_ndx_file,
                                 '-t',trajectory,
                                 '--nthreads',self.thread,
                                 '--plot-thickness', out_file,
                                 '--thickness-cutoff',cutoff]
                                )
        else:
            a = subprocess.call(['fatslim',
                                 'thickness',
                                 '-c',self.gro,
                                 '-n',self.headgroups_ndx_file,
                                 '-t',trajectory,
                                 '--nthreads',self.thread,
                                 '--plot-thickness',out_file],
                                )


    # def raw_AreaPerLipid(self,
    #                      trajectory,
    #                      cutoff, 
    #                      out_file):

    #     """Execute the APL command of fatslim,
    #     which compute the area per lipid thickness 
    #     of lower,upper leaflet and the entire membrane
    #     in a''.xvg file.
    #     In this case we compute the raw values for 
    #     each frame of the trajectory, creating multiple
    #     files.
    
    #     Parameters
    #     ----------
    #     trajectory : str 
    #         Name of the .xtc file.
    #     out_file : str 
    #             Name of the output file.
    #     """
    #     cutoff = float(cutoff)
    #     cutoff =str(cutoff)
    #     if cutoff != 3.0:
    #         a = subprocess.call(['fatslim', 'apl',
    #                              '-c',self.gro,
    #                              '-n',self.headgroups_ndx_file,
    #                              '-t',trajectory,
    #                              '--nthreads',self.thread,
    #                              '--export-apl-raw',out_file,
    #                              '--apl-cutoff',cutoff],
    #                             )
    #     else:
    #         a = subprocess.call(['fatslim', 'apl',
    #                              '-c',self.gro,
    #                              '-n',self.headgroups_ndx_file,
    #                              '-t',trajectory,
    #                              '--nthreads',self.thread,
    #                              '--export-apl-raw',out_file]
    #                             )


    def AreaPerLipid(self,
                     trajectory,
                     cutoff, 
                     out_file):

        """Execute the APL command of fatslim,
        which compute the average area per lipid 
        of lower,upper leaflet and the entire membrane 
        in a ''.xvg file along the trajectory.
    
        Parameters
        ----------
        trajectory : str 
            Name of the .xtc file.
        out_file : str 
                Name of the output file.
        """
        cutoff = float(cutoff)
        cutoff =str(cutoff)
        if cutoff != 3.0:

            a = subprocess.call(['fatslim', 'apl',
                                 '-c',self.gro,
                                 '-n',self.headgroups_ndx_file,
                                 '-t',trajectory,
                                 '--nthreads',self.thread,
                                 '--plot-apl',out_file,
                                 '--apl-cutoff', cutoff],
                                )
        else:
            a = subprocess.call(['fatslim', 'apl',
                                 '-c',self.gro,
                                 '-n',self.headgroups_ndx_file,
                                 '-t',trajectory,
                                 '--nthreads',self.thread,
                                 '--plot-apl',out_file],
                                )


###############################################################################################################################################################    
   
    
class TrajCommandModified:

    def __init__(self,
                 trajectory,
                 topology,
                 headgroups_ndx_file,
                 ):
    
        self.trajectory = trajectory
        self.topology = topology
        self.headgroups_ndx_file = headgroups_ndx_file


    def leaflets(self):
        
        
        # It uses fatslim program to find upper and lower leaflet  and 
        # make an .ndx that cointains both it then split them in specific 
        # index that will be used later
        

        
        if os.path.isfile('index_lower_leaflet_res.ndx'):
            print()
        else:
            # Lower Leaflet index
            make_ndx = tools.Make_ndx(f = self.topology, \
                                      n = 'true_bilayer.ndx',
                                      o = 'index_lower_leaflet_res.ndx', \
                                      input = ('del 1',
                                               'name 0 lower_leaflet', \
                                               'splitres 0', \
                                               'del 0',
                                               'q')
                                     )
            make_ndx.run()


        # Upper leaflet index
        if os.path.isfile('index_upper_leaflet_res.ndx'):
            print()
        else:
            make_ndx = tools.Make_ndx(f = self.topology, \
                                      n = 'true_bilayer.ndx',
                                      o = 'index_upper_leaflet_res.ndx', \
                                      input = ('del 0',
                                               'name 0 upper_leaflet', \
                                               'splitres 0', \
                                               'del 0',
                                               'q')
                                      )
            make_ndx.run()


    def get_lipids_indexes(self,
                           leaflet_ndx):

        """ Method to create, for each single lipid molecule
        constituting upper and lower leaflet, a corresponding
        ''.ndx file , from the main ''.ndx file.

        Parameters
        ----------
        leaflet_ndx : ''.ndx file
                Upper or lower leaflet ndx
        """ 
        

        # Create different txt file as many as the different lipid / res in the ndx file
        
        data = {}
        with open(leaflet_ndx, "r") as f:
            keep_parsing = True
            for line in f:
                if "[" in line:
                    print(line)
                    header = line.strip("\n").strip("[").strip("]").strip(" ").split("_")[2] \
                    + "_" +  line.strip("\n").strip("[").strip("]").strip(" ").split("_")[3]
                    keep_parsing = False
                    header = header.replace('*', '')
                else:
                    if not keep_parsing:
                        data[header] = line
                    else:
                        data[header] += line
                    keep_parsing = True


        # Creates all the different files

        for header in data.keys():
            with open(header +".ndx", "w") as out:
                    out.write("[ %s ]\n%s" % (header, data[header]))


    def get_xvg_lipids(self,
                       dir_name
                       ):
        
        """Method to produce the different xvg
        file using the ndxs produced by the 
        get_lipids_indexes().
        One xvg for each index, to extract the
        X,Y,Z coordinates of each single lipid residue
        Parameter
        --------- 
        dir_name : str
            Name of the folder in which the method operate 
        """

        starting_directory = os.getcwd() 

        dirct = os.getcwd()+'/'+ dir_name # get the path

        os.chdir(dirct)
      
        for i in os.listdir(dirct):

            if i.endswith('.ndx'): #Execute the traj command for each file in the folder
                
                traj = tools.Traj(f = starting_directory+'/'+ self.trajectory,
                                  s = starting_directory+'/'+ self.topology,
                                  ox = i.split('.')[0] + '.xvg',
                                  com = True,
                                  n = dirct + '/' + i)
                traj.run()
        

        # Modifies the xvg to remove the first  25 lines
        # and substitute the first line with the name of
        # the lipid residue and its number
        os.system("for i in *.xvg ;\
                   do var=`echo $i | sed 's/.xvg//'` ; \
                   sed -i '1,25d' $var.xvg ; \
                   perl -pi -e 's/^(.*)/'$var'/ \
                   if $.==1' $var.xvg ; \
                   done"
                  )

        # If in the directory is present the either the word 
        # upper or lower then create and append to the Merged 
        # file each xvg, in order to have a novel xvg file
        # to be plotted 
        if "upper" in dir_name : 
            os.system('printf  "# File created on $(date)\n# Created by $(whoami) \n# For the analysis of membranes\n" \
                   > Merged_coord_upper_leaflet.txt'  )
            os.system("cat *.xvg >> Merged_coord_upper_leaflet.txt")
            shutil.move("Merged_coord_upper_leaflet.txt",starting_directory)
            shutil.rmtree(dirct)
        elif "lower" in dir_name :
            os.system('printf  "# File created on $(date)\n# Created by $(whoami) \n# For the analysis of membranes\n" \
                   > Merged_coord_lower_leaflet.txt'  )
            os.system("cat *.xvg >> Merged_coord_lower_leaflet.txt")
            shutil.move("Merged_coord_lower_leaflet.txt ",starting_directory)
            shutil.rmtree(dirct)


