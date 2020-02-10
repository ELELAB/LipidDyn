#!/usr/bin/env python3
"""
#****************************************************************************
#* Author : Simone Scrima
#* Group : Computational Biology Laboratory ; Danish Cancer Society Research 
#*         Center ; CPH ; Denmark 
#*   PLATFORM : Linux-64
#* Year : 2020
#****************************************************************************

The pipeline exploits a free-available softwares for the automatic 
computation of different biophysical parameters for lipids membranes 
(i.e. GROMACS, FATSLiM).

It allows to compute : 
1) 2D density maps 
2) Thickness and Area Per lipid 
3) Order Parameter 
4) "Movements" ( .txt file containing the X Y Z coordinates 
    for the invidual lipids molecule)
 
"""
# import packages
import time
import errno
import subprocess
import os
import argparse
import shutil
from pathlib import Path
import glob

#Import the classes and functions from the pipeline_tools
import pipeline_tools

#Import MDanalysis
import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder

# Import GROMACS-related modules
import gromacs
import gromacs.tools as tools
import gromacs.setup as setup

# set-up GROMACS
gromacs.config.setup()

# Import Loguru for log files
from loguru import logger
logger.debug("Debugger,let's hope!")
logger.add("Log_criteria.log")





_start_time = time.time()

def tic():
    global _start_time 
    _start_time = time.time()


def tac():
    t_sec = round(time.time() - _start_time)
    (t_min, t_sec) = divmod(t_sec,60)
    (t_hour,t_min) = divmod(t_min,60) 
    print('Time passed: {}hour:{}min:{}sec'.format(t_hour,t_min,t_sec))

    

def module_fatslim(trajectory_file,
                   topology_file,
                   index_headgroups,
                   directory_name,
                   apl_cutoff,
                   tck_cutoff,
                   raw,
                   ncore
                   ):
    
    print('---------------------------------------------------\n')
    
    """Fuction consisting of fatslim analysis 
    1) Thickness (raw +xvg)
    2) APL (raw+xvg)

    Parameters
    ----------
    trajectory_file : str
            Name of the processed .xtc file.
    topology_file : str 
            Name of the last producte .gro file.
    index_headgroups : str
            Filename of the index containing all the headgroups
    directory_name : str
            Name of the output folder
    apl_cutoff : float
            Number representing the cut-off for the APL
    tck_cutoff : float
            Number representing the cut-off for the thickness
    raw : True/False
            Value to include or not the raw output
    ncore : int
            Number of cores for parallelization             
    """


    
    # Creates different folder in which the output is stored
    
    analysis_fatlism = os.path.abspath(directory_name + '/Fatslim/')

    if not os.path.exists(analysis_fatlism):
        try:
            os.mkdir((analysis_fatlism))
            os.mkdir((analysis_fatlism +'/fatslim_apl'))
            os.mkdir((analysis_fatlism + '/fatslim_thickness'))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    

  
    # Begins of the fatslim analysis 
    
    fatslim = pipeline_tools.FatslimCommands(trajectory_file,
                                             topology_file,
                                             index_headgroups,
                                             ncore,
                                             apl_cutoff,
                                             tck_cutoff
                                             )
    

    fatslim.thickness(out_file = analysis_fatlism +'/thickness.xvg'
                     )
     
    
    fatslim.AreaPerLipid(out_file = analysis_fatlism + '/apl.xvg')

    if raw : # if the user required the raw data
        fatslim.raw_AreaPerLipid(out_file = analysis_fatlism +
                                 '/fatslim_apl/raw_apl.csv')      
        fatslim.raw_thickness(out_file = analysis_fatlism + \
                              '/fatslim_thickness/raw_thickness.csv'
                              )
    else: 
        pass

   

def module_densmap(trajectory_file,
                   directory_name):
    
    """ Compute the 2D maps of the system exploiting densmap 
    command of gromacs

    Parameters
    ----------
    trajectory_file : str
            Name of the processed .xtc file.
    directory_name : str
            Name of the output folder
    """


    folder = os.path.abspath(directory_name +'/2D_maps/')

    if not os.path.exists(folder):
        try:
            os.mkdir((folder))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    

    # Lower Leaflet , input == 0
    densmap = tools.G_densmap(f = trajectory_file,\
                              n = 'true_bilayer.ndx',\
                              o = folder +'/lower.xpm',\
                              input = ('0')
                             )
    densmap.run()

    # Upper Leaflet , input == 1
    densmap = tools.G_densmap(f = trajectory_file,\
                              n = 'true_bilayer.ndx',\
                              o = folder + '/upper.xpm',\
                              input = ('1')
                             )
    densmap.run()
    
    # conversion of the xpm file into eps using xpm2ps
    xpm2ps = tools.Xpm2ps(f = folder + '/lower.xpm', \
                          o = folder + '/true_lower.eps', \
                          rainbow = 'red'
                          )
    xpm2ps.run()
    
    xpm2ps = tools.Xpm2ps(f = folder + '/upper.xpm', \
                          o = folder + '/true_upper.eps', \
                          rainbow = 'red'
                          )
    xpm2ps.run()
    
  



def module_movements(trajectory_file,
                     topology_file,
                     index_headgroups,
                     directory_name,
                     ncore
                     ):
     
    print('---------------------------------------------------\n')
    
    # Creates different folder in which the output is stored

    """From the imported module uses these function to extract the 
    trajectories associated with each lipids of both upper and lower 
    leaflet

    Parameters
    ----------
    trajectory_file : str
            Name of the processed .xtc file. 
    topology_file : str 
            Name of the last producte .gro file.
    index_headgroups : str
            Filename of the index containing all the headgroups
    directory_name : str
            Name of the output folder
    ncore : int
            Number of cores for parallelization 
    """

    analysis_traj = os.path.abspath(directory_name + "/Movements/")

    if not os.path.exists(analysis_traj):
        try:
            os.mkdir((analysis_traj))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    

    starting_directory = os.path.abspath(os.getcwd())

    index_lower_leaflet_res = starting_directory + \
                              '/index_lower_leaflet_res.ndx'
    index_upper_leaflet_res = starting_directory + \
                              '/index_upper_leaflet_res.ndx'

    traj_cmd_modified = pipeline_tools.TrajCommandModified(trajectory_file,
                                                           topology_file,
                                                           index_headgroups,
                                                           ncore
                                                           )

    traj_cmd_modified.leaflets() # 
    
    if not os.path.exists("upper_leaflet"):
        try:
            os.mkdir(("upper_leaflet"))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    
    # upper leaflet indexes
    os.chdir("upper_leaflet")
    traj_cmd_modified.get_lipids_indexes(index_upper_leaflet_res) 
    traj_cmd_modified.get_xvg_lipids()  
    
    # Append to the Merged file each .xvg to have a novel xvg file to be 
    # plotted 
    # (upper leaflet)
    os.system('printf  "# File created on $(date)\n\
                        # Created by $(whoami) \n\
                        # For the analysis of membranes\n" \
                        > Merged_coord_upper_leaflet.txt'  )
    os.system("cat *.xvg >> Merged_coord_upper_leaflet.txt")
    shutil.move("Merged_coord_upper_leaflet.txt", analysis_traj)
    os.chdir(starting_directory)

    if not os.path.exists("lower_leaflet"):
        try:
            os.mkdir(("lower_leaflet"))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    
    # lower leaflet indexes
    os.chdir("lower_leaflet")
    traj_cmd_modified.get_lipids_indexes(index_lower_leaflet_res) 
    traj_cmd_modified.get_xvg_lipids() 
    
    # Append to the Merged file each .xvg to have a novel xvg file to be 
    # plotted  
    #(lower leaflet)
    os.system('printf  "# File created on $(date)\n\
                        # Created by $(whoami) \n\
                        # For the analysis of membranes\n" \
                        > Merged_coord_lower_leaflet.txt'  )
    os.system("cat *.xvg >> Merged_coord_lower_leaflet.txt")
    shutil.move("Merged_coord_lower_leaflet.txt",analysis_traj)
    os.chdir(starting_directory)
    
    
    for dir_name in ["upper_leaflet","lower_leaflet"] :
        shutil.rmtree(dir_name) # remove the folders with all .xvgs
        


def module_order_parameter(trajectory_file,
                           topology_file,
                           directory_name
                           ):
    
    """

    Parameters
    ----------
    trajectory_file : str
            Name of the processed .xtc file. 
    topology_file : str 
            Name of the last producte .gro file.
    directory_name : str
            Name of the output folder
    """

    print('---------------------------------------------------\n')
    
    # Creates different folder in which the output is stored

    order_folder = directory_name + '/Order_Parameter/'

    if not os.path.exists(order_folder):
        try:
            os.mkdir((order_folder))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    
    """ 
    Execute the python script for the order parameter 
    To work it uses the information in the format :
      
    OP_name1 Residue Carbon_name Hydrogen_name
    OP_name2 Residue Carbon_name Hydrogen_name
    ------------------------------------------
    Example (CHARMM36):
 
    beta1 POPC C12 H12A
    beta2 POPC C12 H12B
    alpha1 POPC C11 H11A
    alpha2 POPC C11 H11B
    g3_1 POPC C1 HA
    g3_2 POPC C1 HB
    g2_1 POPC C2 HS
    g1_1 POPC C3 HX
    g1_2 POPC C3 HY 
    
    present in  pipeline_tools.py
    
    Parameters
    ----------
    trajectory_file : str
                    Processed .xtc file. 
    topology_file : str 
            Name of the x.gro file.

    """
    

    for def_file in glob.iglob('definitions_files/*.def'):
        
        lipid_name = os.path.splitext(def_file)[0].split("/")[1]
        ordPars = pipeline_tools.parse_op_input(os.path.abspath(def_file))

        pipeline_tools.read_trajs_calc_OPs(ordPars,
                                           topology_file,
                                           [trajectory_file]
                                           )

        for op in ordPars.values():
            (op.avg, op.std, op.stem) = op.get_avg_std_stem_OP
             
        # writes the output file in .dat file 
        try:
            with open('Order_Parameter_'+ lipid_name +'.dat',"w") as f:
                for op in ordPars.values():
                    f.write( "{:20s} {:7s} {:5s} {:5s} {: 2.5f} {: 2.5f} {: 2.5f} \n".format(
                              op.name, op.resname, op.atAname,
                              op.atBname,op.avg, op.std, op.stem)
                           )
        except:
            print("Problem writing the main outputfile")
        

        with open('Order_Parameter_'+ lipid_name +'.dat','r') as g :
            if "nan" in g.read():
                os.remove('Order_Parameter_'+ lipid_name +'.dat')
            else:
                shutil.move('Order_Parameter_'+ lipid_name +'.dat',
                             order_folder) 

         
def cleaning_all():

    # Clean all the temporary files which were created 


    print('---------------------------------------------------\n')
    for filename in glob.glob("bilayer_*.ndx"):
            os.remove(filename)
    os.remove('system_no_solvent.ndx')
    os.remove('index_headgroups.ndx')
    os.remove('index_lower_leaflet_res.ndx')
    os.remove('index_upper_leaflet_res.ndx')
    os.remove('true_bilayer.ndx')
    os.remove('index_headgroups_corrected.ndx')
    





def main():

    
    # Main function 
    
    description= "LypidDyn pipeline for calculating different parameters \
                  of a lipid membrane Molecular Dynamics simulation"

    usage = "python3 pipeline.py \
             -t ["".xtc/trr] \
             -f ["".gro/tpr] \
             -all [True/False] \
             -c -fatslim \
             -2d \
             -mov \
             -ordpar \
             -d [directory name]"
             

    parser = argparse.ArgumentParser(description=description, usage=usage)

    parser.add_argument('-t',
                        '--trajectory',
                        dest='trajectory',
                        type=str,
                        required=True,
                        metavar='',
                        help='Trajectory files [<.xtc/.trr/>] ',
                        )

    parser.add_argument('-f',
                        '--topology_file',
                        dest='topology',
                        type=str,
                        required=True,
                        metavar='',
                        help=' Molecular topology file [.gro/.tpr]',
                        )


    parser.add_argument('-d',
                        '--dir_name',
                        default = 'Analysis',
                        dest='directory',
                        type=str,
                        required=False,
                        metavar='',
                        help='The directory name for outputs of the pipeline',
                        )

    parser.add_argument('-c',
                        '--clean',
                        action='store_true',
                        dest='clean',
                        required = False,
                        help='Clean all the temporary files '
                        )

    parser.add_argument('-all',
                        '--all',
                        dest='all_modules',
                        action='store_true',
                        required = False,
                        help='Execute all module of the pipeline'
                        )

    parser.add_argument('-fatslim',
                        '--fatslim',
                        action='store_true',
                        dest='mod_fat',
                        required = False,
                        help='Execute only the fatslim module of the pipeline'
                        )

    parser.add_argument('-apl_cutoff',
                        '--apl-cutoff',
                        default= 3.0,
                        dest='apl_cutoff',
                        required = False,
                        help='Cutoff distance (in nm) used to approximate \
                              planar region (default: 3.0)'
                        )

    parser.add_argument('-thickness_cutoff',
                        '--thickness-cutoff',
                        default= 6.0,
                        dest='tck_cutoff',
                        required = False,
                        help='Cutoff distance (in nm) used to identify \
                              inter-leaflet neighbors (default: 6.0)'
                        )

    parser.add_argument('-2d',
                        '--2d_maps',
                        action='store_true',
                        dest='mod_maps', 
                        required = False,
                        help='Execute only the density maps module of the \
                              pipeline'
                                                )

    parser.add_argument('-mov',
                        '--movements',
                        action='store_true',
                        dest='mod_mov',
                        required = False,
                        help='Execute only the movement module of the pipeline'
                       )

    parser.add_argument('-ordpar',
                        '--order_parameter',
                        action='store_true',
                        dest='mod_ord',
                        required = False,
                        help='Execute only the order parameter module of \
                              the pipeline'
                        )
    
    parser.add_argument('-prot',
                        '--protein',
                        action='store_true',
                        dest='prot',
                        required = False,
                        help='Specify if a protein is embedded \
                              in the membrane'
                       )

    # parser.add_argument('-r',
    #                    '--resnam',
    #                     dest='resname',
    #                     required = True,
    #                     nargs= '+'
    #                     help='Lipid residues that compose the \
    #                           membrane of interest'
    #                   )

    parser.add_argument('-ncore',
                        '--mpi',
                        nargs ='?',
                        default = 2,
                        dest='mpi',
                        help='Specify on how many cores use for \
                              parallelization'
                       )

    parser.add_argument('-raw',
                        '--raw_data',
                        action='store_true',
                        dest='raw',
                        required = False,
                        help='Get the raw value of thickness and apl for \
                              each frame'
                       )



    args = parser.parse_args()

    if args is None:
        print("Are you sure you followed the instruction? \
               Please follow them and retry")
        quit()
    else:
        
        trajectory_file = os.path.abspath(args.trajectory)
        topology_file = os.path.abspath(args.topology)
        directory_name = args.directory
        #resnames = args.resname
        ncore = args.mpi
        cleaning = args.clean
        apl_cutoff = float(args.apl_cutoff)
        tck_cutoff = float(args.tck_cutoff)
        raw = args.raw

        tic()
        starting_directory = os.getcwd()

        # setup of the MDA universe
        u = mda.Universe(topology_file,trajectory_file)
        
        
        if not os.path.exists(directory_name):
            try:
                os.mkdir((directory_name))
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

        

        
 
        # means a protein is embedded in the bilayer
        if args.prot :

            # index for the fatslim analysis
            # P is for all phopholipids 
            # resname CHOL to check for cholesterol
            # resname ERG to check for ergosterol
            g = u.select_atoms("name P") + \
                u.select_atoms("resname ERG and name O3") +  \
                u.select_atoms("resname CHL1 and name O3")
            p = u.select_atoms("protein and not name H*")

            with mda.selections.gromacs.SelectionWriter(
                                            'index_headgroups_corrected.ndx',
                                             mode='a') as ndx:
                ndx.write(g, name="headgroups", frame=0)
                ndx.write(p, name="protein", frame=0)
    
        else:
            
            # index for the fatslim analysis
            # P is for all phopholipids 
            # resname CHOL to check for cholesterol
            # resname ERG to check for ergosterol
            g = u.select_atoms("name P") + \
                u.select_atoms("resname ERG and name O3") +  \
                u.select_atoms("resname CHL1 and name O3")
            
            with mda.selections.gromacs.SelectionWriter(
                                            'index_headgroups_corrected.ndx',
                                             mode='w') as ndx:
                ndx.write(g, name="headgroups", frame=0)

        # Find leaflets
        L =  LeafletFinder(u,g)
        upper_leaflet = L.groups(0)
        lower_leaflet = L.groups(1)

        # index file of headgroups
        index_headgroups = starting_directory + \
                           '/index_headgroups_corrected.ndx'

        # Initiate Fatslim class
        fatslim = pipeline_tools.FatslimCommands(trajectory_file,
                                                 topology_file,
                                                 index_headgroups,
                                                 ncore,
                                                 apl_cutoff,
                                                 tck_cutoff
                                                )

        # Create an index file in which the different layers of the membrane
        # are listed 

        with mda.selections.gromacs.SelectionWriter('true_bilayer',
                                                     mode='w') as ndx:
                ndx.write(upper_leaflet, name="upper_leaflet")
                ndx.write(lower_leaflet, name="lower_leaflet")
        
        
        

        # The user selected the -all flag for all the analysis

        if args.all_modules :
            
            print("Starting now with the calculation, please stand by")
            print(u"\u2622")

            module_fatslim(trajectory_file,
                           topology_file,
                           index_headgroups,
                           directory_name,
                           apl_cutoff,
                           tck_cutoff,
                           raw,
                           ncore
                           )

            module_densmap(trajectory_file,
                           directory_name)

            module_movements(trajectory_file,
                             topology_file,
                             index_headgroups,
                             directory_name,
                             ncore
                            )

            module_order_parameter(trajectory_file,
                                   topology_file,
                                   directory_name
                                   )
            
        else:

            print("Starting now with the calculation, please stand by")
            print(u"\u2622")

            if args.mod_fat : 
                module_fatslim(trajectory_file,
                               topology_file,
                               index_headgroups,
                               directory_name,
                               apl_cutoff,
                               tck_cutoff,
                               raw,
                               ncore
                               )

            if args.mod_maps :
                module_densmap(trajectory_file,
                               directory_name)

            if args.mod_mov :
                module_movements(trajectory_file,
                                 topology_file,
                                 index_headgroups,
                                 directory_name,
                                 ncore)
            if args.mod_ord :
                module_order_parameter(trajectory_file,
                                       topology_file,
                                       directory_name
                                       )

    

         
        if args.clean :
            cleaning_all()
            tac()
        else:
            tac()

if __name__ == "__main__":
    main()