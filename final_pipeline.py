#!/usr/bin/env python3
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

import time
import errno
import subprocess
import os
import argparse
import shutil
from pathlib import Path
import pipeline_tools
import MDAnalysis as mda
from pathlib import Path


# Import GROMACS-related modules
import gromacs
import gromacs.tools as tools
import gromacs.setup as setup

# set-up GROMACS
gromacs.config.setup()


_start_time = time.time()

def tic():
    global _start_time 
    _start_time = time.time()


def tac():
    t_sec = round(time.time() - _start_time)
    (t_min, t_sec) = divmod(t_sec,60)
    (t_hour,t_min) = divmod(t_min,60) 
    print('Time passed: {}hour:{}min:{}sec'.format(t_hour,t_min,t_sec))



def data_pre_processing(u,
	                trajectory_file,
                        tpr_file,
                        starting_directory,
                        checkpoint_file,
                        skip,
                        begin,
                        end,
                        directory_name
                        ):
    
    """Pre-processing of the file.
    It centers the sistem and depending on the user
    choice skips ps on the trajectory or cuts it in
    specific frame

    Parameters
    ----------
    trajectory_file : str 
            Name of the .xtc file.
    tpr_file : str
            Name of the .tpr file.
    starting_directory : str
            Name of the starting directory
    checkpoint_file : str 
            Name of the .cpt file.
    skip : float
            Number of ps to skip 
    begin : float
            Number wich represent the starting point
    end : float 
          Number wich represent the end point
    """

    # If the root directory Analysis is not present, make it.
    # Needed further for storage purpose. 
    
    

    if not os.path.exists(directory_name):
        try:
            os.mkdir((directory_name))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    else:
        # If the same dir exist remove it and the create it again 
        shutil.rmtree(starting_directory +'/'+ directory_name)
        os.mkdir((directory_name))


    # Create an index of only the lipids of  the system .  
    solv = u.select_atoms(" resname TIP3") + \
           u.select_atoms(" resname HT") + \
           u.select_atoms(" resname HT") + \
           u.select_atoms(" resname HX") + \
           u.select_atoms(" resname OT") + \
           u.select_atoms(" resname OX") + \
           u.select_atoms(" resname LIT")+ \
           u.select_atoms(" resname SOD")+ \
           u.select_atoms(" resname MG") + \
           u.select_atoms(" resname POT")+ \
           u.select_atoms(" resname CAL")+ \
           u.select_atoms(" resname RUB")+ \
           u.select_atoms(" resname CES")+ \
           u.select_atoms(" resname BAR")+ \
           u.select_atoms(" resname ZN") + \
           u.select_atoms(" resname CAD")+ \
           u.select_atoms(" resname CLA")

    system = u.select_atoms("all") 

    memb = system.difference(solv) # (everything but protein )- (all solvent + solute) = only membrane
                                   # and in case of presence of protein , the protein 

    memb.write(filename= "index_no_solvent.ndx", name = "All_but_solvent")

    # If the user choose to skip ps on the trajectory or prefer to analyze 
    # its entirety ,cutting or not via --begin or --end
    # Create the .gro associated with the trajectory with only the last frame in it 
    processed_traj = 'processed_traj.xtc'
    gro_file = 'last_frame.gro'
    
    if os.path.isfile(processed_traj) and os.path.isfile(gro_file):
        print() # path exists
    else:
        if skip == None and (begin == None and end == None):
        # There is no ps to skip
            trjconv = tools.Trjconv(f = trajectory_file, \
		                            s = tpr_file, \
		                            ur = 'compact', \
		                            n = 'index_no_solvent.ndx', \
		                            pbc = 'mol', \
		                            o = processed_traj, \
		                            )
            trjconv.run()

        elif skip == None and (begin is not None and end is not None): 
        # Only cutting the traj and creating a new processed .xtc file
            trjconv = tools.Trjconv(f = trajectory_file, \
		                            s = tpr_file, \
		                            ur = 'compact', \
		                            n = 'index_no_solvent.ndx', \
		                            pbc = 'mol', \
		                            o = processed_traj, \
		                            b = begin, \
		                            e = end 
		                            )
            trjconv.run()

        elif skip is not None and (begin is not None and end is not None): 
        #skipping and cutting the traj and again create a new processed .xtc file
            trjconv = tools.Trjconv(f = trajectory_file, \
		                            s = tpr_file, \
		                            ur = 'compact', \
		                            n = 'index_no_solvent.ndx', \
		                            pbc = 'mol', \
		                            o = processed_traj, \
		                            b = begin, \
		                            e = end ,\
		                            skip = skip
		                            )
            trjconv.run()

        elif skip is not None and (begin == None and end == None): 

        # The user is only skipping ps in the traj
            trjconv = tools.Trjconv(f = trajectory_file, \
		                            s = tpr_file, \
		                            ur = 'compact', \
		                            n = 'index_no_solvent.ndx', \
		                            pbc = 'mol', \
		                            o = processed_traj, \
		                            skip = skip
		                            )
            trjconv.run()

        trjconv = tools.Trjconv(f = checkpoint_file, \
                                s = tpr_file, \
                                ur = 'compact', \
                                n = 'index_no_solvent.ndx', \
                                pbc = 'mol', \
                                o = gro_file, \
                                )
        trjconv.run()

    
    
    

  
   
        

def module_fatslim(processed_traj_file,
                   gro_file,
                   index_headgroups,
                   directory_name,
                   apl_cutoff,
                   tck_cutoff
                   ):
    
    print('---------------------------------------------------\n')
    
    """Fuction consisting of fatslim analysis 
    1) Thickness (raw +xvg)
    2) APL (raw+xvg)

    Parameters
    ----------
    processed_traj_file : str
                            Name of the processed .xtc file.
    gro_file : str 
            Name of the last producte .gro file.
    """


    
    # Creates different folder in which the output is stored
    
    analysis_fatlism = directory_name + '/Fatslim/'

    if not os.path.exists(analysis_fatlism):
        try:
            os.mkdir((analysis_fatlism))
            os.mkdir((analysis_fatlism +'/fatslim_apl'))
            os.mkdir((analysis_fatlism + '/fatslim_thickness'))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    

  
    # Begins of the fatslim analysis 
    


    fatslim = pipeline_tools.FatslimCommands(gro_file,
                                             index_headgroups,
                                             '2'
                                            )
    

    fatslim.thickness(trajectory = processed_traj_file,
                      cutoff = tck_cutoff,
    	              out_file = analysis_fatlism + \
                                  '/thickness.xvg'

                     )
     
    # fatslim.raw_thickness(trajectory = processed_traj_file,
    #                       cutoff = tck_cutoff,
    # 	                  out_file = analysis_fatlism + \
    #                              '/fatslim_thickness/raw_thickness.csv'
    #                       )

    fatslim.AreaPerLipid(trajectory = processed_traj_file,
                         cutoff = apl_cutoff,
    	                 out_file = analysis_fatlism + \
                                   '/apl.xvg')

    # fatslim.raw_AreaPerLipid(trajectory = processed_traj_file,
    #                          cutoff = apl_cutoff,
    # 	                     out_file = analysis_fatlism +
    #                                     '/fatslim_apl/raw_apl')

    

    
   
    



def module_densmap(processed_traj_file,
                   directory_name):
    
    """ Compute the 2D maps of the system exploiting densmap command
    of gromacs

    Parameters
    ----------
    processed_traj_file : str
                                Name of the processed .xtc file.
    """


    folder = directory_name +'/2D_maps/'

    if not os.path.exists(folder):
        try:
            os.mkdir((folder))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


    

    # Lower Leaflet , input == 0
    densmap = tools.G_densmap(f = processed_traj_file,\
                              n = 'true_bilayer.ndx',\
                              o = folder +'lower.xpm',\
                              input = ('0')
                             )
    densmap.run()

    # Upper Leaflet , input == 1
    densmap = tools.G_densmap(f = processed_traj_file,\
                              n = 'true_bilayer.ndx',\
                              o = folder + 'upper.xpm',\
                              input = ('1')
                             )
    densmap.run()
    
    #conversion of the xpm file into eps using xpm2ps
    xpm2ps = tools.Xpm2ps(f = folder + 'lower.xpm', \
                          o = folder + 'true_lower.eps', \
                          rainbow = 'red'
                          )
    xpm2ps.run()
    
    xpm2ps = tools.Xpm2ps(f = folder + 'upper.xpm', \
                          o = folder + 'true_upper.eps', \
                          rainbow = 'red'
                          )
    xpm2ps.run()
    
  



def module_movements(processed_traj_file,
                     gro_file,
                     index_headgroups,
                     directory_name
	                 ):
     
    print('---------------------------------------------------\n')
    
    # Creates different folder in which the output is stored

    analysis_traj = directory_name + "/Movements/"

    if not os.path.exists(analysis_traj):
        try:
            os.mkdir((analysis_traj))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    """From the imported module uses these function to extract the trajectories
    associated with each lipids of both upper and lower leaflet

    Parameters
    ----------
    processed_traj_file : str
                            Name of the processed .xtc file. 
    gro_file : str 
            Name of the last producte .gro file.
    tpr_file : str
            Name of the .tpr file.
    index_headgroups : str
            Name of the index file
    """

    starting_directory = os.getcwd()

    index_lower_leaflet_res = starting_directory + '/' + 'index_lower_leaflet_res.ndx'
    index_upper_leaflet_res = starting_directory + '/' + 'index_upper_leaflet_res.ndx'

    traj_cmd_modified = pipeline_tools.TrajCommandModified(processed_traj_file,
                                                           gro_file,
                                                           index_headgroups,
                                                           "4"
                                                           )
    traj_cmd_modified.leaflets() 
    os.mkdir('upper_leaflet') # create the folder
    os.chdir('upper_leaflet') # enter the folder
    traj_cmd_modified.get_lipids_indexes(index_upper_leaflet_res) # upper leaflet indexes
    os.chdir(starting_directory)
    traj_cmd_modified.get_xvg_lipids("upper_leaflet") 

    os.chdir(starting_directory) # return to the home folder
    
    os.mkdir('lower_leaflet') # create the folder
    os.chdir('lower_leaflet') # enter the folder
    traj_cmd_modified.get_lipids_indexes(index_lower_leaflet_res) # lower leaflet indexes
    os.chdir(starting_directory)
    traj_cmd_modified.get_xvg_lipids("lower_leaflet") 

    os.chdir(starting_directory)
    
                  
  
    a = subprocess.call(['mv',
                         'Merged_coord_lower_leaflet.txt',
                          analysis_traj],
                        )
    
    a = subprocess.call(['mv',
                         'Merged_coord_upper_leaflet.txt',
                          analysis_traj],
                        )




def module_order_parameter(processed_traj_file,
	                       gro_file,
                           directory_name
	                       ):

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
    Execute the python script for the order parameter for 
    all the .def files
    To work it uses the .def files which contains all the 
    information in the format :
      
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
    
    Parameters
    ----------
    processed_traj_file : str
                    Name of the processed .xtc file. 
    gro_file : str 
            Name of the last producte .gro file.

    """
    

    for i in pipeline_tools.d:

        ordPars = pipeline_tools.parse_op_input(pipeline_tools.d[i])

        pipeline_tools.read_trajs_calc_OPs(ordPars,
                                           gro_file,
                                           [processed_traj_file]
                                           )

        for op in ordPars.values():
            (op.avg, op.std, op.stem) = op.get_avg_std_stem_OP
             
        # writes the output file in .dat file 
        try:
            with open('Order_Parameter_'+ i +'.dat',"w") as f:
                f.write("# OP_name    resname    atom1    atom2    OP_mean   OP_stddev  OP_stem\n")
                f.write("#--------------------------------------------------------------------\n")
                for op in ordPars.values():
                    f.write( "{:20s} {:7s} {:5s} {:5s} {: 2.5f} {: 2.5f} {: 2.5f} \n".format(
                             op.name, op.resname, op.atAname, op.atBname,
                             op.avg, op.std, op.stem)
                           )
        except:
            print("Problem writing the main outputfile")
        

        with open('Order_Parameter_'+ i +'.dat','r') as g :
            if "nan" in g.read():
                subprocess.call(['rm','Order_Parameter_'+ i +'.dat'])
            else:
                a = subprocess.call(['mv','Order_Parameter_'+ i +'.dat',
                              order_folder])               

         
def cleaning_all():

# Clean all the temporary files which were created 


    print('---------------------------------------------------\n')

    a = subprocess.call(['rm', 'index_no_solvent.ndx'],)
    a = subprocess.call(['rm', 'index_headgroups.ndx'],)
    a = subprocess.call(['rm', 'index_lower_leaflet_res.ndx'],)
    a = subprocess.call(['rm', 'index_upper_leaflet_res.ndx'],)
    a = subprocess.call(['rm bilayer_*'], shell = True)
    a = subprocess.call(['rm', 'true_bilayer.ndx'],)
    a = subprocess.call(['rm \#_*'], shell = True)






def main():

    
    # Main function 
    
    description= "Pipeline that calculates four different parameter for a MD simulation \
    of a lipid membrane coming from gromacs"

    usage = "python3 pipeline.py -t ["".xtc/trr] \
                                 -f ["".tpr] \
                                 -g [*.cpt]\
                                 -s [skip ps] \
                                 -b [begin]\
                                 -e[end]\
                                 -all [True/False]\
                                 -clean  \
                                 -fatslim\
                                 -2d\
                                 -mov\
                                 -ordpar\
                                 -d [directory name] "
             

    parser = argparse.ArgumentParser(description='', usage=usage)

    parser.add_argument('-t',
                        '--trajectory',
                        dest='trajectory',
                        type=str,
                        required=True,
                        metavar='',
                        help='Trajectory files [<.xtc/.trr/>] ',
                        )

    parser.add_argument('-f',
                        '--tpr file',
                        dest='topology',
                        type=str,
                        required=True,
                        metavar='',
                        help=' molecular topology and all the simulation parameters ',
                        )

    parser.add_argument('-g',
                        '--checkpoint file',
                        dest='checkpoint',
                        type=str,
                        required=True,
                        metavar='',
                        help='<.cpt> file',
                        )

    parser.add_argument('-d',
                        '--dir_name',
                        default = 'Analysis',
                        dest='directory',
                        type=str,
                        required=False,
                        metavar='',
                        help='The directory name for the analysis',
                        )

    parser.add_argument('-s',
                        '--skip',
                        default = None,
                        dest='skip',
                        type=str,
                        required = False,
                        help='How much ps are intended to skip in the trajectory',
                        metavar='',
                        )

    parser.add_argument('-b',
                        '--begin',
                        default = None,
                        dest='begin',
                        type=str,
                        required = False,
                        help='Customization of the trajectories',
                        metavar='',
                        )

    parser.add_argument('-e',
                        '--end',
                        default = None,
                        dest='end',
                        type=str,
                        required = False,
                        help='Customizatione of the trajectories',
                        metavar='',
                        )

    parser.add_argument('-c',
                        '--clean',
                        action='store_true',
                        dest='clean',
                        required = False,
                        help='Clean all the temporary .ndx files '
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
                        help='Cutoff distance (in nm) used to approximate planar\
                        region (default: 3.0)'
                        )

    parser.add_argument('-thickness_cutoff',
                        '--thickness-cutoff',
                        default= 6.0,
                        dest='tck_cutoff',
                        required = False,
                        help='Cutoff distance (in nm) used to identify inter-leaflet\
                        neighbors (default: 6.0)'
                        )

    parser.add_argument('-2d',
                        '--2d_maps',
                        action='store_true',
                        dest='mod_maps', 
                        required = False,
                        help='Execute only the density maps module of the pipeline'
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
                        help='Execute only the order parameter module of the pipeline'
                        )
    
    parser.add_argument('-prot',
                        '--protein',
                        action='store_true',
                        dest='prot',
                        required = False,
                        help='Specify if a protein is embedded in the membrane'
                       )

    args = parser.parse_args()

    if args is None:
        print("Are you sure you followed the instruction? Please follow them and retry")
        pass
    else:
        
        trajectory_file = args.trajectory
        tpr_file = args.topology
        checkpoint_file = args.checkpoint
        directory_name = args.directory
        begin = args.begin  
        end = args.end  
        skip = args.skip
        starting_directory = os.getcwd()
        cleaning = args.clean
        apl_cutoff = args.apl_cutoff
        tck_cutoff = args.tck_cutoff
        

        tic()
        
        u = mda.Universe(tpr_file,trajectory_file)

        data_pre_processing(u,
        	            trajectory_file,
                            tpr_file,
                            starting_directory,
                            checkpoint_file,
                            skip,
                            begin,
                            end,
                            directory_name
                            )

        
        # If the trajectory was skipped then the processed xtc file exist,
        # therefore it is used by the different modules, if not used the 
        # xtc file with all the trajectory 

        # P is for all phopholipids 
        # resname CHOL to check for cholesterol
        # resname ERG to check for ergosterol

        

        processed_traj_file = 'processed_traj.xtc'
        gro_file = 'last_frame.gro'
        
        # means a protein is embedded in the bilayer
        if args.prot :

        	# index for the fatslim analysis
            g = u.select_atoms("name P") + \
                u.select_atoms("resname ERG and name O3") +  \
                u.select_atoms("resname CHL1 and name O3")
            p = u.select_atoms("protein and not name H*")
            g.write(filename= "index_headgroups.ndx",name="headgroups",mode ="a")
            p.write(filename= "index_headgroups.ndx",name="protein",mode = "a")


            # Sometimes Fatslim has problem with ndx written by MDAnalysis
            # so we use make_ndx to reconvert it, we use MDAnalysis because
            # we can have more control over the selections of headgroups 
            make_ndx = tools.Make_ndx(f = tpr_file, \
            	                      n = "index_headgroups.ndx",\
                                      o = "index_headgroups_corrected.ndx", \
                                      input = ('q')
                                      )
            make_ndx.run()
        else:
            g = u.select_atoms("name P") + \
                u.select_atoms("resname ERG and name O3") +  \
                u.select_atoms("resname CHL1 and name O3")
            g.write(filename= "index_headgroups.ndx",name="headgroups")


            # same applies here
            make_ndx = tools.Make_ndx(f = tpr_file, \
            	                      n = "index_headgroups.ndx",\
                                      o = "index_headgroups_corrected.ndx", \
                                      input = ('q')
                                      )

            make_ndx.run()

        
        # index file of headgroups
        index_headgroups = 'index_headgroups_corrected.ndx'

        # Call Fatslim class
        fatslim = pipeline_tools.FatslimCommands(gro_file,
                                                 index_headgroups,
                                                 '2'
                                                )

        # Create an index file in which the different layers of the membrane
        # are listed 

        home_dir = os.getcwd()

        fatslim.membranes(out_file = home_dir +"/bilayer.ndx")
        os.rename('bilayer_0000.ndx', 'true_bilayer.ndx')
        a = subprocess.call(['rm bilayer_*.ndx'], shell = True)
        
        

        # means that the user selected the -all flag 
        # and he/she wants all the pipeline to be
        # executed

        if args.all_modules :

            module_fatslim(processed_traj_file,
                           gro_file,
                           index_headgroups,
                           directory_name,
                           apl_cutoff,
                           tck_cutoff
                           )

            module_densmap(processed_traj_file,
                           directory_name)

            module_movements(processed_traj_file,
                             gro_file,
                             index_headgroups,
                             directory_name
                            )

            module_order_parameter(processed_traj_file,
                                   gro_file,
                                   directory_name
                                   )
            
        else:


            if args.mod_fat : 
            	module_fatslim(processed_traj_file,
                               gro_file,
                               index_headgroups,
                               directory_name,
                               apl_cutoff,
                               tck_cutoff
                               )

            if args.mod_maps :
            	module_densmap(processed_traj_file,
                               directory_name)

            if args.mod_mov :
            	module_movements(processed_traj_file,
                                 gro_file,
                                 index_headgroups,
                                 directory_name
                                 )
            if args.mod_ord :
            	module_order_parameter(processed_traj_file,
                                       gro_file,
                                       directory_name
                                       )

            else:
            	print()

         
        if args.clean :
            cleaning_all()
            tac()
        else:
            tac()
            print()

if __name__ == "__main__":
    main()
