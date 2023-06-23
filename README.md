
Cancer Structural Biology, Danish Cancer Society Research Center, 2100, Copenhagen, Denmark 

Cancer Systems Biology, Health and Technology Department, Section for Bioinformatics, 2800, Lyngby, Denmark


# LipidDyn: A computational microscope to scrutinize membrane properties at the organelle-level

If you use the tool, please cite:<br/>

Unravelling membrane properties at the organelle-level with LipidDyn

Simone Scrima, Matteo Tiberti, Alessia Campo, Elisabeth Corcelle-Termeau, Delphine Judith, Mads Moller Foged, Knut Kristoffer Bungaard Clemmens,  Sharon A Tooze, Marja Jaattela, Kenji Maeda, Matteo Lambrughi*, Elena Papaleo*

Comput. Struct. Biotechnol. J. 2022 Jun 30;20:3604-3614. doi: 10.1016/j.csbj.2022.06.054 PMID: 35860415 PMCID: PMC9283888 <br/>

contacts: tiberti-at-cancer.dk, matl-at-cancer.dk, elenap-at-cancer.dk, elpap-at-dtu.dk

We have built an automated computational pipeline, implemented in Python 3.8, LipidDyn,
for the accurate analysis of the structural properties and dynamics of lipid bilayers
simulations that can be used to validate the ensembles against experimental data.<br/>
<br/>
The framework is divided into different modules that can be run separately, 
depending on the user choice.<br/> Each performs a different set of analysis on both full-atom and coarse-grained systems.<br/>
Moreover, the pipeline is able to account for embedded proteins into the membrane. <br/>
<br/>
At the present state is designed to work with GROMACS files.
In order to be installed and work, LipidDyn requires Python version >= 3.8. 
The workflow was successfully installed and tested using Ubuntu v18.04 Bionic Beaver -XXX and MacOs 10.15 MacOs 12.1 Monterey.

## Set of Analysis 

### Thickness and Area per Lipid

The membrane thickness and area per lipid are two parameters used to 
validate or compare molecular dynamics simulations of lipids [1] <br/>
In this step we employed FatSlim software that relies on the calculation of local normals
for the estimation of the parameters [2], making it very efficient in terms 
of both execution speed and memory consumption.<br/>


### 2D density maps

2D density maps are a visual representation of how the density of the membrane change
along the simulation time.
These are computed for each leaflet constituting the bilayer,
giving insight on the system phase. <br/>


### Enrichment/Depletion Analysis
 
When proteins are embedded in a lipid bilayer they modulates their local lipid environment.<br/> 
This happens through enrichment or depletion of specific lipid residues,that may result 
in thickness and curvature modifications.<br/>
This step was implemented as presented in [3].

### Diffusion movements

With this analysis the motions of the lipid residues of the system are investigated. 
All the coordinates are extracted to define “maps” that can be used to explore specific
lipid clusters.<br/>

### Order parameter

The deuterium order parameter is a measure for the orientational mobility of the bonds 
between the carbon and hydrogen atoms of the acyl chain of the lipids.
It is  used for estimating the overall order of the membrane and
details of the conformations that the atoms in the lipid tails adopt. <br/>
The calculation of the Order Parameter is done for full-atom systems taking into account the SCH parameter based on the algorithm originally developed by J. Melcr. with the contribution from  H. Antila for NMRlipids project and readapted for the purpose of this work [4].
For coarse-grained systems, the Order Parameter SCC is computed. Its implementation is based on the usage of the Lipyphilic class SCC() [5]. <br/>
 
### Membrane Curvature
This analysis tool allows the investigation of how lipid composition and membrane proteins generate membrane curvature. It provides visual representation of how the shape and curvature of the membrane changes along the simulation time. The surface, mean curvature and Gaussian curvature values are computed for each leaflet constituting the bilayer, and the graphical outputs consist of average plots for all frames in the trajectory as well as single- and multi-frame graphical outputs. The latter reveal the frame-by-frame alterations in membrane curvature.

### Protein-lipid interaction
Tracks and measures the contacts between the defined protein group and the lipids, calculating the depletion-enrichment factor, that compares the concentration of lipid species around the protein with that in the system.

### Lipid scrambling
Tracks and measures the z-coordinate of lipids, comparing it to an estimation of the membrane surface. Allows to easily identify lipids changing from one leaflet to another and tracks their contacts with the protein (if present). Both estimation of the surface and tracking of z-coordinate are performed based on the lipid headgroup position.


## Installing LipidDyn UNIX users

These instructions will get you a copy of the project running on your machine.


### Setup of the virtual environment

In a directory of the user's choice open a terminal and type:
```
virtualenv "LipidDyn_env" -p /usr/bin/python3.8
```

Once created activate the environment with:
```
source LipidDyn_env/bin/activate
```

**N.B.** <br/>
To check if the virtual environment is activated see if 
the name of the enviroment is visible on the terminal

### Prerequisites

Install numpy with:
```
pip install numpy==1.23.5
```
```
pip install packaging
```
Install LipidDyn with:

```
git clone https://github.com/ELELAB/LipidDyn.git
cd LipidDyn/
python setup.py install
cd ../
```
**N.B.**<br/>
All the softwares and libraries required to work are automatically installed when running ```python setup.install``` and they are listed below:
```
FatSlim
MDAnalysis
matplotlib
seaborn
pyyaml
lipyphilic
```


In case you want to deactivate the virtual environment just type the following
command into the terminal.

```
deactivate
``` 
## Installing LipidDyn MacOS users

To install the pipeline on pc with MacOS installed, the user should follow these passages:

### Install Anaconda environment

To install Anaconda environment we refer to the official page of Anaconda with alle step by step instructions:
```
https://docs.anaconda.com/anaconda/install/mac-os/
```
Once installed the user needs to activate the environment on a terminal with :
```
conda activate base
```
### Setup of the virtual environment

Create the virtual environment using :
```
python -m venv "LipidDyn_env"
```

Once created activate the environment with:
```
source LipidDyn_env/bin/activate
```

Then proceed normally as explained above : 
```
pip install numpy
pip install packaging
git clone https://github.com/ELELAB/LipidDyn.git
cd LipidDyn/
python setup.py install
cd ../
```

**N.B.**<br/>
In order to use LipidDyn and all its tools the user must activate the virtual environment each time.

## Configuration file
LipidDyn provides the usage of configuration files designed to work with full-atom and coarse-grained systems. The supported format is ```.yml```.

Besides giving insight on how the lipid categories and atoms selections are handled by the program, it can also be edited by the user to run the analysis according to the system requirements. <br/>

Referring to the configuration files provided by LipidDyn as templates and designed for  full-atom (FA) and coarse-grained (CG) systems, we describe here the structure and the usage of the different levels by which they are defined: <br/>

1) ```lipids:``` the user can add lipid categories if missing in this  section of the config file <br/>
2) ```protein: ``` protein atoms selection can be specified 
3) lipid ```headgroups: ``` for each lipid category the user has to define the lipid headgroup atoms for the definition of the bilayer system within the program;  <br/>
4) lipid ```bulk_ratio```:  for each lipid type the user can insert the ratio in bulk of each lipid for enrichment maps calculation. If not provided it will be computed by the program. Since the ratio in bulk is a value taken only when computing the enrichment analysis it could be specified when a protein is present in the system <br/>
5) lipid ```sn1: ``` and ```sn2: ``` acyl chains: for CG systems only, the user can define the CC-atom pairs to be used for the SCC order parameter calculation. <br/>
6) ```forcefield:``` specify if the system analysed has been simulated in CG or FA forcefield <br/>
7) ```interaction_domains:``` additive list of items that indicate which residues of the protein to use in the protein-lipid interaction and scrambling module. Two dash-separated integers are interpreted as a range of residue indices. Any other string is interpreted as a MDAnalysis selection command. Defaults to None, which indicates to use the full protein. N.B. Remember to define ```protein```, as the selection will be performed on the resulting group.
8) ```interaction_cutoff:``` float indicating interaction cutoff. Atoms within the cutoff distance (in Å) will be regarded as contact for the protein-lipid interaction module. Defaults to 6 Å.
9) ```interaction_by_res:``` any. If defined with a True value (e.g., 1), computes depletion-enrichment factor for each residue in the selection for the protein-lipid interaction module. This can be very cpu intensive and slow for large selections.
10) ```scrambling_lipids``` list of lipid residue names to analyse for scrambling. If not defined, every lipid is analyzed.
11) ```scrambling_domains:``` similar to ```interaction_domains```, employed by the scrambling module to compute lipid-protein contacts.
12) ```scrambling_cutoff:``` similar to ```interaction_cutoff```, employed by the scrambling module to compute lipid-protein contacts.
13) ```scrambling_upper_z:``` and ```scrambling_upper_rel_z:``` float boundary to define scrambling. If at a given frame, a lipid is above the given value, it is considered to be on the upper side of the membrane. If both are defined, both need to be fulfilled (logical and); if neither is defined, and all lipids will be considered. z (absolute) is the coordinate in Å, while relative z takes the leaflets as reference (0 is the lower leaflet, 0.5 the middle, and 1 the upper leaflet)
14) ```scrambling_lower_z:``` and ```scrambling_lower_rel_z:``` float boundary to define scrambling. If at a given frame, a lipid is below the given value, it is considered to be on the lower side of the membrane. If both are defined, both need to be fulfilled (logical and); if neither is defined, and all lipids will be considered. z (absolute) is the coordinate in Å, while relative z takes the leaflets as reference (0 is the lower leaflet, 0.5 the middle, and 1 the upper leaflet)
15) ```scrambling_upper_ratio:``` and ```scrambling_lower_ratio:``` portion of the trajectory (float between 0 and 1) that a lipid needs to be in the upper or lower side of the membrane (see previous parameters) to be considered scrambling. If not defined, defaults to one frame.

**N.B.** <br/> 
To allow LipidDyn to run correctly it is neccessary to keep the selection strings as they are presented in the config file templates provided here. 

## Running the pipeline

This section illustrates how to run some simple analysis. At the present moment only GROMACS trajectory and topology files are supported.<br/>
Prior to run the analysis the user needs to process the trajectory and the topology files for the periodic boundary conditions (and eventually skip frames along the trajectory or use a smaller part over the whole simulation length). 
Both files must have the same number of atoms.<br/>

### Input files for testing the pipeline

This section describes which files to use as test cases for both CG and FA systems. 
To test the pipeline all the required input files can be retrieved in  ```/LipidDyn/tests/data/full_atom/```  and in ```/LipidDyn/tests/data/coarse_grained/``` for FA and CG systems respectively. 

In each of these locations the user can find two subfolders called ```membrane/``` and ```membrane_protein/``` in which the trajectory (file.xtc), topology (file.gro) and configuration (file.yml) files of the test cases with lipid-only and lipid + protein systems are stored respectively.

On each of these folders the user can also find:

1) A bash script ```run_pipeline.sh``` containing all the commands used to run LipidDyn on the test cases, including the visualization tool command-lines
2) The outputs from each analysis
3) The outputs from the visualization tools ready to be visualized

**N.B.**<br/>
The simulations here reported as example are 51 frames long.


### Basic Usage 

```
LipidDyn -f file.xtc/trr -s file.gro -g file.yml -a -n "n" -c
```
This command will run the full set of analysis (-a) and store all the output files in the working directory, using "n" cores and clean (-c) all intermediate files. <br/>
The suggested minimum number of cores to run the analysis is ``` -n 4```, but it could be increased for longer trajectories. <br/>
The output files will be organized in different folders each one representing a set of analysis. By default, the program will skip an analysis module if it finds all its output files already in the output folder. In order to force the program to repeat the analysis and overwrite the output, use the flag ```-ow```.
<br/>


**N.B.**<br/>
If you want to obtain raw data for apl and thickness for each frame use the ```-r``` flag. Raw data will be produced in ```.csv``` format  and stored inside two subfolders under   the ```Fatslim/``` directory and called  ```fatslim_apl``` and  ```fatslim_thickness```.  <br/>


In the case of a protein embedded in the lipid bilayer use the flag ```-p```:

```
LipidDyn -f file.xtc/trr -s file.gro -g file.yml -a -n "n" -c -p
```

### Single analysis usage

This section lists the command-lines to run LipidDyn on each single analysis <br/>

1) APL and Thickness calculation with fatslim:  
```
LipidDyn -f file.xtc/.trr -s file.gro -g file.yml -fatslim -n "n"  -c
```
2) 2Density maps calculation for upper and lower leaflets:
```
LipidDyn -f file.xtc/.trr -s file.gro -g file.yml -2d -n "n"  -c
```
3) Enrichment maps calculation for each lipid category in upper and lower leaflets:
```
LipidDyn -f file.xtc/.trr -s file.gro -g file.yml -enr -n "n" -p -c
```
4) Diffusion movements calculation for upper and lower leaflets:
```
LipidDyn -f file.xtc/.trr -s file.gro -g file.yml -mov -n "n" -c
```
5) Order Parameter calculation for each lipid category:
```
LipidDyn -f file.xtc/.trr -s file.gro -g file.yml -op -n "n" -c
```
6) Membrane Curvature calculation for upper and lower leaflets:
```
LipidDyn -f file.xtc/.trr -s file.gro -g file.yml -mc -n "n" -c
```
7) Protein-lipid interaction: calculation of depletion-enrichment factor:
```
LipidDyn -f file.xtc/.trr -s file.gro -g file.yml -int -n "n" -p -c
```
8) Lipid scrambling: characterization of lipids undergoing scrambling:
```
LipidDyn -f file.xtc/.trr -s file.gro -g file.yml -int -n "n" -c
```


### Visualization of data

For the visualization of the output data LipidDyn includes a set of tools for graphical representation.
All the output are stored in the working directory

#### 1) Thickness and Area per Lipid

All the output of this step are stored in the directoyry ```Fatslim/```. It contains 
the .xvg files ```apl.xvg```  and ```thickness.xvg``` with the values 
for upper, lower and the entire membrane along the simulation time. <br/>
We can use the ```profiler``` tool to plot the data. 
Access the folder with the terminal and run :   

```
profiler -p Fatslim -out thickness.pdf -d thickness_distribution.pdf -plot thickness -spe
profiler -p Fatslim -out apl.pdf -d apl_distribution.pdf -plot apl -spe
```
N.B.<br/>
In some cases, there could be thickness or apl uncommon values (i.e. too low or negative values) computed for one or more frames. This could compromise the overall plot layout. The user can handle this by inserting the flag ```-th``` that specifies the threshold below which discard the values in the plotting.
Furthermore, the user is able to make visualization for each lipid species and a combination there of by adding the argument ```-spe```

#### 2) Density Maps

The output of this step can be found in the directory ```2Dmaps/```. 
The folder contains ```lower_leaflet_2dmap.dat``` and ```upper_leaflet_2dmap.dat``` density files. <br/>
We can use the ```dmaps``` tool to plot the data. 
Access the folder with the terminal and run: 

```
dmaps -i 2D_maps/lower_leaflet_2dmap.dat -o 2D_maps/lower_leaflet_2dmap.pdf
dmaps -i 2D_maps/upper_leaflet_2dmap.dat -o 2D_maps/upper_leaflet_2dmap.pdf
```


#### 3) Enrichment/Depletion 

In the case of a embedded protein the output of this step can be found in the directory ```Enrichment/```. 
The folder contains multiple file ```.dat``` ( as in the 2dmaps) named after the lipid residues constituting
the constituting the upper and lower leaflets **i.e** "POPC_enrich_lower_leaflet_enrich.dat", "SSM_upper_leaflet_enrich.dat" etc...  
We can use again the ```dmaps``` tool to plot the data with ```-enr``` flag to warn the tool that is an erichment
plot.
Access the folder with the terminal and run : 

```
dmaps -i lipid_residue_upper_leaflet_enrich.dat -o custom_name.pdf -enr
dmaps -i lipid_residue_lower_leaflet_enrich.dat -o custom_name.pdf -enr
```
**N.B.** <br/>
When running 2Dmaps and enrichment calculations on long simulation systems it is suggested to use a small subset of frames (i.e. last 1 us from the entire trajectory) according to the user requirments, so as to obtain a more clear output map.

**N.B.**<br/>
If you want to gain a better visualization of the densities of the system for both 2Dmaps and Enrichment analysis, you can adjust the upper and lower limits of the densities specifying the flag ```-min``` and ```-max```

#### 4) Diffusion Movements

The output of this step can be found in the directory ```Diffusion_movements```. 
The folder contains ```Lower_leaflet_coordinates.dat``` and ```Upper_leaflet_coordinates.dat```.
These files contains the x and y coordinates of all the lipid residue constituting the bilayer.<br/>
We can use the ```diffusion``` tool to plot the data. Depending on the of bilayer composition we need to select the ```-t``` flag that takes as arguments ```ho``` (homogeneous) or ```he``` (heterogeneous). <br/>
Access the folder with the terminal and run : 

```
diffusion -i Lower_leaflet_coordinates.dat -o suffix_output_name -t  -[he/ho] 
```
**N.B.** <br/>
When running Diffusion Movements analysis on long simulation systems it is suggested to skip some frames to avoid plotting too many coordinate values. <br/>

In the case of heterogeneous systems, the plot will produce a ```.pdf``` file for each lipid category by default. If you want to produce a single plot containing all the lipid categories use the flag ```-m```.

#### 5) Order Parameter

The output of this step can be found in the directory ```Order_parameter```. 
The folder contains ```.csv``` files named after the composition of the membrane **i.e** "Order_parameter_POPC.csv", "Order_parameter_SSM.csv".<br/>
According to the force-field used (CG or FA) it is required to specify the ```-s``` flag which takes the argument "sch" in case of full-atom systems and "scc" in case of coarse-grained systems.<br/>
Access the folder with the terminal and run : 

```
ordpar -i Order_Parameter_lipid_residue.csv -o custom_name.pdf -s [scc/sch]
```

**N.B.**<br/>
If you want to customize the range of order parameter values to visualize in the ```ordpar``` output plot, you can adjust the upper and lower limits of the SCH/SCC values specifying the flag ```-min``` and ```-max```

#### 6) Membrane Curvature 

The output of this step can be found in the directory ```curv/```.
The folder contains multiple file ```.dat``` named after the lipid residues constituting the upper and lower leaflets **i.e** "up_Avg_surface.dat", "low_AF_mean_curvature.dat" etc... "Avg" refers to the average data taken over n frames in the trajectory. "AF" stands for all frames and contains the curvature values for each individual frame.  
We can use the ```curvature``` tool to plot the data with ```-plot``` flag to decide which type of plot to produce.  
Access the folder with the terminal and run :

For average plots:
The ```-plot``` flag can be used to plot the basic surface, smooth surface, mean curvature, and Gaussian curvature as well as surface, mean and Gaussian curvature side-by-side in one figure for each leaflet.
```
curvature -l1 up_Avg_surface.dat -l2 low_Avg_surface.dat -o custom_name.pdf -plot basic
curvature -l1 up_Avg_surface.dat -l2 low_Avg_surface.dat -o custom_name.pdf -plot smooth 
curvature -l1 up_Avg_mean_curvature.dat -l2 low_Avg_mean_curvature.dat -o custom_name.pdf -plot mean 
curvature -l1 up_Avg_gaussian_curvature.dat -l2 low_Avg_gaussian_curvature.dat -o custom_name.pdf -plot gaussian 
curvature -l1 low_Avg_surface.dat -l2 low_Avg_mean_curvature.dat -l3 low_Avg_gaussian_curvature.dat -o custom_name.pdf -plot 3_curvatures -t 'custom title lower leaflet'
curvature -l1 up_Avg_surface.dat -l2 up_Avg_mean_curvature.dat -l3 up_Avg_gaussian_curvature.dat -o custom_name.pdf -plot 3_curvatures -t 'custom title upper leaflet' 

```

For single frame plots:
The ```-frame``` flag can be used to plot a specific frame in the trajectory.
```
curvature -l1 up_AF_surface.dat -l2 low_AF_surface.dat -o custom_name.pdf -plot basic -frame [int] 
curvature -l1 up_AF_surface.dat -l2 curv/low_AF_surface.dat -o custom_name.pdf -plot smooth -frame [int]
curvature -l1 up_AF_gaussian_curvature.dat -l2 low_AF_gaussian_curvature.dat -o custom_name.pdf -plot gaussian -frame [int] 
curvature -l1 up_AF_mean_curvature.dat -l2 low_AF_mean_curvature.dat -o custom_name.pdf -plot mean -frame [int]

```

For multi-frame plots:
The ```-plot``` flag can be used to plot the first, middle and last six frames of the trajectory for surface. 
```
curvature -l1 up_AF_surface.dat -l2 low_AF_surface.dat -o custom_name_upper.pdf -o2 custom_name_lower.pdf -plot first_6
curvature -l1 up_AF_surface.dat -l2 low_AF_surface.dat -o custum_name_upper.pdf -o2 custom_name_lower.pdf -plot middle_6
curvature -l1 up_AF_surface.dat -l2 low_AF_surface.dat -o custom_name_upper.pdf -o2 custom_name_lower.pdf -plot last_6

```

**N.B.** <br/>
When running membrane curvature calculations on long simulation systems it is suggested to use a small subset of frames for the average plots (i.e. last 1 us from the entire trajectory) according to the user requirments, so as to obtain a more clear output where the curvature has not evened out to much.

**N.B.**<br/>
If you want to gain a better visualization of the curvature of the system for mean and Gaussian curvature plots, you can adjust the upper and lower limits of the plots color bar by specifying the flag ```-lim```. 

**N.B.**<br/>
The user can deviate from default settings by using the ```-t``` argument to define the title of the plot, ```-c``` for the color spectrum on the colorbar, and ```-label``` for the label on the colorbar. For the outputs of mean and Gaussian curvature, the user can also specify the number of contour lines via ```-level```.

#### 7) Protein-lipid interaction

The output of this step can be found in the directory ```Interaction/```.
The folder contains ```d_e.tsv```, and, if specified, ```d_e_by_res.tsv```. Contains the depletion-enrichment factor for each lipid with respect to the protein, per-frame.

Access the folder with the terminal and run:

```
interaction -i d_e.tsv -l <lipid1> <lipid2> ... -o <ouput_prefix>
```

Output files include <ouput_prefix>_line.pdf and <ouput_prefix>_box.pdf with the evolution and distribution of the depletion enrichment factor, respectively. Note that the visualization of by-residue data is not implemented.

#### 8) Lipid scrambling

The output of this step can be found in the directory ```Scrambling/```.
The folder contains ```z_data.tsv``` with the data of each lipid. It is used to filter the lipids that undergo scrambling. To be loaded as a numpy.ndarray, must be reshaped to (n_frames, n_lipids, 4) to form a 3D array (per-frame, per-lipid headgroup). It contains time (ps), estimated z coordinate for upper and lower leaflet, and z coordinate of the headgroup.
Further results of the module are in the subfolder ```lipids/```, where there is a file for every lipid found to be scrambling, according to the input parameters (<resname>_<resid>.tsv). This files contain the fraction of ```z_data.tsv``` corresponding to a single lipid.

To plot the results, access the folder with the terminal and run:

```
scrambling -i lipids/* -o <ouput_name>
scrambling -i lipids/<resname>_<resid>.tsv -o <ouput_name> -c
```

-a displays the absolute z coordinate instead the relative z to the leaflets
-l hides the leaflets positions in the plot
-c creates a heatmap in the background representing the contacts of the lipid with the protein (if the analysis had a defined protein). Only available for one lipid.
-r allows to specify a range of protein residues if ```-c``` is called. Useful to zoom in or out for better visualization


## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details

## Citation

If you use this software, please cite it as follows:

Unravelling membrane properties at the organelle-level with LipidDyn

Simone Scrima, Matteo Tiberti, Alessia Campo, Elisabeth Corcelle-Termeau, Delphine Judith, Mads Moller Foged, Knut Kristoffer Bungaard Clemmens, Sharon A Tooze, Marja Jaattela, Kenji Maeda, Matteo Lambrughi*, Elena Papaleo*

Comput. Struct. Biotechnol. J. 2022 Jun 30;20:3604-3614. doi: 10.1016/j.csbj.2022.06.054 PMID: 35860415 PMCID: PMC9283888`

```


## References

[1] Piggot T.J. et al.  (2012) Molecular dynamics simulations of phosphatidylcholine membranes: a comparative force field study. J. Chem. Theory Comput., 8, 4593–4609. https://doi.org/10.1021/ct3003157<br/>

[2] Sébastien Buchoux, FATSLiM: a fast and robust software to analyze MD simulations of membranes, Bioinformatics, Volume 33, Issue 1, 1 January 2017, Pages 133–134, https://doi.org/10.1093/bioinformatics/btw563 <br/>

[3] Noemi Jiménez-Rojo et al. Conserved function of ether lipids and sphingolipids in the early secretory 15 pathway (2020), Current Biology,30, 19, P3775-3787.E7. https://doi.org/10.1101/2019.12.19.881094 <br/>

[4] https://github.com/NMRLipids/MATCH/scripts

[5] Smith, Paul and Lorenz, Christian D. (2021), LiPyphilic: A Python Toolkit for the Analysis of Lipid Membrane Simulations, Journal of Chem. Theory Comput., 17,9, 5907-5919. https://doi.org/10.1021/acs.jctc.1c00447






