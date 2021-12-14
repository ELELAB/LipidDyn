# LipidDyn: A computational microscope to scrutinize membrane properties at the organelle-level

We have built an automated computational pipeline, implemented in Python 3, LipidDyn,
for the accurate analysis of the structural properties and dynamics of lipid bilayers
simulations that can be used to validate the ensembles against experimental data.<br/>
<br/>
The framework is divided into different modules that can be run separately, 
depending on the user choice.<br/> Each performs a different set of analysis on both full-atom and coarse-grained systems.<br/>
Moreover, the pipeline is able to account for embedded proteins into the membrane. <br/>
<br/>
At the present state is designed to work with GROMACS files.


## Set of Analysis 

### Thickness and Area per Lipid

The membrane thickness and area per lipid are two parameters that are used to 
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
This happen through enrichment or depletion of specific lipid residues,that may result 
in thickness and curvature modifications.<br/>
This step was implemented as presented in [3].

### Diffusion movements

With this analysis the motions of the lipid residues of system are investigated. 
All the coordinates are extracted to defines “maps” that can be used to explore specific
lipid clusters.<br/>

### Order parameter

The deuterium order parameter is a measure for the orientational mobility of the bonds 
between the carbon and hydrogen atom of the acyl chain of the lipid.
It is  used for estimating the overall order of the membrane and
details of the conformations that the atoms in the lipid tails adopt. <br/>
The calculation of the Order Parameter is done for full-atom system taking into account the SCH parameter based on the algorithm originally developed by J. Melcr. with the contribution from  H. Antila 
for NMRlipids project and readapted for the purpose of this work [4].
For coarse-grained systems, the Order Parameter SCC is computed. Its implementation is based on the usage of the Lipyphilic class SCC() [5] 
 <br/>

## Installing LipidDyn

These instructions will get you a copy of the project running on your machine.

### Prerequisites

Softwares and libraries required to work:

```
FatSlim
MDAnalysis
matplotlib
numpy
seaborn
pyyaml
lipyphilic
```

### Setup of the virtual environment

In a directory of the user's choice open a terminal and type:
```
virtualenv "LipidDyn_env" -p /usr/bin/python3.8
```

Once created activate the environment with:
```
source LipidDyn_env/bin/activate
```

**N.B** <br/>
To check if the virtual environment is activated see if 
the name of the enviroment is visible on the terminal


Install LipidDyn with:

```
git clone https://github.com/ELELAB/LipidDyn.git
cd LipidDyn/
python setup.py install
cd ../
```

**N.B**<br/>
The setup.py will also install all the python library required by LipidDyn to work

In case you want to deactivate the virtual environment just tipe the following
command into the terminal.

```
deactivate
``` 

**N.B**<br/>
In order to use LipidDyn and all its tools the user must activate the virtual environment each time.

## Configuration file
LipidDyn provides the usage of configuration files designed to work with FA and CG systems. The supported format is .yml.

Besides giving insight on how the system lipid categories and atoms definition are handled by the program, it can be edited by the user to run the analysis according to the system requirments. <br/>

Referring to the configuration files  provided by LipidDyn as templates and designed for FA and CG systems, we describe here the structure and the usage of the different levels by which is defined: <br/>
1) ```lipids:``` the user can add lipid categories if missing in this  section of the config file <br/>
2) ```protein: ``` protein atoms selection can be specified 
3) lipid ```headgroups: ``` for each lipid category the user has to define the lipid headgroup atoms for the definition of the bilayer system within the program;  <br/>
4) lipid ```bulk_ratio```:  for each lipid type the user can insert the ratio in bulk of each lipid for enrichment maps calculation. If not provided it will be computed by the program <br/>
5) lipid ```sn1: ``` and ```sn2: ``` acyl chains: for CG systems only, the user can define the CC-atom pairs to be used for the SCC order parameter calculation. <br/>
6) ```forcefield:``` specify if the system analysed has been simulated in CG or FA forcefield <br/> 

**N.B.** <br/> 
To allow LipidDyn to run correctly it is neccessary to keep the selection strings as they are presented by the config templates provided here. 

## Running the pipeline

This section illustrates how to run some simple analysis. At the present moment only 
GROMACS trajectory and topology filesare supported.<br/>
Prior to run the analysis the user needs to processed the trajectory and the topology 
file for the periodic boundary conditions (and eventually skip frames along the trajectory or use a smaller part over the whole simulation length). 
Both files must have the same number of atoms.<br/>

 

### Basic Usage 

```
LipidDyn -t ["".xtc/.trr] -f ["".gro/.tpr] -g ["".yml] -all -ncore "n"  -c
```
This command will run the full set of analysis (-all) and store all the output files in the working directory, using "n" cores and clean (-c) all intermediate files. <br/> 
The output files will be organized in different folders each one representing a set of analysis. 
<br/>
In the case of a protein embedded in the lipid bilayer use the flag -prot:

```
LipidDyn -t ["".xtc/.trr] -f ["".gro/.tpr] -g ["".yml] -all  -ncore "n" -prot -c
```

### Single analysis usage
```
LipidDyn -t ["".xtc/.trr] -f ["".gro/.tpr] -g ["".yml] -fatslim -ncore "n"  -c
```
```
LipidDyn -t ["".xtc/.trr] -f ["".gro/.tpr] -g ["".yml] -2d -ncore "n"  -c
```
```
LipidDyn -t ["".xtc/.trr] -f ["".gro/.tpr] -g ["".yml] -enr -ncore "n" -prot -c
```

```
LipidDyn -t ["".xtc/.trr] -f ["".gro/.tpr] -g ["".yml] -mov -ncore "n" -c
```
```
LipidDyn -t ["".xtc/.trr] -f ["".gro/.tpr] -g ["".yml] -ordpar -ncore "n" -c
```
### Visualization of data

For the visualization of the output data LipidDyn includes a set of tools for graphical representation.
All the output are stored in the working directory

#### 1) Thickness and Area per Lipid

All the output of this step are stored in the directoyry ```Fatslim/```. It contains 
the .xvg files ```apl.xvg```  and ```thickness.xvg``` with the values 
for upper, lower and the entire membrane along simulation time. <br/>
We can use the ```profiler``` tool to plot the data. 
Access the folder with the terminal and run :   

```
profiler -i apl.xvg -out area_per_lipid.pdf -d distribution.pdf -plot apl
profiler -i thickness.xvg -out thickness.pdf -d distribution.pdf -plot thick 
```
N.B.<br/>
In some cases, there could be thickness or apl uncommon values (i.e. too low or negative values) computed for one or more frames. This could compromise the overall plot layout. The user can handle this by inserting the flag ```-th``` that specifies the threshold below which discard the values in the plotting

#### 2) Density Maps

The output of this step can be found in the directory ```2Dmaps/```. 
The folder contains ```lower_leaflet_2dmap.dat``` and ```upper_leaflet_2dmap.dat``` 
density files. <br/>
We can use the ```dmaps``` tool to plot the data. 
Access the folder with the terminal and run : 

```
dmaps -i lower_leaflet_2dmap.dat -o lower_leaflet_2dmap.pdf 
dmaps -i upper_leaflet_2dmap.dat -o upper_leaflet_2dmap.pdf
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
When running 2Dmaps and enrichment calculation on long simulation systems it is suggested to use a small subset of frames (i.e. last 1 us from the entire trajectory) according to the user requirments so as to obtain a more clear output map

#### 4) Diffusion Movements

The output of this step can be found in the directory ```Diffusion_movments```. 
The folder contains ```Lower_leaflet_coordinates.dat``` and ```Upper_leaflet_coordinates.dat```.
These files contains the x and y coordinates of all the lipid residue constituting the bilayer.<br/>
We can use the ```diffusion``` tool to plot the data. Depending on the of bilayer composition we need to 
select the ```-ho``` (homogeneous) or ```-he``` (heterogeneous) flags.<br/>
Access the folder with the terminal and run : 

```
diffusion -i lipid_residue_enrich.dat -o custom_name.pdf -enr
```
**N.B.** <br/>
When running Diffusion Movments analysis on long simulation systems it is suggested to skip some frames to avoid plotting too many coordinates values 

#### 5) Order Parameter

The output of this step can be found in the directory ```Order_parameter```. 
The folder contains ```.csv``` files named after the composition of the membrane 
**i.e** "Order_parameter_POPC.csv", "Order_parameter_SSM.csv".<br/>
According to the force-field used (CG or FA) it is required to specify the -s flag which takes the argument "sch" in case of full-atom systems and "scc" in case of coarse-grained systems.<br/>
Access the folder with the terminal and run : 

```
ordpar -i Order_Parameter_lipid_residue.csv -o custom_name.pdf -s [scc/sch]
```

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details



## References

[1] Piggot T.J. et al.  (2012) Molecular dynamics simulations of phosphatidylcholine membranes: a comparative force field study. J. Chem. Theory Comput., 8, 4593–4609. https://doi.org/10.1021/ct3003157<br/>

[2] Sébastien Buchoux, FATSLiM: a fast and robust software to analyze MD simulations of membranes, Bioinformatics, Volume 33, Issue 1, 1 January 2017, Pages 133–134, https://doi.org/10.1093/bioinformatics/btw563 <br/>

[3] Noemi Jiménez-Rojo et al. Conserved function of ether lipids and sphingolipids in the early secretory 15 pathway (2020), Current Biology,30, 19, P3775-3787.E7. https://doi.org/10.1101/2019.12.19.881094 <br/>

[4] https://github.com/NMRLipids/MATCH/scripts

[5] Smith, Paul and Lorenz, Christian D. (2021), LiPyphilic: A Python Toolkit for the Analysis of Lipid Membrane Simulations, Journal of Chem. Theory Comput., 17,9, 5907-5919. https://doi.org/10.1021/acs.jctc.1c00447






