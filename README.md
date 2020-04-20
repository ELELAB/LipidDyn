# LypidDyn: A computational microscope to scrutinize membrane properties at the organelle-level

We have built an automated computational pipeline, implemented in Python 3, LypidDyn,
for the accurate analysis of the structural properties and dynamics of lipid bilayers
simulations that can be used to validate the ensembles against experimental data.<br/>
<br/>
The framework is divided into different modules that can be run separately, 
depending on the user choice.<br/> Each performs a different set of analysis.<br/>
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
This step was implemented as presented in [3]

### Diffusion movements

With this analysis the motions of the lipid residues of system are investigated. 
All the coordinates are extracted to defines “maps” that can be used to explore specific
lipid clusters.<br/>

### Order parameter

The deuterium order parameter is a measure for the orientational mobility of the bonds 
between the carbon and hydrogen atom of the acyl chain of the lipid.
It is  used for estimating the overall order of the membrane and
details of the conformations that the atoms in the lipid tails adopt. <br/>
The original algorithm was originally developed by J. Melcr. with the contribution from  H. Antila 
for NMRlipids project and readapted for the purpose of this work [4] <br/>




## Installing LypidDyn

These instructions will get you a copy of the project running on your machine.

### Prerequisites

Softwares and libraries required to work:

```
FatSlim
MDAnalysis
matplotlib
numpy
```

### Setup of the virtual environment

In a directory of the user's choice open a terminal and type:
```
virtualenv "LypidDyn_env" -p /usr/bin/python3.6
```

Once created activate the environment with:
```
source LypidDyn_env/bin/activate
```

**N.B** 
To check if the virtual environment is activated see if 
the name of the enviroment is visible on the terminal

Install numpy library with :

```
pip install numpy
```

Install FatSlim software:

```
pip install fatslim
```

Install LypidDyn with:

```
git clone https://github.com/ELELAB/LipidDyn.git
cd LipidDyn/
python setup.py install
cd ../
```

**N.B**
The setup.py will also install all the python library required by LypidDyn to work

In case you want to deactivate the virtual environment just tipe the following
command into the terminal.

```
deactivate
``` 

**N.B**
In order to use LypidDyn and all its tools the user must activate the virtual environment each time.



## Running the pipeline

This section illustrates how to run some simple analysis. At the present moment only 
GROMACS trajectory and topology filesare supported.<br/>
Prior to run the analysis the user needs to processed the trajectory and the topology 
file for the periodic boundary conditions (and eventually skip frames along the trajectory). 
Both files must have the same number of atoms.<br/>

### Basic Usage 

```
LypidDyn -t ["".xtc/.trr] -f ["".gro/.tpr] -all -d Analysis -ncore "n"  -c
```
This command will run the full set of analysis (-all) and store all the output files in the directory
"Analysis", using "n" cores and clean (-c) all intermediate files. <br/> 
Inside this folder there will be different folders each one representing a set of analysis. 
<br/>
In case there of a protein embedded in the lipid bilayer specify use the flag -prot:

```
LypidDyn -t ["".xtc/.trr] -f ["".gro/.tpr] -all -d Analysis -ncore "n" -prot -c
```

### Visualization of data

For the visualization of the output data LypidDyn includes a set of tools for for graphical representation.
All the output are stored under  ```Analysis/``` (or under the directory with the user custom name).

#### 1) Thickness and Area per Lipid

All the output of this step are stored in the directoyry ```Fatslim/```. It contains 
the .xvg files ```apl.xvg```  and ```thickness.xvg``` with the values 
for upper, lower and the entire membrane along simulation time. <br/>
We can use the ```profiler``` tool to plot the data. 
Access the folder with the terminal and run :   

```
profiler -i apl.xvg -out area_per_lipid.png -plot apl -upp -low -memb
profiler -i thickness.xvg -out thickness.png -plot thick -upp -low -memb
```


#### 2) Density Maps

The output of this step can be found in the directory ```2Dmaps/```. 
The folder contains ```lower_leaflet_2dmap.dat``` and ```upper_leaflet_2dmap.dat``` 
density files. <br/>
We can use the ```dmaps``` tool to plot the data. 
Access the folder with the terminal and run : 

```
dmaps -i lower_leaflet_2dmap.dat -o lower_leaflet_2dmap.png 
dmaps -i upper_leaflet_2dmap.dat -o upper_leaflet_2dmap.png
```


#### 3) Enrichment/Depletion 

In the case of a embedded protein the output of this step can be found in the directory ```Enrichment/```. 
The folder contains multiple file ```.dat``` ( as in the 2dmaps) named after the lipid residues constituting
the bilayer **i.e** "POPC_enrich.dat", "SSM_enrich.dat" etc...  
We can use again the ```dmaps``` tool to plot the data with ```-enr``` flag to warn the tool that is an erichment
plot.
Access the folder with the terminal and run : 

```
dmaps -i lipid_residue_enrich.dat -o custom_name.png -enr
```

#### 4) Diffusion Movements

The output of this step can be found in the directory ```Diffusion_moments```. 
The folder contains ```Lower_leaflet_coordinates.dat``` and ```Upper_leaflet_coordinates.dat```.
These files containt the x and y coordinates of all the lipid residue constituting the bilayer.<br/>
We can use the ```diffusion``` tool to plot the data. Depending on the of bilayer composition we need to 
select the ```-ho``` (homogeneous) or ```-he``` (heterogeneous) flags.<br/>
Access the folder with the terminal and run : 

```
diffusion -i lipid_residue_enrich.dat -o custom_name.png -enr
```

#### 5) Order Parameter

The output of this step can be found in the directory ```Order_parameter```. 
The folder contains ```.csv``` files named after the composition of the membrane 
**i.e** "Order_parameter_POPC.csv", "Order_parameter_SSM.csv".<br/>
Access the folder with the terminal and run : 

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details

```
ordpar -i Order_Parameter_"".csv -o "".png
```

## References

[1] Piggot T.J. et al.  (2012) Molecular dynamics simulations of phosphatidylcholine membranes: a comparative force field study. J. Chem. Theory Comput., 8, 4593–4609. https://doi.org/10.1021/ct3003157<br/>

[2] Sébastien Buchoux, FATSLiM: a fast and robust software to analyze MD simulations of membranes, Bioinformatics, Volume 33, Issue 1, 1 January 2017, Pages 133–134, https://doi.org/10.1093/bioinformatics/btw563 <br/>

[3] Noemi Jiménez-Rojo et al. Conserved function of ether lipids and sphingolipids in the early secretory
15 pathway, BiorXiv ,https://doi.org/10.1101/2019.12.19.881094 <br/>

[4] https://github.com/NMRLipids/MATCH/scripts









