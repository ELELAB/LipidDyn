# LypidDyn: A computational microscope to scrutinize membrane properties at the organelle-level

We have built an automated computational pipeline, implemented in Python 3,LypidDyn,
for the accurate analysis of the structural properties and dynamics of lipid bilayers
simulations. <br/>
The framework is divided into four independent modules that can be run separately, 
depending on the user choice. Each performs a different analysis able to account 
for embedded proteins into the membrane. <br/>
LipidDyn has been developed to include parameters that can be used to validate 
simulation ensembles against experimental data[x].<br/>
One module estimates the thickness of the membrane ( usually defined as the 
distance between phosphorus atoms in two leaflets ) and the area per lipid 
(the surface of the cross-section of the cylindrical hydrocarbon part of 
the lipid), providing information about the fluidity of the system. <br/>
The software employed is FatSlim [x].<br/>
A second module computes 2D density maps, a visual representation of how 
the density of the membrane changes, on both the leaflets constituting the bilayer,
giving insight on the system phase. <br/>
The third module investigates the diffusion motions of the system. It extracts all
the coordinates of each lipid residue of the system and defines “maps” that can be 
used to explore specific lipid clusters.<br/>
Ultimately a fourth module estimates the deuterium order parameter, a measure for 
the orientational mobility of the bonds between the carbon and hydrogen atom of the
acyl chain of the lipid, used for estimating the overall order of the membrane and
details of the conformations that the atoms in the lipid tails adopt. <br/>
This computational platform has currently employed to study, how different 
compositions in sphingolipids affects the structural and dynamical properties 
of organelle-like membrane models, in collaboration with our colleagues in the
Unit of Cell Death and Metabolism[x]. <br/>


## Getting Started


These instructions will get you a copy of the project up and running on your 
local machine. 
The pipeline was developed using Python3.6.

### Prerequisites

Softwares and libraries required to work

```
python3.6
Gromacs5.X
FatSlim
GromacsWrapper
MDAnalysis
```

### Installing

Step by step series of examples to setup the environment

Setup of the virtualenv for python3.6

```
virtualenv "LypidDyn" -p /usr/bin/python3.6
```

Activate the virtual environment and install all the packages into 
the environment 

```
source LypidDyn/bin/activate
```

Install the FatSlim [x]

```
git clone https://github.com/FATSLiM/fatslim.git
cd fatslim/
python setup.py install
```

Install GromacsWrapper [x] with pip

```
pip install GromacsWrapper
```

Install MDAnalysis [x] with pip

```
pip install MDAnalysis
```

For the installation of Gromacs software we refer to [x].
To deactivate the virtual environment just tipe the following
command into the terminal

```
deactivate
``` 

## Running the tests

This section illustrates how to run the pipeline for two identical
system i.e ratio of lipid residues,the first pure membrane bilayer,
the other with a transmenbrane protein embedded.


## 1) Source the virtual environment 

```
source LypidDyn/bin/activate
```

### Running the pipeline for the first 100ns


```
python final_pipeline.py -t trajectory.xtc -f topology.tpr -g state.cpt -d Analysis_100ns -b 0 -e 100000 -all -nthread 4 -raw -c
```

This command runs the complete set of analysis of the pipeline (-all) for the first 100ns for the trajectory ( -b 0 -e 100000) and run on 4 different cores (-nthread 4) and stores all the output files in the Analysis_100ns/ directory using the trajectory file (-f trajectory.xtc),topology file (-f topology.tpr), checkpoint file (state.cpt). <br/>
The -raw and -c are used to produce the raw output for the thickness and apl calculation for further processing and to clean all the intermediate files (i.e ndx file) 

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.


## Authors

* **Simone Scrima** - *Initial work* - [SimoneScrima](https://github.com/SimoneScrima)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc




