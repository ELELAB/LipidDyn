The lipid2MD tools was designed with the aim to connect experimental lipidomics data
to in silico data. With this tool we are able to connect categories of lipids from 
lipidomics data to the definitions of the lipids available in the all-atom force field charmm36. 
The first input file should be in the following format:

"species","ER"
"Cer 32:0;2",0
"Cer 32:1;2",0.0176653725206329
"Cer 34:0;2",0.00176727834288462
"Cer 34:1;2",0.158168168483941
"Cer 34:2;2",0.0287423986260666
"Cer 34:3;2",0
....

In one column we have the category of the lipid while on the other its concentration.
The second input file is the config2MD.yaml available in the current folder in which
the definitions of the lipids of charmm36 are stored.
To run the tool:

lipid2MD -i dataset_species_ER_median.csv -c config2MD.yaml -fa -o trial -og trial

The flag "-o" and "-og" specify the name of the output files. Only one or both
flag can be selected.

The "-o" flag specify an output in the following format :

---------------------------------------
Species,ID,name,abbrev,abbrev_chains,Charmm36,concentration
CL 68:2,LMGP12010006,"CL(1'-[16:0/18:1(11Z)],3'-[16:0/18:1(11Z)])",CL 68:2,CL 16:0_16:0_18:1_18:1,PVCL2,0.0
CL 68:2,LMGP12010081,"CL(1'-[16:0/18:1(9Z)],3'-[16:0/18:1(9Z)])",CL 68:2,CL 16:0_16:0_18:1_18:1,POCL1,0.0
CL 68:2,LMGP12010081,"CL(1'-[16:0/18:1(9Z)],3'-[16:0/18:1(9Z)])",CL 68:2,CL 16:0_16:0_18:1_18:1,POCL2,0.0
CL 72:8,LMGP12010001,"CL(1'-[18:2(9Z,12Z)/18:2(9Z,12Z)],3'-[18:2(9Z,12Z)/18:2(9Z,12Z)])",CL 72:8,CL 18:2_18:2_18:2_18:2,TLCL1,0.0
CL 72:8,LMGP12010001,"CL(1'-[18:2(9Z,12Z)/18:2(9Z,12Z)],3'-[18:2(9Z,12Z)/18:2(9Z,12Z)])",CL 72:8,CL 18:2_18:2_18:2_18:2,TLCL2,0.0
PA 32:0,LMGP10010027,PA(16:0/16:0),PA 32:0,PA 16:0_16:0,DPPA,0.0
PA 34:1,LMGP10010032,PA(16:0/18:1(9Z)),PA 34:1,PA 16:0_18:1,POPA,0.0256606363735779
....
---------------------------------------

Here are stored the information coming from the original lipidomics dataset i.e column "Species", the lipidmaps ID of the lipid,
the name, the category of the lipid, the abbreviations of the chains, the force field name and ultimately the concentration from
the original dataset.

The "-og" group all the results from above and organized them in the following manner:


---------------------------------------
Species,ID,charmm36,concentration
PC 34:1,['LMGP01010005'],['POPC'],6.73572203030144
PC 32:1,['LMGP01010566'],['PYPC'],4.38083838085729
PC 36:4,['LMGP01010937'],"['DUPC', 'DLiPC']",0.595177427560921
.....
---------------------------------------

Here there are stored the information from the original dataset i.e column "Species", the lipid ID of the lipids, the force field
name of the lipid and the concentration.
In case of more than one species more name will be stored in the list i.e è ['DUPC', 'DLiPC]
