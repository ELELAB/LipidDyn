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
Species,abbrev,lm_id,name,Exceptions,Charmm36,concentration
CL 68:2,CL 68:2,LMGP12010081,"CL(1'-[16:0/18:1(9Z)],3'-[16:0/18:1(9Z)])",False,POCL2,0.0
CL 68:2,CL 68:2,LMGP12010081,"CL(1'-[16:0/18:1(9Z)],3'-[16:0/18:1(9Z)])",False,POCL1,0.0
Cer 42:1;2,Cer 42:1;O2,LMSP02010012,Cer(d18:1/24:0),True,CER240,0.0675962090464608
Cer 42:2;2,Cer 42:2;O2,LMSP02010009,Cer(d18:1/24:1(15Z)),True,CER241,0.134158391912655
PA 32:2,PA 32:2,LMGP10010969,PA(16:1(9Z)/16:1(9Z)),False,DYPA,0.0283727439444567
PA 34:1,PA 34:1,LMGP10010032,PA(16:0/18:1(9Z)),False,POPA,0.0526526227726728
....
CL 68:5,NA,NA,NA,False,NA,0.0
---------------------------------------

Here are stored the information coming from the original lipidomics dataset i.e column "Species",the abbreviations of the chains, the lipidmaps ID of the lipid,the name, the category of the lipid, the force field name and ultimately the concentration from
the original dataset.

The "-og" group all the results from above and organized them in the following manner:


---------------------------------------
Species,ID,charmm36,concentration
CL 64:3,NA,['NA'],['NA'],0.0
CL 64:4,"CL(1'-[16:1(9Z)/16:1(9Z)],3'-[16:1(9Z)/16:1(9Z)])",['NA'],"['TYCL2', 'TYCL1']",0.0
CL 66:2,CL 66:2,"['LMGP12010655', 'LMGP12010025', 'LMGP12010115', 'LMGP12010010']",['NA'],0.0
.....
---------------------------------------

Here there are stored the information from the original dataset i.e column "Species", the lipid ID of the lipids, the force field
name of the lipid and the concentration.
In case of more than one species more name will be stored in the list i.e è ['DUPC', 'DLiPC].

In the case of not matching species values are replaced with NA both in the non-grouped and grouped output.
