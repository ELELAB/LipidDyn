# Import the environment prior to run this script

#LipidDyn -t heterogeneous_membrane_protein.xtc -f heterogeneous_membrane_protein.gro -g full_atom.yml -prot -all -ncore 4 -c -raw

# Plotting 

# Plot density maps
dmaps -i 2D_maps/lower_leaflet_2dmap.dat -o 2D_maps/lower_leaflet_2dmap.pdf
dmaps -i 2D_maps/upper_leaflet_2dmap.dat -o 2D_maps/upper_leaflet_2dmap.pdf

#Enrichment
dmaps -i Enrichment/POPC_lower_leaflet_enrich.dat -o Enrichment/POPC_lower_leaflet_enrich.pdf -enr
dmaps -i Enrichment/POPC_upper_leaflet_enrich.dat -o Enrichment/POPC_upper_leaflet_enrich.pdf -enr
dmaps -i Enrichment/SSM_lower_leaflet_enrich.dat -o Enrichment/SSM_lower_leaflet_enrich.pdf -enr
dmaps -i Enrichment/SSM_upper_leaflet_enrich.dat -o Enrichment/SSM_upper_leaflet_enrich.pdf -enr

# Plot APL and Thicknes
profiler -i Fatslim/thickness.xvg -out Fatslim/thickness.pdf -d Fatslim/distribution_thickness.pdf -plot thick
profiler -i Fatslim/apl.xvg -out Fatslim/apl.pdf -d Fatslim/distribution_apl.pdf -plot apl

# Diffusion
diffusion -i Diffusion_movements/Lower_leaflet_coordinates.dat -o Diffusion_movements/Lower_leaflet_coordinates -t he
diffusion -i Diffusion_movements/Lower_leaflet_coordinates.dat -o Diffusion_movements/Lower_leaflet_coordinates -t he -m 
diffusion -i Diffusion_movements/Upper_leaflet_coordinates.dat -o Diffusion_movements/Upper_leaflet_coordinates -t he
diffusion -i Diffusion_movements/Upper_leaflet_coordinates.dat -o Diffusion_movements/Upper_leaflet_coordinates -t he -m

# Order Parameter
ordpar -i Order_Parameter/Order_Parameter_SSM.csv -o Order_Parameter/Order_Parameter_SSM.pdf -s sch
ordpar -i Order_Parameter/Order_Parameter_POPC.csv -o Order_Parameter/Order_Parameter_POPC.pdf -s sch
