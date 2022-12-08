# Import the environment prior to run this script

LipidDyn -f heterogeneous_membrane_protein.xtc -s heterogeneous_membrane_protein.gro -g full_atom.yml -p -a -n 4 -c -spe

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
profiler -p Fatslim -out thickness.pdf -d thickness_distribution.pdf -plot thickness -s
profiler -p Fatslim -out apl.pdf -d apl_distribution.pdf -plot apl -s

# Diffusion
diffusion -i Diffusion_movements/Lower_leaflet_coordinates.dat -o Diffusion_movements/Lower_leaflet_coordinates -t he
diffusion -i Diffusion_movements/Lower_leaflet_coordinates.dat -o Diffusion_movements/Lower_leaflet_coordinates -t he -m 
diffusion -i Diffusion_movements/Upper_leaflet_coordinates.dat -o Diffusion_movements/Upper_leaflet_coordinates -t he
diffusion -i Diffusion_movements/Upper_leaflet_coordinates.dat -o Diffusion_movements/Upper_leaflet_coordinates -t he -m

# Order Parameter
ordpar -i Order_Parameter/Order_Parameter_SSM.csv -o Order_Parameter/Order_Parameter_SSM.pdf -s sch
ordpar -i Order_Parameter/Order_Parameter_POPC.csv -o Order_Parameter/Order_Parameter_POPC.pdf -s sch
