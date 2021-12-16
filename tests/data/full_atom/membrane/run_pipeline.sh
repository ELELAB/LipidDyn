# Import the environment prior to run this script

LipidDyn -t heterogeneous_membrane.xtc -f heterogeneous_membrane.gro -g full_atom.yml -all -ncore 4 -c -raw

# Plotting 

# Plot density maps
dmaps -i 2D_maps/lower_leaflet_2dmap.dat -o 2D_maps/lower_leaflet_2dmap.pdf
dmaps -i 2D_maps/upper_leaflet_2dmap.dat -o 2D_maps/upper_leaflet_2dmap.pdf

# Plot APL and Thicknes
profiler -i Fatslim/thickness.xvg -out Fatslim/thickness.pdf -d Fatslim/distribution_thickness.pdf -plot thick
profiler -i Fatslim/apl.xvg -out Fatslim/apl.pdf -d Fatslim/distribution_apl.pdf -plot apl

# Diffusion
diffusion -i Diffusion_movements/Lower_leaflet_coordinates.dat -o Diffusion_movements/Lower_leaflet_coordinates -t he
diffusion -i Diffusion_movements/Lower_leaflet_coordinates.dat -o Diffusion_movements/Lower_leaflet_coordinates -t he -m
diffusion -i Diffusion_movements/Upper_leaflet_coordinates.dat -o Diffusion_movements/Upper_leaflet_coordinates -t he
diffusion -i Diffusion_movements/Upper_leaflet_coordinates.dat -o Diffusion_movements/Upper_leaflet_coordinates -t he -m

# Order Parameter
ordpar -i Order_Parameter/Order_Parameter_LSM.csv -o Order_Parameter/Order_Parameter_LSM.pdf -s sch
ordpar -i Order_Parameter/Order_Parameter_NSM.csv -o Order_Parameter/Order_Parameter_NSM.pdf -s sch
ordpar -i Order_Parameter/Order_Parameter_PSM.csv -o Order_Parameter/Order_Parameter_PSM.pdf -s sch
ordpar -i Order_Parameter/Order_Parameter_POPC.csv -o Order_Parameter/Order_Parameter_POPC.pdf -s sch
