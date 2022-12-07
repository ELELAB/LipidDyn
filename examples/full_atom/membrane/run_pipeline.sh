# Import the environment prior to run this script

LipidDyn -f heterogeneous_membrane.xtc -s heterogeneous_membrane.gro -g full_atom.yml -n 4 -c -r -fatslim -mov -op -2d 

# Plotting 

# Plot density maps
dmaps -i 2D_maps/lower_leaflet_2dmap.dat -o 2D_maps/lower_leaflet_2dmap.pdf
dmaps -i 2D_maps/upper_leaflet_2dmap.dat -o 2D_maps/upper_leaflet_2dmap.pdf

# Plot APL and Thicknes
profiler -p Fatslim -out thickness.pdf -d thickness_distribution.pdf -plot thickness -s
profiler -p Fatslim -out apl.pdf -d apl_distribution.pdf -plot apl -s

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
