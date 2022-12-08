# Import the environment prior to run this script

LipidDyn -f heterogeneous_membrane.xtc -s heterogeneous_membrane.gro -g coarse_grained.yml -fatslim -op -2d -mov -enr -n 4 -c -spe

# Plotting 

# Plot density maps
dmaps -i 2D_maps/lower_leaflet_2dmap.dat -o 2D_maps/lower_leaflet_2dmap.pdf
dmaps -i 2D_maps/upper_leaflet_2dmap.dat -o 2D_maps/upper_leaflet_2dmap.pdf

# Plot APL and Thicknes
profiler -p Fatslim -out thickness.pdf -d thickness_distribution.pdf -plot thickness -s
profiler -p Fatslim -out apl.pdf -d apl_distribution.pdf -plot apl -s

# Diffusion
diffusion -i Diffusion_movements/Lower_leaflet_coordinates.dat -o Diffusion_movements/Lower_leaflet_coordinates -t he -m
diffusion -i Diffusion_movements/Lower_leaflet_coordinates.dat -o Diffusion_movements/Lower_leaflet_coordinates -t he 
diffusion -i Diffusion_movements/Upper_leaflet_coordinates.dat -o Diffusion_movements/Upper_leaflet_coordinates -t he -m 
diffusion -i Diffusion_movements/Upper_leaflet_coordinates.dat -o Diffusion_movements/Upper_leaflet_coordinates -t he 

# Order Parameter
ordpar -i Order_Parameter/Order_Parameter_DOPA.csv -o Order_Parameter/Order_Parameter_DOPA.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_DOPE.csv -o Order_Parameter/Order_Parameter_DOPE.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_DVPC.csv -o Order_Parameter/Order_Parameter_DVPC.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_POPC.csv -o Order_Parameter/Order_Parameter_POPC.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_PAPI.csv -o Order_Parameter/Order_Parameter_PAPI.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_PIPI.csv -o Order_Parameter/Order_Parameter_PIPI.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_POPI.csv -o Order_Parameter/Order_Parameter_POPI.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_DOPC.csv -o Order_Parameter/Order_Parameter_DOPC.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_DPSM.csv -o Order_Parameter/Order_Parameter_DPSM.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_PAPE.csv -o Order_Parameter/Order_Parameter_PAPE.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_PIPE.csv -o Order_Parameter/Order_Parameter_PIPE.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_PNSM.csv -o Order_Parameter/Order_Parameter_PNSM.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_POPG.csv -o Order_Parameter/Order_Parameter_POPG.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_POPS.csv -o Order_Parameter/Order_Parameter_POPS.pdf -s scc
   
