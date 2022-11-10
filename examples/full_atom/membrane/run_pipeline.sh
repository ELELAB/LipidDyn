# Import the environment prior to run this script

LipidDyn -f heterogeneous_membrane.xtc -s heterogeneous_membrane.gro -g full_atom.yml -n 4 -c -r -fatslim -mc -mov -op -2d 

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

# Membrane Curvature 
curvature -l1 curv/up_Avg_surface.dat -l2 curv/low_Avg_surface.dat -o curv/Avg_basic_surface.pdf -plot basic
curvature -l1 curv/up_Avg_surface.dat -l2 curv/low_Avg_surface.dat -o curv/Avg_smooth_surface.pdf -plot smooth
curvature -l1 curv/up_Avg_mean_curvature.dat -l2 curv/low_Avg_mean_curvature.dat -o curv/Avg_mean_curvature.pdf -plot mean
curvature -l1 curv/up_Avg_gaussian_curvature.dat -l2 curv/low_Avg_gaussian_curvature.dat -o curv/Avg_gaussian_curvature.pdf -plot gaussian
# Membrane Curvature - Avg. Plots side-by-side for each leaflet
curvature -l1 curv/low_Avg_surface.dat -l2 curv/low_Avg_mean_curvature.dat -l3 curv/low_Avg_gaussian_curvature.dat -o curv/3_Avg_lower_plots.pdf -plot 3_curvatures -t 'Lower Leaflet'
curvature -l1 curv/up_Avg_surface.dat -l2 curv/up_Avg_mean_curvature.dat -l3 curv/up_Avg_gaussian_curvature.dat -o curv/3_Avg_upper_plots.pdf -plot 3_curvatures -t 'Upper Leaflet'
# Membrane Curvature - Single and Multiframe Plots 
# (to be added)