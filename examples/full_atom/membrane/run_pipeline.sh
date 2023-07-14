# Import the environment prior to run this script

LipidDyn -f heterogeneous_membrane.xtc -s heterogeneous_membrane.gro -g full_atom.yml -n 4 -c -fatslim -mc -mov -op -2d -spe

# Plotting 

# Plot density maps
dmaps -i 2D_maps/lower_leaflet_2dmap.dat -o 2D_maps/lower_leaflet_2dmap.pdf
dmaps -i 2D_maps/upper_leaflet_2dmap.dat -o 2D_maps/upper_leaflet_2dmap.pdf

# Plot APL and Thicknes
profiler -p Fatslim -out thickness.pdf -d thickness_distribution.pdf -plot thickness -spe
profiler -p Fatslim -out apl.pdf -d apl_distribution.pdf -plot apl -spe

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
# Membrane Curvature - Single Frame Plots
curvature -l1 curv/up_AF_surface.dat -l2 curv/low_AF_surface.dat -o curv/AF_basic_surface.pdf -plot basic -frame 5
curvature -l1 curv/up_AF_surface.dat -l2 curv/low_AF_surface.dat -o curv/AF_smooth_surface.pdf -plot smooth -frame 5
curvature -l1 curv/up_AF_gaussian_curvature.dat -l2 curv/low_AF_gaussian_curvature.dat -o curv/AF_gaussian_curvature.pdf -plot gaussian -frame 5 -level 40
curvature -l1 curv/up_AF_mean_curvature.dat -l2 curv/low_AF_mean_curvature.dat -o curv/AF_mean_curvature.pdf -plot mean -frame 5 -level 40
# Membrane Curvature - Multi Frame (6 frame) Plots
curvature -l1 curv/up_AF_surface.dat -l2 curv/low_AF_surface.dat -o curv/6f_up_smooth_surface.pdf -o2 curv/6f_low_smooth_surface.pdf -plot first_6
curvature -l1 curv/up_AF_surface.dat -l2 curv/low_AF_surface.dat -o curv/6m_up_smooth_surface.pdf -o2 curv/6m_low_smooth_surface.pdf -plot middle_6
curvature -l1 curv/up_AF_surface.dat -l2 curv/low_AF_surface.dat -o curv/6l_up_smooth_surface.pdf -o2 curv/6l_low_smooth_surface.pdf -plot last_6

# Scrambling plot (if scrambling lipids are detected)
scrambling -i Scrambling/lipids/ -o Scrambling/z_trajectory_all_relative.pdf
scrambling -i Scrambling/lipids/ -o Scrambling/z_trajectory_all.pdf -a