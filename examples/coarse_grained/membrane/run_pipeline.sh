# Import the environment prior to run this script

LipidDyn -f heterogeneous_membrane.xtc -s heterogeneous_membrane.gro -g coarse_grained.yml -fatslim -op -2d -mc -mov -enr -n 4 -c -spe


# Plotting 

# Plot density maps
dmaps -i 2D_maps/lower_leaflet_2dmap.dat -o 2D_maps/lower_leaflet_2dmap.pdf
dmaps -i 2D_maps/upper_leaflet_2dmap.dat -o 2D_maps/upper_leaflet_2dmap.pdf

# Plot APL and Thicknes
profiler -p Fatslim -out thickness.pdf -d thickness_distribution.pdf -plot thickness -spe
profiler -p Fatslim -out apl.pdf -d apl_distribution.pdf -plot apl -spe

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

# Membrane Curvature - Avg. Plots
curvature -l1 curv/up_Avg_surface.dat -l2 curv/low_Avg_surface.dat -o curv/Avg_basic_surface.pdf -plot basic
curvature -l1 curv/up_Avg_surface.dat -l2 curv/low_Avg_surface.dat -o curv/Avg_smooth_surface.pdf -plot smooth 
curvature -l1 curv/up_Avg_mean_curvature.dat -l2 curv/low_Avg_mean_curvature.dat -o curv/Avg_mean_curvature.pdf -plot mean 
curvature -l1 curv/up_Avg_gaussian_curvature.dat -l2 curv/low_Avg_gaussian_curvature.dat -o curv/Avg_gaussian_curvature.pdf -plot gaussian 
# Membrane Curvature - Avg. Plots side-by-side for each leaflet
curvature -l1 curv/low_Avg_surface.dat -l2 curv/low_Avg_mean_curvature.dat -l3 curv/low_Avg_gaussian_curvature.dat -o curv/3_Avg_lower_plots.pdf -plot 3_curvatures -t 'Lower Leaflet' -level 60
curvature -l1 curv/up_Avg_surface.dat -l2 curv/up_Avg_mean_curvature.dat -l3 curv/up_Avg_gaussian_curvature.dat -o curv/3_Avg_upper_plots.pdf -plot 3_curvatures -t 'Upper Leaflet' -level 60
# Membrane Curvature - Single Frame Plots
curvature -l1 curv/up_AF_surface.dat -l2 curv/low_AF_surface.dat -o curv/AF_basic_surface.pdf -plot basic -frame 5
curvature -l1 curv/up_AF_surface.dat -l2 curv/low_AF_surface.dat -o curv/AF_smooth_surface.pdf -plot smooth -frame 5
curvature -l1 curv/up_AF_gaussian_curvature.dat -l2 curv/low_AF_gaussian_curvature.dat -o curv/AF_gaussian_curvature.pdf -plot gaussian -frame 5 -level 40
curvature -l1 curv/up_AF_mean_curvature.dat -l2 curv/low_AF_mean_curvature.dat -o curv/AF_mean_curvature.pdf -plot mean -frame 5 -level 40
# Membrane Curvature - Multi Frame (6 frame) Plots
curvature -l1 curv/up_AF_surface.dat -l2 curv/low_AF_surface.dat -o curv/6f_up_smooth_surface.pdf -o2 curv/6f_low_smooth_surface.pdf -plot first_6
curvature -l1 curv/up_AF_surface.dat -l2 curv/low_AF_surface.dat -o curv/6m_up_smooth_surface.pdf -o2 curv/6m_low_smooth_surface.pdf -plot middle_6
curvature -l1 curv/up_AF_surface.dat -l2 curv/low_AF_surface.dat -o curv/6l_up_smooth_surface.pdf -o2 curv/6l_low_smooth_surface.pdf -plot last_6
