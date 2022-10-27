# Import the environment prior to run this script

LipidDyn -f heterogeneous_membrane_protein.xtc -s heterogeneous_membrane_protein.gro -g coarse_grained.yml -p -a -mc -n 4 -c -r

# Plotting 

# Plot density maps
dmaps -i 2D_maps/lower_leaflet_2dmap.dat -o 2D_maps/lower_leaflet_2dmap.pdf
dmaps -i 2D_maps/upper_leaflet_2dmap.dat -o 2D_maps/upper_leaflet_2dmap.pdf

#Enrichment
dmaps -i Enrichment/CHOL_lower_leaflet_enrich.dat -o Enrichment/CHOL_lower_leaflet_enrich.pdf -enr
dmaps -i Enrichment/CHOL_upper_leaflet_enrich.dat -o Enrichment/CHOL_upper_leaflet_enrich.pdf -enr
dmaps -i Enrichment/DPSM_lower_leaflet_enrich.dat -o Enrichment/DPSM_lower_leaflet_enrich.pdf -enr
dmaps -i Enrichment/DPSM_upper_leaflet_enrich.dat -o Enrichment/DPSM_upper_leaflet_enrich.pdf -enr
dmaps -i Enrichment/POPC_lower_leaflet_enrich.dat -o Enrichment/POPC_lower_leaflet_enrich.pdf -enr
dmaps -i Enrichment/POPC_upper_leaflet_enrich.dat -o Enrichment/POPC_upper_leaflet_enrich.pdf -enr
  

# Plot APL and Thicknes
profiler -i Fatslim/thickness.xvg -out Fatslim/thickness.pdf -d Fatslim/distribution_thickness.pdf -plot thick
profiler -i Fatslim/apl.xvg -out Fatslim/apl.pdf -d Fatslim/distribution_apl.pdf -plot apl

# Diffusion
diffusion -i Diffusion_movements/Lower_leaflet_coordinates.dat -o Diffusion_movements/Lower_leaflet_coordinates -t he 
diffusion -i Diffusion_movements/Lower_leaflet_coordinates.dat -o Diffusion_movements/Lower_leaflet_coordinates -t he -m 
diffusion -i Diffusion_movements/Upper_leaflet_coordinates.dat -o Diffusion_movements/Upper_leaflet_coordinates -t he -m
diffusion -i Diffusion_movements/Upper_leaflet_coordinates.dat -o Diffusion_movements/Upper_leaflet_coordinates -t he

# Order Parameter
ordpar -i Order_Parameter/Order_Parameter_DPSM.csv -o Order_Parameter/Order_Parameter_DPSM.pdf -s scc
ordpar -i Order_Parameter/Order_Parameter_POPC.csv -o Order_Parameter/Order_Parameter_POPC.pdf -s scc

# Membrane Curvature
curvature -u curv/up_Avg_surface.dat -l curv/low_Avg_surface.dat -o curv/Avg_basic_surface.pdf -plot basic
curvature -u curv/up_Avg_surface.dat -l curv/low_Avg_surface.dat -o curv/Avg_smooth_surface.pdf -plot smooth
curvature -u curv/up_Avg_mean_curvature.dat -l curv/low_Avg_mean_curvature.dat -o curv/Avg_mean_curvature.pdf -plot mean
curvature -u curv/up_Avg_gaussian_curvature.dat -l curv/low_Avg_gaussian_curvature.dat -o curv/Avg_gaussian_curvature.pdf -plot gaussian
