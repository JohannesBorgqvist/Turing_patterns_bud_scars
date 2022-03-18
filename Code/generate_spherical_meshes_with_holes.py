# =================================================================================
# =================================================================================
# Script:"generate_spherical_meshes_with_holes"
# Date: 2022-03-18
# Implemented by: Johannes Borgqvist and Carl Lundholm
# Description:
# This is the script where we automate the generation of the spherical meshes with
# holes in them.  This function calls the only function that is stored in the
# script called "toolbox_generate_meshes.py". 
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import toolbox_generate_meshes # Home-made
# =================================================================================
# =================================================================================
# Generate the meshes
# =================================================================================
# =================================================================================
# Generate a list of points defining the points on the sphere
c0_lists = [[(0,0,0)], [(0,0,0)], [(0,0,0)], [(0,0,0)], [(0,0,0)]]
c1_lists = [[(-1,1,0)], [(-1,1,0)], [(-1,1,0)], [(-1,1,0)], [(-1,1,0)]]
# Generate a list of hole radii
hole_radii_list = [[0.0],[0.3],[0.4],[0.5],[0.6]]
# Loop over all points and radii to generate the corresponding meshes
for list_index in range(len(hole_radii_list)):
    toolbox_generate_meshes.generate_spherical_mesh_with_holes(c0_lists[list_index],c1_lists[list_index],hole_radii_list[list_index])
