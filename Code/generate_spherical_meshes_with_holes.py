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
import numpy as np # For numpy functionalities
# =================================================================================
# =================================================================================
# Generate the meshes
# =================================================================================
# =================================================================================
# Create the hole radii in an automated fashion
hole_radii_list = [[r] for r in np.arange(0,6.5,0.5)]
# Allocate the lists of points defining the points on the sphere
c0_lists = []
c1_lists = []
# Generate a list of points defining the points on the sphere
for index in range(hole_radii_list):
    c0_lists.append([(0,0,0)])
    c1_lists.append([(-1,1,0)])
# Loop over all points and radii to generate the corresponding meshes
for list_index in range(len(hole_radii_list)):
    toolbox_generate_meshes.generate_spherical_mesh_with_holes(c0_lists[list_index],c1_lists[list_index],hole_radii_list[list_index])
