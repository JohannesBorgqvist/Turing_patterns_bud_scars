# =================================================================================
# =================================================================================
# Script:"toolbox_generate_meshes"
# Date: 2022-03-18
# Implemented by: Johannes Borgqvist and Carl Lundholm
# Description:
# This script contains the sole function which generates the spherical meshes with a single
# hole in them. The holes are generated by jamming a cylinder through the sphere and then
# remove the intersection between the two surfaces.
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import gmsh # For generating the meshes
import sys  # Needed by gmsh to launch the GUI, i.e. open the window with the plot
import math# For using mathematical functions
# =================================================================================
# =================================================================================
# Functions for generating the meshes
# =================================================================================
# =================================================================================
#------------------------------------------------------------------
# Function 1: "generate_spherical_mesh_with_holes"
# The function takes the following inputs:
# 1. The list of start points of the cylinders c0_list,
# 2. The list of end points of the cylinders c1_list,
# 3. The list of radii of the cylinder hole_radii,
# The function saves a mesh in the meshes folder named after
# the number of holes as well as their respective radii.
#------------------------------------------------------------------
def generate_spherical_mesh_with_holes(c0_list,c1_list,hole_radii):
    #---------------------------------------------------------------
    # Part 1 out of 7: Initialise Gmsh
    #---------------------------------------------------------------
    gmsh.initialize(sys.argv)
    #---------------------------------------------------------------
    # Part 2 out of 7: Add the larger unit sphere
    #---------------------------------------------------------------
    # Create sphere and get its boundary
    v0 = gmsh.model.occ.addSphere(0,0,0, 1)
    gmsh.model.occ.synchronize()
    s0 = gmsh.model.getBoundary([(3, v0)])
    gmsh.model.occ.remove([(3, v0)])
    #---------------------------------------------------------------
    # Part 3 out of 7: Add the holes
    #---------------------------------------------------------------
    # Initialise the rest of the sphere
    rest = [s0]
    # Loop of the holes and add them
    for hole_index in range(len(hole_radii)):
        # Extract the points on the cylindes
        c0 = c0_list[hole_index]
        c1 = c1_list[hole_index]
        # Extract the hole radius
        hole_radius = hole_radii[hole_index]
        if hole_radius > 0:
            # Create the cylinder
            v1 = gmsh.model.occ.addCylinder(c0[0],c0[1],c0[2], c1[0],c1[1],c1[2], hole_radius)
            # Create the hole being the intersection between the unit sphere and the cylinder
            the_hole = gmsh.model.occ.intersect(rest[0], [(3, v1)], removeObject=False)
            # Remove the hole
            rest = gmsh.model.occ.cut(rest[0], the_hole[0], removeObject=True)    
    #---------------------------------------------------------------
    # Part 4 out of 7: Finalise the mesh
    #---------------------------------------------------------------
    # fragment all surfaces to make everything conformal
    gmsh.model.occ.synchronize()
    # set mesh size
    gmsh.option.setNumber('Mesh.MeshSizeMax', 0.055)
    #---------------------------------------------------------------
    # Part 5 out of 7: Add physical regions
    #---------------------------------------------------------------
    # Add the rest of the sphere
    rest_of_sphere = gmsh.model.addPhysicalGroup(2,[rest[0][0][1]],1)
    gmsh.model.setPhysicalName(2,rest_of_sphere,"The sphere with a hole with radius, r=" + str(hole_radius))
    # -------------------------------------------------------------------
    # Part 6 out of 7: Add colours to the mesh
    # -------------------------------------------------------------------
    # For nice colours, please see (https://colorbrewer2.org/)...
    #-----------------------------------------------------------------------------
    # THE REST OF THE SPHERE
    gmsh.model.setColor([(2,rest[0][0][1])], 2, 56, 88)  # Darkest blue
    #---------------------------------------------------------------
    # Part 7 out of 7: Generate the mesh
    #---------------------------------------------------------------
    # Generate the two dimensional mesh (as we work with surfaces) 
    gmsh.model.mesh.generate(2)
    # Allocate a file name for the mesh
    if len(hole_radii) == 1 and hole_radii[0]== 0:
        mesh_name = "../Meshes/s_h_0.msh"
    else:
        mesh_name = "../Meshes/s_h_" + str(len(hole_radii)) + "_"
        # Loop over the hole radii and add these to the file names
        for hole_radius in hole_radii:
            mesh_name += "r_" + str(round(hole_radius,3)).replace(".","p") + "_"
        # Add the suffix
        mesh_name += ".msh"
    # Write the mesh to our file
    gmsh.write(mesh_name.replace("_.msh",".msh"))    
    # Shut down gmsh
    gmsh.finalize()
