# =================================================================================
# =================================================================================
# Script:"generate_mesh_s_h_1"
# Date: 2022-02-24
# Implemented by: Johannes Borgqvist and Carl Lundholm
# Description:
# The script generates a spherical FEM-mesh with
# one hole on the surface and an adjacent region.
# It uses the occ-part of Gmsh where numerous
# standard geometrical shapes can be implemented.
# This script is entirely based on code that
# Professor Christophe Geuzaine <cgeuzaine@uliege.be>
# (one of the creators of Gmsh) so kindly provided
# Johannes Borgqvist with in their e-mail exchange.
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import gmsh # For generating the meshes
import sys  # Needed by gmsh to launch the GUI, i.e. open the window with the plot
import math# For using mathematical functions
import numpy as np # We need the sine function to convert from cylindrical to geodesic radius
# =================================================================================
# =================================================================================
# Creating the mesh
# =================================================================================
# =================================================================================
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
# Part 3 out of 7: Add the hole
#---------------------------------------------------------------
# Cut with a cylinder to generate first spherical patch
# Add start and end points to the cylinder
c0 = (0,0,0) # Cylinder starting point
c1 = (0,0,-1) # First cylinder end point at the south pole
# Add the radius of the cylinder used to generate the hole
hole_radius = np.sin(0.3)
# Generate the hole by removing the intersection between the cylinder and the sphere
if hole_radius == 0:
    rest = [s0]
else:
    v1 = gmsh.model.occ.addCylinder(c0[0],c0[1],c0[2], c1[0],c1[1],c1[2], hole_radius)
    the_hole = gmsh.model.occ.intersect(s0, [(3, v1)], removeObject=False)
    rest = gmsh.model.occ.cut(s0, the_hole[0], removeObject=True)
#---------------------------------------------------------------
# Part 4 out of 7: Finalise the mesh
#---------------------------------------------------------------
# Fragment all surfaces to make everything conformal
gmsh.model.occ.synchronize()
# Set mesh size
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
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
