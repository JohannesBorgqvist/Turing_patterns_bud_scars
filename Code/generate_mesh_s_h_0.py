# =================================================================================
# =================================================================================
# Script:"generate_mesh_s_h_0"
# Date: 2021-09-09
# Implemented by: Johannes Borgqvist
# Description:
# The script generates a spherical FEM-mesh.
# It uses the occ-part of Gmsh where numerous
# standard geometrical shapes can be implemented.
# This script is entirely based on code that
# Professor Christophe Geuzaine <cgeuzaine@uliege.be>
# (one of the creators of Gmsh) so kindly provided
# me with in our e-mail exchange.
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import gmsh # For generating the meshes
import sys  # Needed by gmsh to launch the GUI, i.e. open the window with the plot
import math# For using mathematical functions
#---------------------------------------------------------------
# Part 1: Initialise Gmsh
#---------------------------------------------------------------
gmsh.initialize(sys.argv)
#---------------------------------------------------------------
# Part 2: Add the larger unit sphere
#---------------------------------------------------------------
# Create sphere and get its boundary
v0 = gmsh.model.occ.addSphere(0,0,0, 1)
gmsh.model.occ.synchronize()
#s0 = gmsh.model.getBoundary([(3, v0)])[0][1]
s0 = gmsh.model.getBoundary([(3, v0)])
gmsh.model.occ.remove([(3, v0)])
# Define the rest of the sphere
rest = s0[0][1]
#---------------------------------------------------------------
# Part 4: Finalise the mesh
#---------------------------------------------------------------
gmsh.model.occ.synchronize()
# set mesh size
#gmsh.option.setNumber('Mesh.MeshSizeMax', 0.1)
#gmsh.option.setNumber('Mesh.MeshSizeMax', 0.075)
gmsh.option.setNumber('Mesh.MeshSizeMax', 0.055)
#---------------------------------------------------------------
# Part 5: Add physical regions
#---------------------------------------------------------------
# Add the rest of the sphere
rest_of_sphere = gmsh.model.addPhysicalGroup(2,[rest],1)
gmsh.model.setPhysicalName(2,rest_of_sphere,"Rest of the sphere")
# -------------------------------------------------------------------
# Part 6: Add colours to the mesh
# -------------------------------------------------------------------
# For nice colours, please see (https://colorbrewer2.org/)...
#-----------------------------------------------------------------------------
# THE REST OF THE SPEHRE
gmsh.model.setColor([(2,rest)], 2, 56, 88)  # Darkest blue
#---------------------------------------------------------------
# Part 6: Generate the mesh
#---------------------------------------------------------------
# Generate the two dimensional mesh (as we work with surfaces) 
gmsh.model.mesh.generate(2)
# Write the mesh to our file
gmsh.write("../Meshes/s_h_0.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
