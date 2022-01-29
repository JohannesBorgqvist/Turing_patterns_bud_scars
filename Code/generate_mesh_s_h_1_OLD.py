# =================================================================================
# =================================================================================
# Script:"generate_mesh_s_h_1"
# Date: 2021-09-09
# Implemented by: Johannes Borgqvist
# Description:
# The script generates a spherical FEM-mesh with
# one hole on the surface and an adjacent region.
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
# =================================================================================
# =================================================================================
# Creating the mesh
# =================================================================================
# =================================================================================
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
#---------------------------------------------------------------
# Part 3: Add the hole
#---------------------------------------------------------------
# Cut with a cylinder to generate first spherical patch
v1 = gmsh.model.occ.addCylinder(0,0,0, 0,1,1, 0.25)
the_hole = gmsh.model.occ.intersect([(2, s0[0][1])], [(3, v1)], removeObject=False)
#hole_boundary = gmsh.model.getBoundary([the_hole[0][0]])
#---------------------------------------------------------------
# Part 4: Add the adjacent region
#---------------------------------------------------------------
# Create a wire in the parametric plane of the spherical surface
# [-pi,pi]x[-pi/2,pi/2] to create the other spherical patch
rest = gmsh.model.occ.addRectangle(math.pi/2-0.5,0.5,0, 1,0.6)
gmsh.model.occ.synchronize()
b3 = gmsh.model.getBoundary([(2, rest)])
w3 = gmsh.model.occ.addWire([p[1] for p in b3])
adjacent = gmsh.model.occ.addTrimmedSurface(s0[0][1], [w3])
#adj_curve = gmsh.model.getBoundary([(2, adjacent)])
gmsh.model.occ.remove([(2, rest)], recursive=True)
#---------------------------------------------------------------
# Part 4: Finalise the mesh
#---------------------------------------------------------------
# fragment all surfaces to make everything conformal
gmsh.model.occ.removeAllDuplicates()
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
# Add the adjacent region
adjacent_region = gmsh.model.addPhysicalGroup(2,[adjacent],2)
gmsh.model.setPhysicalName(2,adjacent_region,"Adjacent region")
# Add the hole
hole = gmsh.model.addPhysicalGroup(the_hole[0][0][0],[the_hole[0][0][1]],3)
gmsh.model.setPhysicalName(the_hole[0][0][0],hole,"Hole")
# -------------------------------------------------------------------
# Part 6: Add colours to the mesh
# -------------------------------------------------------------------
# For nice colours, please see (https://colorbrewer2.org/)...
#-----------------------------------------------------------------------------
# THE REST OF THE SPEHRE
gmsh.model.setColor([(2,rest)], 2, 56, 88)  # Darkest blue
# THE HOLE
gmsh.model.setColor([(2,the_hole[0][0][1])], 255, 247, 251)  # Light blue
# THE ADJACENT REGION
#gmsh.model.setColor([(1,circle_arc[3]) for circle_arc in list_of_lists[1]], 103, 0, 31)  # Darkest magenta
gmsh.model.setColor([(2,adjacent)], 116, 169, 207)  # In between blue
#---------------------------------------------------------------
# Part 6: Generate the mesh
#---------------------------------------------------------------
# Generate the two dimensional mesh (as we work with surfaces) 
gmsh.model.mesh.generate(2)
# Write the mesh to our file
gmsh.write("../Meshes/s_h_1.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()

