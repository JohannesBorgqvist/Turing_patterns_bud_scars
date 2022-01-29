# =================================================================================
# =================================================================================
# Script:"generate_mesh_s_h_2"
# Date: 2022-01-28
# Implemented by: Johannes Borgqvist and Carl Lundholm
# Description:
# The script generates a spherical FEM-mesh with
# two holes on the surface and an adjacent region for
# each hole. It uses the occ-part of Gmsh where numerous
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
# Part 3: Add the two holes
#---------------------------------------------------------------
# Cut with cylinders to generate spherical patches
c0 = (0,0,0) # Cylinder starting point (use for both holes)
hole_radius = 0.25
# HOLE 1
#c1 = (0,1,1) # Original first cylinder end point
c1 = (-1,1,0) # First cylinder end point
vh1 = gmsh.model.occ.addCylinder(c0[0],c0[1],c0[2], c1[0],c1[1],c1[2], hole_radius)
hole_1 = gmsh.model.occ.intersect(s0, [(3, vh1)], removeObject=False)
rest = gmsh.model.occ.cut(s0, hole_1[0], removeObject=True)
# HOLE 2
#c2 = (-0.75,-0.15,1.7) # Original second cylinder end point
c2 = (-1,-1,0) # Second cylinder end point
vh2 = gmsh.model.occ.addCylinder(c0[0],c0[1],c0[2], c2[0],c2[1],c2[2], hole_radius)
hole_2 = gmsh.model.occ.intersect(rest[0], [(3, vh2)], removeObject=False)
rest = gmsh.model.occ.cut(rest[0], hole_2[0], removeObject=False)
#---------------------------------------------------------------
# Part 4: Add the adjacent regions
#---------------------------------------------------------------
# Cut with cylinders to generate spherical patches
#adjacent_region_thickness = 0.5*hole_radius
adjacent_region_thickness = 0.125
adjacent_region_radius = hole_radius + adjacent_region_thickness
# ADJACENT REGION 1
va1 = gmsh.model.occ.addCylinder(c0[0],c0[1],c0[2], c1[0],c1[1],c1[2], adjacent_region_radius)
adjacent_region_1 = gmsh.model.occ.intersect(rest[0], [(3, va1)], removeObject=False)
# ADJACENT REGION 2
va2 = gmsh.model.occ.addCylinder(c0[0],c0[1],c0[2], c2[0],c2[1],c2[2], adjacent_region_radius)
adjacent_region_2 = gmsh.model.occ.intersect(rest[0], [(3, va2)], removeObject=False)
### OLD WEIRD RECTANGULAR REGIONS ####
# ADJACENT REGION 1
# Create a wire in the parametric plane of the spherical surface
# [-pi,pi]x[-pi/2,pi/2] to create the other spherical patch
#rest_1 = gmsh.model.occ.addRectangle(math.pi/2-0.5,0.5,0, 1,0.6)
#gmsh.model.occ.synchronize()
#b3_1 = gmsh.model.getBoundary([(2, rest_1)])
#w3_1 = gmsh.model.occ.addWire([p[1] for p in b3_1])
#adjacent_1 = gmsh.model.occ.addTrimmedSurface(s0[0][1], [w3_1])
# ADJACENT REGION 2
# Create a wire in the parametric plane of the spherical surface
# [-pi,pi]x[-pi/2,pi/2] to create the other spherical patch
#rest_2 = gmsh.model.occ.addRectangle(-9*math.pi/8-0.12,0.82,1, 1.5,0.6)
#gmsh.model.occ.synchronize()
#b3_2 = gmsh.model.getBoundary([(2, rest_2)])
#w3_2 = gmsh.model.occ.addWire([p[1] for p in b3_2])
#adjacent_2 = gmsh.model.occ.addTrimmedSurface(s0[0][1], [w3_2])
# Update the adjacent region
#adjacent_1 = rest_2
# Remove the squares that are projected onto the circle
#gmsh.model.occ.remove([(2, rest_1)], recursive=True)
#gmsh.model.occ.remove([(2, rest_2)], recursive=True)
#---------------------------------------------------------------
# Part 5: "Add" the rest of the sphere
#---------------------------------------------------------------
# Weird trick from old code (See Part 4) to get rest of sphere 
rest_weird = gmsh.model.occ.addRectangle(math.pi/2-0.5,0.5,0, 1,0.6)
gmsh.model.occ.remove([(2, rest_weird)], recursive=True)
#---------------------------------------------------------------
# Part 6: Finalise the mesh
#---------------------------------------------------------------
# fragment all surfaces to make everything conformal
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()
# set mesh size
#gmsh.option.setNumber('Mesh.MeshSizeMax', 0.1)
#gmsh.option.setNumber('Mesh.MeshSizeMax', 0.075)
gmsh.option.setNumber('Mesh.MeshSizeMax', 0.055)
#---------------------------------------------------------------
# Part 7: Add physical regions
#---------------------------------------------------------------
# Add the rest of the sphere
#rest_of_sphere = gmsh.model.addPhysicalGroup(2,[rest_1, rest_2],1)
rest_of_sphere = gmsh.model.addPhysicalGroup(2,[rest_weird],1)
gmsh.model.setPhysicalName(2,rest_of_sphere,"Rest of the sphere")
# Add the adjacent regions
adjacent_regions = gmsh.model.addPhysicalGroup(2,[adjacent_region_1[0][0][1], adjacent_region_2[0][0][1]],2)
gmsh.model.setPhysicalName(2,adjacent_regions,"Adjacent regions")
# JUST LUMP IN ALL HOLES INTO THE SAME DOMAIN
# Add holes
#hole = gmsh.model.addPhysicalGroup(hole_1[0][0][0],[hole_1[0][0][1], hole_2[0][0][1]],3)
#gmsh.model.setPhysicalName(hole_1[0][0][0],hole,"Hole")
# Add hole 1 separately
#hole = gmsh.model.addPhysicalGroup(hole_1[0][0][0],[hole_1[0][0][1]],3)
#gmsh.model.setPhysicalName(hole_1[0][0][0],hole,"Hole 1")
# Add hole 2 separately
#hole = gmsh.model.addPhysicalGroup(hole_2[0][0][0],[hole_2[0][0][1]],4)
#gmsh.model.setPhysicalName(hole_1[0][0][0],hole,"Hole 2")
# -------------------------------------------------------------------
# Part 8: Add colours to the mesh
# -------------------------------------------------------------------
# For nice colours, please see (https://colorbrewer2.org/)...
#-----------------------------------------------------------------------------
# THE REST OF THE SPHERE
gmsh.model.setColor([(2,rest_weird)], 2, 56, 88)  # Darkest blue
# THE HOLE
#gmsh.model.setColor([(2,hole_1[0][0][1]),(2,hole_2[0][0][1])], 255, 247, 251)  # Light blue
# THE ADJACENT REGIONS
gmsh.model.setColor([(2,adjacent_region_1[0][0][1]),(2,adjacent_region_2[0][0][1])], 116, 169, 207) # In between blue
#---------------------------------------------------------------
# Part 9: Generate the mesh
#---------------------------------------------------------------
# Generate the two dimensional mesh (as we work with surfaces) 
gmsh.model.mesh.generate(2)
# Write the mesh to our file
gmsh.write("../Meshes/s_h_2.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()


