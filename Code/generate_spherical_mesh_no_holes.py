# =================================================================================
# =================================================================================
# Script:"generate_spherical_mesh"
# Date: 2021-08-25
# Implemented by: Johannes Borgqvist
# Description:
# The script generates a spherical FEM-mesh with
# one hole on the surface. It uses the different
# functions that are stored in the help script
# "functions_generate_sphere_with_holes.py".
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import gmsh  # For generating the meshes
import math  # For using mathematical functions
import sys  # Needed by gmsh to launch the GUI, i.e. open the window with the plot
import numpy as np  # Needed to create linspace arrays and other nice functions
import functions_generate_spherical_mesh_with_holes  # Home-made

# =================================================================================
# =================================================================================
# Initialising Gmsh
# =================================================================================
# =================================================================================
# Initialise the Gmsh framework
gmsh.initialize()
# Add a model entity
gmsh.model.add("sphere_with_holes")
# We start by defining some points and some lines. To make the code shorter we can redefine a namespace:
hole_on_sphere = gmsh.model.geo
# =================================================================================
# =================================================================================
# Assembling the mesh
# =================================================================================
# =================================================================================
# -------------------------------------------------------------------
# Allocating memory for all lists and vectors
# -------------------------------------------------------------------
# Define a centimetre
cm = 1e-02
Lc1 = 0.01
R = 1 # Big radius
#--------------------------------------------------------------------
# ADD THE MAIN CENTRE
#--------------------------------------------------------------------
c1 = (0,0,0) # The main centre of the unit ball
label_counter = 1# Define a label counter
main_centre = (c1[0],c1[1],c1[2],Lc1,label_counter)
hole_on_sphere.addPoint(main_centre[0], main_centre[1], main_centre[2], main_centre[3], main_centre[4])
label_counter += 1
#--------------------------------------------------------------------
# DEFINE THE TWO POLES
#--------------------------------------------------------------------
# Define the north pole which will be our start point every time
north_pole = (0,0,R,Lc1,label_counter)
# Add the point to the mesh
hole_on_sphere.addPoint(north_pole[0], north_pole[1], north_pole[2], north_pole[3], north_pole[4])        
# Increment the point counter
label_counter += 1
# Define the south pole which will be our start point every time
south_pole = (0,0,-R,Lc1,label_counter)
# Add the point to the mesh
hole_on_sphere.addPoint(south_pole[0], south_pole[1], south_pole[2],south_pole[3], south_pole[4])        
# Increment the point counter
label_counter += 1
#---------------------------------------------------------------
# Calculate the rest of the sphere
#---------------------------------------------------------------    
# SEGMENT BELOW HOLE
# Define phi
phi_below = (0,math.pi)
# Define theta
theta_below = [0,2*math.pi]
# Caclulate the rest of the sphere
new_counter, list_of_lists_sphere = functions_generate_spherical_mesh_with_holes.add_sphere_without_hole(R,phi_below,theta_below,label_counter,north_pole,south_pole,main_centre)
# Update the label_counter
label_counter = new_counter+1
#---------------------------------------------------------------
# Add the sphere
#---------------------------------------------------------------
# Add all points
for segment in list_of_lists_sphere[0]:
    for point in segment:
        if point[4]!= 2 and point[4]!= 3:
            hole_on_sphere.addPoint(point[0], point[1], point[2], point[3], point[4])
# Add all lines
for segment in list_of_lists_sphere[1]:
    for circle_arc in segment:
        hole_on_sphere.addCircleArc(circle_arc[0], circle_arc[2], circle_arc[1], circle_arc[3])
# Add all curve loops
for segment in list_of_lists_sphere[2]:
    for curve_loop in segment:
        hole_on_sphere.addCurveLoop(curve_loop[0],curve_loop[1])
# Add all surfaces
for segment in list_of_lists_sphere[3]:
    for surface in segment:
        hole_on_sphere.addPlaneSurface(surface[0],surface[1])
# =================================================================================
# =================================================================================
# Finishing up: Writing the mesh to a file and marking physical groups
# =================================================================================
# =================================================================================
# -------------------------------------------------------------------
# SYNCHRONIZE THE MESH
# -------------------------------------------------------------------
# Synchronise everything with Gmsh
hole_on_sphere.synchronize()
# -------------------------------------------------------------------
# DEFINE PHYSICAL REGIONS
# -------------------------------------------------------------------
surface_rest_of_sphere = gmsh.model.addPhysicalGroup(2,[list_of_lists_sphere[3][i][0][1] for i in range(len(list_of_lists_sphere[3]))])
gmsh.model.setPhysicalName(2, surface_rest_of_sphere, "Surface the rest of the sphere")
# -------------------------------------------------------------------
# DEFINE COLOURS FOR THE MESH
# -------------------------------------------------------------------
# For nice colours, please see (https://colorbrewer2.org/)...
#-----------------------------------------------------------------------------
# THE REST OF THE SPEHRE
gmsh.model.setColor([(2,list_of_lists_sphere[3][i][0][1]) for i in range(len(list_of_lists_sphere[3]))], 2, 56, 88)  # Darkest blue
# -------------------------------------------------------------------
# GENERATING THE MESH
# -------------------------------------------------------------------
# Generate the two dimensional mesh (as we work with surfaces) 
gmsh.model.mesh.generate(2)
# Write the mesh to our file
gmsh.write("../Meshes/sphere_with_no_holes.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
