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
r = 0.1*R # Small radius
# -------------------------------------------------------------------
# Define holes
# -------------------------------------------------------------------
# SPHERE 1: Centre
c1 = (0,0,0) # with radius R defined above
# SPHERE 1: Centre
x2 = 0.1
y2 = -0.2
z2 = math.sqrt(R**2-x2**2-y2**2)
c2 = (x2,y2,z2) # with radius R defined above
# Add all centres to a list
list_of_centres = [c1, c2]
# Define a label counter
label_counter = 1
# Define the step size
step_size = math.pi/20
#step_size = math.pi/100
#--------------------------------------------------------------------
# ADD THE MAIN CENTRE
#--------------------------------------------------------------------
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
# Compute all holes in the mesh
#---------------------------------------------------------------
# Calculate holes
new_label, list_of_lists, theta_hole, phi_hole = functions_generate_spherical_mesh_with_holes.add_circular_holes_on_sphere(list_of_centres,r,R,label_counter,main_centre)
# Update the label counter
label_counter = new_label+1
#---------------------------------------------------------------
# Add the hole to the mesh
#---------------------------------------------------------------
# Add all points
for point in list_of_lists[0]:
    hole_on_sphere.addPoint(point[0], point[1], point[2], point[3], point[4])
# Add all circle arcs    
for circle_arc in list_of_lists[1]:
    hole_on_sphere.addCircleArc(circle_arc[0], circle_arc[1], circle_arc[2], circle_arc[3])
# Add all circle arcs    
for circle_arc in list_of_lists[2]:
    hole_on_sphere.addCircleArc(circle_arc[0], circle_arc[1], circle_arc[2], circle_arc[3])    
# Add all curve loops adjacent region
for curve_loop in list_of_lists[3]:
    hole_on_sphere.addCurveLoop(curve_loop[0], curve_loop[1])    
# Add all surfaces for the hole
for hole in list_of_lists[4]:
    hole_on_sphere.addPlaneSurface(hole[0], hole[1])
# Add all surfaces for the hole
for adjacent_hole in list_of_lists[5]:
    hole_on_sphere.addPlaneSurface(adjacent_hole[0], adjacent_hole[1])
#---------------------------------------------------------------
# Calculate the rest of the sphere
#---------------------------------------------------------------    
# SEGMENT BELOW HOLE
# Define phi
#phi_below = (phi_hole[1],math.pi)
epsilon = 0.1
phi_below = (0,math.pi)
#phi_below = (epsilon,math.pi)
#phi_below = (0,math.pi-epsilon)
#phi_below = (epsilon,math.pi-epsilon)
# Define theta
#theta_below = theta_hole
theta_below = [theta_hole[0], theta_hole[0]+math.pi/24]
# Caclulate the rest of the sphere
new_counter, list_of_lists_sphere_below = functions_generate_spherical_mesh_with_holes.add_sphere_without_hole(R,phi_below,theta_below,label_counter,north_pole,south_pole,main_centre)
# Update the label_counter
label_counter = new_counter+1
# SEGMENT ABOVE HOLE
# Define phi
phi_above = (0,phi_hole[0])
# Define theta
theta_above = theta_hole
# Caclulate the rest of the sphere
new_counter, list_of_lists_sphere_above = functions_generate_spherical_mesh_with_holes.add_sphere_without_hole(R,phi_above,theta_above,label_counter,north_pole,south_pole,main_centre)
# Update the label_counter
label_counter = new_counter+1
# SEGMENT ON THE LEFT
# Define phi
phi_left = (0,math.pi)
# Define theta
theta_left = [-np.pi, theta_hole[0]]
# Caclulate the rest of the sphere
new_counter, list_of_lists_sphere_left = functions_generate_spherical_mesh_with_holes.add_sphere_without_hole(R,phi_left,theta_left,label_counter,north_pole,south_pole,main_centre)
# Update the label_counter
label_counter = new_counter+1
# SEGMENT ON THE RIGHT
# Define phi
phi_right = (0,math.pi)
# Define theta
theta_right = [theta_hole[1],np.pi]
# Caclulate the rest of the sphere
new_counter, list_of_lists_sphere_right = functions_generate_spherical_mesh_with_holes.add_sphere_without_hole(R,phi_right,theta_right,label_counter,north_pole,south_pole,main_centre)
# Update the label_counter
label_counter = new_counter+1
#---------------------------------------------------------------
# MERGE THE LISTS:
#---------------------------------------------------------------
# Allocate memory for the final list which we will plot
list_of_lists_sphere = []
# Save the number of fragments of the spehre which
# we loop over. This number is four where we add
# the surfaces above the hole, below the hole, to the
# left of the hole and to the right of the hole.
number_of_sides = 4 # Above, below, right and left
# Allocate four lists the four types of objects
# that we have namely points, circle arcs, curve loops
# and surfaces.
points_rest_of_sphere = []
circle_arcs_rest_of_sphere = []
curve_loops_rest_of_sphere = []
surfaces_rest_of_sphere = []
# Loop over our four sides and add all objects into the
# respective lists
for side in range(0,0):
#for side in range(number_of_sides):
#for side in range(0,1):#Above
#for side in range(1,2):#Below
#for side in range(2,3):# Left
#for side in range(3,4):#Right    
    # Allocate memory for a data list
    data_list = []
    # The number of segments in each of our four sides
    # determines how long the subsequent nested for-loop
    # is
    if side == 0:# Above
        number_of_segments = len(list_of_lists_sphere_above[0])
        data_list = list_of_lists_sphere_above
    elif side == 1: # Below
        number_of_segments = len(list_of_lists_sphere_below[0])
        data_list = list_of_lists_sphere_below        
    elif side == 2: # Left
        number_of_segments = len(list_of_lists_sphere_left[0])
        data_list = list_of_lists_sphere_left        
    elif side == 3: # Right
        number_of_segments = len(list_of_lists_sphere_right[0])
        data_list = list_of_lists_sphere_right        
    # Loop over the number of segments and add them all to our four
    # separate lists
    for segment_index in range(number_of_segments):
        # Add the points
        points_rest_of_sphere.append(data_list[0][segment_index])
        # Add the points
        circle_arcs_rest_of_sphere.append(data_list[1][segment_index])
        # Add the points
        curve_loops_rest_of_sphere.append(data_list[2][segment_index])
        # Add the points
        surfaces_rest_of_sphere.append(data_list[3][segment_index])        
#---------------------------------------------------------------
# Add the rest of the sphere
#---------------------------------------------------------------
# Add all points
for segment in points_rest_of_sphere:
    for point in segment:
        if point[4]!= 2 and point[4]!= 3:
            hole_on_sphere.addPoint(point[0], point[1], point[2], point[3], point[4])
# Add all lines
for segment in circle_arcs_rest_of_sphere:
    for circle_arc in segment:
        hole_on_sphere.addCircleArc(circle_arc[0], circle_arc[2], circle_arc[1], circle_arc[3])
# Add all curve loops
for segment in curve_loops_rest_of_sphere:
    for curve_loop in segment:
        hole_on_sphere.addCurveLoop(curve_loop[0],curve_loop[1])
# Add all surfaces
for segment in surfaces_rest_of_sphere:
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
# PHYSICAL GROUPS FOR CURVES:
circle_arc = gmsh.model.addPhysicalGroup(1,[circle_arc[3] for circle_arc in list_of_lists[1]])
gmsh.model.setPhysicalName(1,circle_arc,"Boundary curve for the hole")
# PHYSICAL GROUPS FOR SURFACES:
# The hole itself
surface_hole = gmsh.model.addPhysicalGroup(2,[hole[1] for hole in list_of_lists[4]])
gmsh.model.setPhysicalName(2,surface_hole,"Surface hole")
# Adjacent hole
surface_adjacent_hole = gmsh.model.addPhysicalGroup(2,[adjacent_hole[1] for adjacent_hole in list_of_lists[5]])
gmsh.model.setPhysicalName(2, surface_adjacent_hole, "Surface adjacent hole")
# The rest of the sphere
surface_rest_of_sphere = gmsh.model.addPhysicalGroup(2,[surfaces_rest_of_sphere[i][0][1] for i in range(len(surfaces_rest_of_sphere))])
gmsh.model.setPhysicalName(2, surface_rest_of_sphere, "Surface the rest of the sphere")
# -------------------------------------------------------------------
# DEFINE COLOURS FOR THE MESH
# -------------------------------------------------------------------
# For nice colours, please see (https://colorbrewer2.org/)...
#-----------------------------------------------------------------------------
# THE REST OF THE SPEHRE
gmsh.model.setColor([(2,surfaces_rest_of_sphere[i][0][1]) for i in range(len(surfaces_rest_of_sphere))], 2, 56, 88)  # Darkest blue
# THE HOLE AND THE ADJACENT REGIONS
gmsh.model.setColor([(2,hole[1]) for hole in list_of_lists[4]], 255, 247, 251)  # Light blue
gmsh.model.setColor([(1,circle_arc[3]) for circle_arc in list_of_lists[1]], 103, 0, 31)  # Darkest magenta
gmsh.model.setColor([(2,adjacent[1]) for adjacent in list_of_lists[5]], 116, 169, 207)  # In between blue
# -------------------------------------------------------------------
# GENERATING THE MESH
# -------------------------------------------------------------------
# Generate the two dimensional mesh (as we work with surfaces) 
gmsh.model.mesh.generate(2)
# Write the mesh to our file
gmsh.write("../Meshes/sphere_with_1_hole.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
