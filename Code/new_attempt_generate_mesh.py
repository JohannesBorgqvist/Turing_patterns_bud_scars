# =================================================================================
# =================================================================================
# Script:"generate_spherical_mesh"
# Date: 2021-06-28
# Implemented by: Johannes Borgqvist
# Description:
# The program generates a spherical mesh with
# a given radius, and it generates a given
# number of holes in the mesh determined by
# the input.
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import gmsh  # For generating the meshes
import math  # For using mathematical functions
import sys  # Needed by gmsh to launch the GUI, i.e. open the window with the plot
import numpy as np  # Needed to create arrays
# =================================================================================
# =================================================================================
# Code for generating the mesh
# =================================================================================
# =================================================================================
# Initialise the Gmsh framework
gmsh.initialize()
# Add a model entity
gmsh.model.add("sphere_with_holes")
# We start by defining some points and some lines. To make the code shorter we
# can redefine a namespace:
sphere_with_holes = gmsh.model.geo
# -------------------------------------------------------------------
# Allocating memory for all lists and vectors
# -------------------------------------------------------------------
# Define a centimetre
cm = 1e-02
Lc1 = 0.01
r = 1 * cm
# Define the theta angle
theta_vec = np.linspace(0,2*math.pi,100)
#theta_vec = np.linspace(0,math.pi/10,3)
# Define the phi angle
phi_vec = np.linspace(0,math.pi,100)
# Define all points that we want to add
point_list = []
# Define a point counter
point_counter = 1
# Define a line counter
line_counter = 1
# Define a curve counter
curve_counter = 1
# Define a list for the lines
line_list = []
# Define a list for the curves
curve_list = []
#--------------------------------------------------------------------
# DEFINE THE TWO POLES
#--------------------------------------------------------------------
# Define the north pole which will be our start point every time
north_pole = (0,0,r,point_counter)
# Add the point to the mesh
sphere_with_holes.addPoint(north_pole[0], north_pole[1], north_pole[2], Lc1, north_pole[3])        
# Increment the point counter
point_counter += 1
# Define the south pole which will be our start point every time
south_pole = (0,0,-r,point_counter)
# Add the point to the mesh
sphere_with_holes.addPoint(south_pole[0], south_pole[1], south_pole[2], Lc1, south_pole[3])        
# Increment the point counter
point_counter += 1
# -------------------------------------------------------------------
# DEFINING ALL POINTS IN THE MESH
# -------------------------------------------------------------------
# Loop over the theta vector (x,y-axes)
for theta in theta_vec:
    # Allocate a list for each line segment
    temp_points = []
    # Loop over the phi angles (z-axis)
    for phi in phi_vec:
        if phi == 0:
            # Add the point to our list
            temp_points.append(north_pole)
        elif phi == math.pi:
            # Add the point to our list
            temp_points.append(south_pole)
        else:
            # Define the 3D point by its cartesian coordinates
            x = r * np.cos(theta) * np.sin(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(phi)
            # Add the point to the mesh
            sphere_with_holes.addPoint(x, y, z, Lc1, point_counter)        
            # Add the point to our list
            temp_points.append((x,y,z,point_counter))
            # Increment the point counter
            point_counter += 1
    # Save the points for the first line segment
    point_list.append(temp_points)
# -------------------------------------------------------------------
# DEFINING ALL LINES IN THE MESH
# -------------------------------------------------------------------
# Introduce a line counter
line_counter = point_counter
#print("The line counter:")
#print(line_counter)
#print("The point_list contains %d lists"%(int(len(point_list))))
#print(line_counter)
# Divide the lines into two classes,
# where one has one direction of its lines
# and the other has the opposite direction
line_classifier = 0
# Add all lines
for line_segments in point_list:
    # Allocate an empty list for the lines in this segment
    temp_lines = []
    # Loop over all points in our line segments
    # and combine them to lines
    for i in range(len(line_segments)):
        # Make sure that we loop until the second last point
        if i < len(line_segments)-1:
            # For the even ones we add lines in one direction
            if line_classifier%2 == 0:
                # Add the point
                sphere_with_holes.addLine(line_segments[i][3], line_segments[i+1][3], line_counter)
                # Append the new line to our list of lines
                temp_lines.append((line_segments[i][3],line_segments[i+1][3],line_counter))
            else:
                # Add the point
                sphere_with_holes.addLine(line_segments[i+1][3],line_segments[i][3], line_counter)
                # Append the new line to our list of lines
                temp_lines.append((line_segments[i+1][3],line_segments[i][3],line_counter))                
            # Update the counter for the lines
            line_counter += 1
    # Save all lines in our line list
    if line_classifier%2 == 0: # Even lines
        line_list.append(temp_lines) # Standard direction
    else: # Odd lines
        line_list.append(temp_lines[::-1]) # Reverse direction
    # Lastly, increment the line classifier
    line_classifier += 1
# -------------------------------------------------------------------
# Add lines between each arch and define curves
# -------------------------------------------------------------------
# Allocate memory for the inbetween lines
in_between_lines = []
# Define a list with all curves
curve_list = []
# Define a counter for all curves
curve_counter = 1
# Loop through all points and add in between lines
for arc_index in range(len(line_list)-1):
    # Create a between index
    between_index = -1
    # Allocate a temporary list
    in_between_line = []
    #Loop over all lines in the arc
    for line_index in range(len(line_list[arc_index])):
        # Start with the even lines
        if arc_index %2 == 0:
            # Take the standard surface elements in the middle    
            if line_index != 0 and line_index != len(line_list[arc_index])-1: # Normal case
                # Extract the first point on curve 1
                point_1_on_curve_1 = line_list[arc_index][line_index][0]
                # Extract the second point on curve 1
                point_2_on_curve_1 = line_list[arc_index][line_index][1]            
                # Extract the first point on curve 2
                point_1_on_curve_2 = line_list[arc_index+1][len(line_list[arc_index])-(line_index+1)][0]
                # Extract the second point on curve 2
                point_2_on_curve_2 = line_list[arc_index+1][len(line_list[arc_index])-(line_index+1)][1]            
                # Add the first in-between line
                sphere_with_holes.addLine(point_2_on_curve_1,point_1_on_curve_2, line_counter)
                # Save the first in between line as a 3-tuple
                in_between_line.append((point_2_on_curve_1,point_1_on_curve_2, line_counter))
                # Increment the line_counter
                line_counter += 1
                # Add the second in-between line
                sphere_with_holes.addLine(point_2_on_curve_2,point_1_on_curve_1, line_counter)
                # Save the second in between line as a 3-tuple
                in_between_line.append((point_2_on_curve_2,point_1_on_curve_1, line_counter))
                # Increment the line_counter
                line_counter += 1
                # Define a mini-curve
                mini_curve = [line_list[arc_index][line_index][2], line_counter-2,line_list[arc_index+1][len(line_list[arc_index])-(line_index+1)][2], line_counter-1]
                # Increment the curve_counter
                curve_counter=line_counter
                # Save the mini-curve as a curve_loop
                sphere_with_holes.addCurveLoop(mini_curve, curve_counter)
                # Append the curve at hand to the curve_list
                curve_list.append(curve_counter)
                # Increment the line_counter
                line_counter += 1
            elif line_index == 0:# Special case, first surface element
                # Extract the second point on curve 1
                point_2_on_curve_1 = line_list[arc_index][line_index][1]            
                # Extract the first point on curve 2
                point_1_on_curve_2 = line_list[arc_index+1][len(line_list[arc_index])-(line_index+1)][0]
                # Add the second in-between line
                sphere_with_holes.addLine(point_2_on_curve_1,point_1_on_curve_2, line_counter)
                # Save the second in between line as a 3-tuple
                in_between_line.append((point_2_on_curve_1,point_1_on_curve_2, line_counter))
                # Increment the line_counter
                line_counter += 1                
                # Define a mini-curve
                mini_curve = [line_list[arc_index][line_index][2], line_counter-1,line_list[arc_index+1][len(line_list[arc_index])-(line_index+1)][2]]
                # Increment the curve_counter
                curve_counter=line_counter
                # Save the mini-curve as a curve_loop
                sphere_with_holes.addCurveLoop(mini_curve, curve_counter)
                # Append the curve at hand to the curve_list
                curve_list.append(curve_counter)
                # Increment the line_counter
                line_counter += 1                
            elif line_index == len(line_list[arc_index])-1:
                # Extract the second point on curve 1
                point_1_on_curve_1 = line_list[arc_index][line_index][0]            
                # Extract the first point on curve 2
                point_2_on_curve_2 = line_list[arc_index+1][len(line_list[arc_index])-(line_index+1)][1]
                # Add the second in-between line
                sphere_with_holes.addLine(point_2_on_curve_2,point_1_on_curve_1, line_counter)
                # Save the second in between line as a 3-tuple
                in_between_line.append((point_2_on_curve_2,point_1_on_curve_1, line_counter))
                # Increment the line_counter
                line_counter += 1                
                # Define a mini-curve
                mini_curve = [line_list[arc_index][line_index][2],line_list[arc_index+1][len(line_list[arc_index])-(line_index+1)][2],line_counter-1]
                # Increment the curve_counter
                curve_counter=line_counter
                # Save the mini-curve as a curve_loop
                sphere_with_holes.addCurveLoop(mini_curve, curve_counter)
                # Append the curve at hand to the curve_list
                curve_list.append(curve_counter)
                # Increment the line_counter
                line_counter += 1                                
            # Increment the between_index
            between_index += 2
        else: # Odd lines
            # Take the standard surface elements in the middle    
            if line_index != 0 and line_index != len(line_list[arc_index])-1: # Normal case
                # Extract the first point on curve 2
                point_1_on_curve_2 = line_list[arc_index+1][line_index][0]
                # Extract the second point on curve 2
                point_2_on_curve_2 = line_list[arc_index+1][line_index][1]
                # Extract the second point on curve 1
                point_2_on_curve_1 = line_list[arc_index][len(line_list[arc_index])-(line_index+1)][1]   
                # Extract the first point on curve 1
                point_1_on_curve_1 = line_list[arc_index][len(line_list[arc_index])-(line_index+1)][0]
                # Add the second in-between line
                sphere_with_holes.addLine(point_2_on_curve_1,point_1_on_curve_2, line_counter)
                # Save the first in between line as a 3-tuple
                in_between_line.append((point_2_on_curve_1,point_1_on_curve_2, line_counter))
                # Increment the line_counter
                line_counter += 1
                # Add the second in-between line
                sphere_with_holes.addLine(point_2_on_curve_2,point_1_on_curve_1, line_counter)
                # Save the second in between line as a 3-tuple
                in_between_line.append((point_2_on_curve_2,point_1_on_curve_1, line_counter))
                # Increment the line_counter
                line_counter += 1
                # Define a mini-curve
                mini_curve = [line_list[arc_index][len(line_list[arc_index])-(line_index+1)][2], line_counter-1,line_list[arc_index+1][line_index][2], line_counter-2]
                # Increment the curve_counter
                curve_counter=line_counter
                # Save the mini-curve as a curve_loop
                sphere_with_holes.addCurveLoop(mini_curve, curve_counter)
                # Append the curve at hand to the curve_list
                curve_list.append(curve_counter)
                # Increment the line_counter
                line_counter += 1
            elif line_index == 0:# Special case, first surface element
                # Extract the second point on curve 1
                point_1_on_curve_1 = line_list[arc_index][len(line_list[arc_index])-(line_index+1)][0]            
                # Extract the first point on curve 2
                point_2_on_curve_2 = line_list[arc_index+1][line_index][1]
                # Add the second in-between line
                sphere_with_holes.addLine(point_2_on_curve_2,point_1_on_curve_1, line_counter)
                # Save the second in between line as a 3-tuple
                in_between_line.append((point_2_on_curve_2,point_1_on_curve_1, line_counter))
                # Increment the line_counter
                line_counter += 1                
                # Define a mini-curve
                mini_curve = [line_list[arc_index+1][line_index][2], line_counter-1,line_list[arc_index][len(line_list[arc_index])-(line_index+1)][2]]
                # Increment the curve_counter
                curve_counter=line_counter
                # Save the mini-curve as a curve_loop
                sphere_with_holes.addCurveLoop(mini_curve, curve_counter)
                # Append the curve at hand to the curve_list
                curve_list.append(curve_counter)
                # Increment the line_counter
                line_counter += 1
            elif line_index == len(line_list[arc_index])-1:
                # Extract the second point on curve 1
                point_2_on_curve_1 = line_list[arc_index][len(line_list[arc_index])-(line_index+1)][1]            
                # Extract the first point on curve 2
                point_1_on_curve_2 = line_list[arc_index+1][line_index][0]
                # Add the second in-between line
                sphere_with_holes.addLine(point_2_on_curve_1,point_1_on_curve_2, line_counter)
                # Save the second in between line as a 3-tuple
                in_between_line.append((point_2_on_curve_1,point_1_on_curve_2, line_counter))
                # Increment the line_counter
                line_counter += 1                
                # Define a mini-curve
                mini_curve = [line_list[arc_index+1][line_index][2],line_list[arc_index][len(line_list[arc_index])-(line_index+1)][2],line_counter-1]
                # Increment the curve_counter
                curve_counter=line_counter
                # Save the mini-curve as a curve_loop
                sphere_with_holes.addCurveLoop(mini_curve, curve_counter)
                # Append the curve at hand to the curve_list
                curve_list.append(curve_counter)
                # Increment the line_counter
                line_counter += 1                                
            # Increment the between_index
            between_index += 2                
        #--------------------------------------------------------------------
        # End of special cases
        #--------------------------------------------------------------------
        # Increment the curve_counter
        #curve_counter=line_counter+1
        # Save the mini-curve as a curve_loop
        #sphere_with_holes.addCurveLoop(mini_curve, curve_counter)
        # Append the curve at hand to the curve_list
        #curve_list.append(curve_counter)    
        # Increment the between_index
        #between_index += 1
        # Increment the line_counter
        #line_counter += 1
    # Save all in-between lines
    in_between_lines.append(in_between_line)        
# -------------------------------------------------------------------
# DEFINING ALL CURVE LOOPS IN THE MESH
# -------------------------------------------------------------------
# Define a list with all curves
#curve_list = []
# Define a counter for all curves
#curve_counter = line_counter+1
# Loop through all arcs and define a curve for each little surface element
#for arc_index in range(len(line_list)-1):
    
# Define a counter controlling the tag of the curve loop
#curve_counter = line_counter+1
#mini_curve = [line_list[0][0][2], line_counter,line_list[1][2][2]]
#print("Print the little curve")
#print(mini_curve)
# Add the curve loop  to the mesh
#sphere_with_holes.addCurveLoop(mini_curve, curve_counter)
# Append our newly defined curve to our curve list
#curve_list.append(curve_counter)

# -------------------------------------------------------------------
# DEFINING ALL SURFACE IN THE MESH
# -------------------------------------------------------------------
#print("All our defined curves are:")
#print(curve_list)
# Define a surface counter then?
surface_counter = curve_counter
# Define a list of surfaces
surface_list = []


# Add all surfaces now
for curve in curve_list:
    # Increment the surface counter
    surface_counter += 1
    # Add the surface
    sphere_with_holes.addPlaneSurface([curve], surface_counter)
    # Save the surface
    surface_list.append(surface_counter)

#sphere_with_holes.addPlaneSurface(curve_list, 666)
    
#gmsh.model.occ.fuse ([(2, surface_list[index]) for index in range(len(surface_list))], [(2, surface_list[len(surface_list)-1])], 666)
# -------------------------------------------------------------------
# GENERATING THE MESH
# -------------------------------------------------------------------
# For nice colours, please see (https://colorbrewer2.org/)
# Synchronise everything with Gmsh
sphere_with_holes.synchronize()
# Set the colours of the surfaces
gmsh.model.setColor([(2, surface) for surface in surface_list], 221, 28, 119)  # Magenta
#----------------------------------
# Finally, we can add some comments by creating a post-processing view
# containing some strings:
#v = gmsh.view.add("comments")

# Add a text string in window coordinates, 10 pixels from the left and 10 pixels
# from the bottom:
#gmsh.view.addListDataString(v, [10, -10], ["Created with Gmsh"])
#
# Add a text string in model coordinates centered at (X,Y,Z) = (0, 0.11, 0),
# with some style attributes:
#gmsh.view.addListDataString(v, [0, 0.11, 0], ["Hole"],
#                            ["Align", "Center", "Font", "Helvetica"])

# If a string starts with `file://', the rest is interpreted as an image
# file. For 3D annotations, the size in model coordinates can be specified after
# a `@' symbol in the form `widthxheight' (if one of `width' or `height' is
# zero, natural scaling is used; if both are zero, original image dimensions in
# pixels are used):
#gmsh.view.addListDataString(v, [0, 0.09, 0], ["file://../t4_image.png@0.01x0"],
#                            ["Align", "Center"])

# The 3D orientation of the image can be specified by proving the direction
# of the bottom and left edge of the image in model space:
#gmsh.view.addListDataString(v, [-0.01, 0.09, 0],
#                            ["file://../t4_image.png@0.01x0,0,0,1,0,1,0"])

# The image can also be drawn in "billboard" mode, i.e. always parallel to
# the camera, by using the `#' symbol:
#gmsh.view.addListDataString(v, [0, 0.12, 0],
#                            ["file://../t4_image.png@0.01x0#"],
#                            ["Align", "Center"])

# The size of 2D annotations is given directly in pixels:
#gmsh.view.addListDataString(v, [150, -7], ["file://../t4_image.png@20x0"])

# These annotations are handled by a list-based post-processing view. For
# large post-processing datasets, that contain actual field values defined on
# a mesh, you should use model-based post-processing views instead, which
# allow to efficiently store continuous or discontinuous scalar, vector and
# tensor fields, or arbitrary polynomial order.

# Views and geometrical entities can be made to respond to double-click
# events, here to print some messages to the console:
#gmsh.option.setString("View[0].DoubleClickedCommand",
#                      "Printf('View[0] has been double-clicked!');")
#gmsh.option.setString(
#    "Geometry.DoubleClickedLineCommand",
#    "Printf('Curve %g has been double-clicked!', "
#    "Geometry.DoubleClickedEntityTag);")
# Generate the two dimensional mesh (as we work with surfaces) 
gmsh.model.mesh.generate(2)
# Write the mesh to our file
gmsh.write("../Meshes/sphere.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
