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
#theta_vec = np.linspace(0,2*math.pi,100)
theta_vec = np.linspace(math.pi/6,math.pi/3,2)
#theta_vec = np.linspace(0,math.pi/10,4)
# Define the phi angle
#phi_vec = np.linspace(0,math.pi,100)
phi_vec = np.linspace(math.pi/10,math.pi,2)
# Define a corresponding mesh
phi_list = phi_vec.tolist()
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
# -------------------------------------------------------------------
# Looping over angles and defining all points, lines and curve loops 
# -------------------------------------------------------------------
# Make a list of previous points
prev_points = []
# Make a list of previous lines
prev_lines = []
# Loop over the theta vector (x,y-axes)
for theta_index in range(len(theta_vec)-1):
    print("\tTheta iter %d out of %d"%(theta_index+1,len(theta_vec)-1))
    # Extract a pair of theta angles
    theta_pair = [theta_vec[theta_index], theta_vec[theta_index +1]]
    print("\tTheta angles:\t%s"%(str(theta_pair)))
    # Previous index
    prev_index = 0
    # Loop over the phi angles (z-axis)
    for phi_index in range(len(phi_vec)-1):    
        print("\t\tPhi iter %d out of %d"%(phi_index+1,len(phi_vec)-1))
        # Extract a pair of phi angles
        phi_pair = [phi_vec[phi_index], phi_vec[phi_index + 1 ]]
        # Now, we will extract four points which will
        # be used to define four lines and then a curve loop
        # using these lines
        #---------------------------------------------------------
        # POINT 1
        if theta_index == 0:
            # Define the 3D point by its cartesian coordinates
            x1 = r * np.cos(theta_pair[0]) * np.sin(phi_pair[0])
            y1 = r * np.sin(theta_pair[0]) * np.sin(phi_pair[0])
            z1 = r * np.cos(phi_pair[0])
            # Add the point to our list
            point_list.append((x1,y1,z1,point_counter))
            # Add the point to our previous list
            prev_points.append((x1,y1,z1,point_counter))            
            # Set the point label
            point_label = point_counter
            # Increase the tag number of the point
            point_counter += 1
        else:
            x1 = prev_points[prev_index][0]
            y1 = prev_points[prev_index][1]
            z1 = prev_points[prev_index][2]
            point_label = prev_points[prev_index][3]
            prev_index += 1
        # Add the point to the mesh
        sphere_with_holes.addPoint(x1, y1, z1, Lc1, point_label)
        #---------------------------------------------------------
        # POINT 2
        # Define the 3D point by its cartesian coordinates
        x2 = r * np.cos(theta_pair[1]) * np.sin(phi_pair[0])
        y2 = r * np.sin(theta_pair[1]) * np.sin(phi_pair[0])
        z2 = r * np.cos(phi_pair[0])
        # Add the point to the mesh
        sphere_with_holes.addPoint(x2, y2, z2, Lc1, point_counter)
        # Add the point to our list
        point_list.append((x2,y2,z2,point_counter))
        # Increase the tag number of the point
        point_counter += 1
        #---------------------------------------------------------
        # POINT 3
        # Define the 3D point by its cartesian coordinates
        x3 = r * np.cos(theta_pair[1]) * np.sin(phi_pair[1])
        y3 = r * np.sin(theta_pair[1]) * np.sin(phi_pair[1])
        z3 = r * np.cos(phi_pair[1])
        # Add the point to the mesh
        sphere_with_holes.addPoint(x3, y3, z3, Lc1, point_counter)
        # Add the point to our list
        point_list.append((x3,y3,z3,point_counter))
        # Increase the tag number of the point
        point_counter += 1
        #---------------------------------------------------------
        # POINT 4
        if theta_index == 0:
            # Define the 3D point by its cartesian coordinates
            x4 = r * np.cos(theta_pair[0]) * np.sin(phi_pair[1])
            y4 = r * np.sin(theta_pair[0]) * np.sin(phi_pair[1])
            z4 = r * np.cos(phi_pair[1])
            # Add the point to our list
            point_list.append((x4,y4,z4,point_counter))
            # Add the point to our previous list
            prev_points.append((x4,y4,z4,point_counter))            
            # Set the point label
            point_label = point_counter
            # Increase the tag number of the point
            point_counter += 1
        else:
            x4 = prev_points[prev_index][0]
            y4 = prev_points[prev_index][1]
            z4 = prev_points[prev_index][2]
            point_label = prev_points[prev_index][3]
            prev_index += 1            
        # Add the point to the mesh
        sphere_with_holes.addPoint(x4, y4, z4, Lc1, point_label)
        #---------------------------------------------------------
        # POINT 1 (AGAIN)
        if theta_index == 0:
            # Define the 3D point by its cartesian coordinates
            x1 = r * np.cos(theta_pair[0]) * np.sin(phi_pair[0])
            y1 = r * np.sin(theta_pair[0]) * np.sin(phi_pair[0])
            z1 = r * np.cos(phi_pair[0])
            # Add the point to our list
            point_list.append((x1,y1,z1,point_counter))
            # Add the point to our previous list
            prev_points.append((x1,y1,z1,point_counter))            
            # Set the point label
            point_label = point_counter-4
            # Increase the tag number of the point
            point_counter += 1
        else:
            x1 = prev_points[prev_index][0]
            y1 = prev_points[prev_index][1]
            z1 = prev_points[prev_index][2]
            point_label = prev_points[prev_index][3]
            prev_index += 1
            # Add the point to our list
            point_list.append((x1,y1,z1,point_label))
            # Add the point to our previous list
            prev_points.append((x1,y1,z1,point_label))                        
        # Add the point to the mesh
        #sphere_with_holes.addPoint(x1, y1, z1, Lc1, point_label)
        #---------------------------------------------------------
        #-----------------------------------------------------------
        # THE LINES
        #-----------------------------------------------------------
        # The line counter
        line_counter = point_counter
        # Loop over the points and draw lines between them
        for i in range(0,len(point_list)):
            # Make sure that we loop until the second last point
            if i < len(point_list)-1:
                # Add the point
                sphere_with_holes.addLine(point_list[i][3], point_list[i+1][3], line_counter)
                # Append the new line to our list of lines
                line_list.append(line_counter)
                # Update the counter for the lines
                line_counter += 1
        #-----------------------------------------------------------
        # THE CURVE LOOPS
        #-----------------------------------------------------------
        # Define a counter controlling the tag of the curve loop
        curve_counter = line_counter+1
        # Add the curve loop  to the mesh
        sphere_with_holes.addCurveLoop(line_list, curve_counter)
        # Append our newly defined curve to our curve list
        curve_list.append(curve_counter)
        #-----------------------------------------------------------
        # RESET THE LISTS FOR THE POINTS AND THE LINES
        #-----------------------------------------------------------
        point_list = []
        line_list = []
        # Lastly, update the point_counter
        point_counter = curve_counter + 1 
        # End of inner for-loop that has to do with phi
# -------------------------------------------------------------------
# Assembling the surfaces over the sphere
# -------------------------------------------------------------------
# Introduce a counter for the surface
surface_counter = point_counter
# List of surfaces
surface_list = [surface_counter]
# Loop over the various curves and add them to the mesh
for curve in curve_list:
    # Add the curve loop as a surface
    sphere_with_holes.addPlaneSurface([curve], surface_counter)
    # Increase the counter keeping track of the surfaces
    surface_counter += 1
    # Save the surface at hand
    surface_list.append(surface_counter)
# -------------------------------------------------------------------
# Generating the mesh
# -------------------------------------------------------------------
# Synchronise everything with Gmsh
sphere_with_holes.synchronize()
# Set the colours of the surfaces
gmsh.model.setColor([(2, surface) for surface in surface_list], 221, 28, 119)  # Magenta
# For nice colours, please see (https://colorbrewer2.org/)
#----------------------------------
# Finally, we can add some comments by creating a post-processing view
# containing some strings:
v = gmsh.view.add("comments")

# Add a text string in window coordinates, 10 pixels from the left and 10 pixels
# from the bottom:
gmsh.view.addListDataString(v, [10, -10], ["Created with Gmsh"])

# Add a text string in model coordinates centered at (X,Y,Z) = (0, 0.11, 0),
# with some style attributes:
gmsh.view.addListDataString(v, [0, 0.11, 0], ["Hole"],
                            ["Align", "Center", "Font", "Helvetica"])

# If a string starts with `file://', the rest is interpreted as an image
# file. For 3D annotations, the size in model coordinates can be specified after
# a `@' symbol in the form `widthxheight' (if one of `width' or `height' is
# zero, natural scaling is used; if both are zero, original image dimensions in
# pixels are used):
gmsh.view.addListDataString(v, [0, 0.09, 0], ["file://../t4_image.png@0.01x0"],
                            ["Align", "Center"])

# The 3D orientation of the image can be specified by proving the direction
# of the bottom and left edge of the image in model space:
gmsh.view.addListDataString(v, [-0.01, 0.09, 0],
                            ["file://../t4_image.png@0.01x0,0,0,1,0,1,0"])

# The image can also be drawn in "billboard" mode, i.e. always parallel to
# the camera, by using the `#' symbol:
gmsh.view.addListDataString(v, [0, 0.12, 0],
                            ["file://../t4_image.png@0.01x0#"],
                            ["Align", "Center"])

# The size of 2D annotations is given directly in pixels:
gmsh.view.addListDataString(v, [150, -7], ["file://../t4_image.png@20x0"])

# These annotations are handled by a list-based post-processing view. For
# large post-processing datasets, that contain actual field values defined on
# a mesh, you should use model-based post-processing views instead, which
# allow to efficiently store continuous or discontinuous scalar, vector and
# tensor fields, or arbitrary polynomial order.

# Views and geometrical entities can be made to respond to double-click
# events, here to print some messages to the console:
gmsh.option.setString("View[0].DoubleClickedCommand",
                      "Printf('View[0] has been double-clicked!');")
gmsh.option.setString(
    "Geometry.DoubleClickedLineCommand",
    "Printf('Curve %g has been double-clicked!', "
    "Geometry.DoubleClickedEntityTag);")
# Generate the two dimensional mesh (as we work with surfaces) 
gmsh.model.mesh.generate(2)
# Write the mesh to our file
gmsh.write("../Meshes/sphere.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
