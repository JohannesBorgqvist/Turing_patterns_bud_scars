# =================================================================================
# =================================================================================
# Script:"generate_spherical_mesh"
# Date: 2021-07-01
# Implemented by: Johannes Borgqvist
# Description:
# The program generates a bud scar on the
# sphere with radius r. 
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import gmsh  # For generating the meshes
import math  # For using mathematical functions
import sys  # Needed by gmsh to launch the GUI, i.e. open the window with the plot
import numpy as np  # Needed to create linspace arrays and other nice functions
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
hole_on_sphere = gmsh.model.geo
# -------------------------------------------------------------------
# Functions
# -------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# Function 1: "add_circular_holes_on_sphere"
# The function adds a circular hole onto the sphere by defining a
# subdomain of the based on the intersection between a big sphere
# and a smaller sphere. This intersection is a circle and it defines
# the hole in this case. 
# The function takes the following four inputs:
# 1. "list_of_centres": a list containing two tuples of the type c=(x,y,z)
# declaring the centre of the big sphere and the smaller intersecting sphere,
# the small intersecting sphere, 
# 2. "r": a value of the radius of the intersecting small sphere,
# 3. "R": a value of the radius of the big sphere with centre in the
# origin (0,0,0),
# 4. "label_counter": a counter defining the numbers of the points, curves
# and surfaces. 
# The function returns a three outputs:
# 1. "list_of_lists": Contains the list of all points, curves, curve loops,
# surfaces for the holes and the surfaces of the adjacent regions. 
# 2."theta_tuples": a tuple of the type (theta_min,theta_max) defining the
# boundaries of the square containing the circular hole,
# 3."phi_tuples": a tuple of the type (phi_min,phi_max) defining the
# boundaries of the square containing the circular hole.
# Using all the lists in "list_of_lists" it is possible to create the part
# of the mesh containing all holes. Using the boundaries it is possible to
# create the rest of the mesh by completing the remainder of the surfaces
# of the sphere.
#-----------------------------------------------------------------------------------
def add_circular_holes_on_sphere(list_of_centres,r,R,label_counter,step_size):
    # Allocate memory for the outputs (what we will return).
    list_of_lists = [] # Points, curves and surfaces
    theta_list = [] # Boundaries of holes in the theta direction
    phi_list = [] # Boundaries of holes in the theta direction
    # Define the thickness of one point
    Lc1 = 0.01
    # Note that the latter two angles will be updated before they are
    # returned.
    # Allocate memory for the lists that the "list_of_lists" will be
    # comprised of.
    list_of_points = [] # Points
    list_of_lines = [] # Lines
    list_of_circle_arcs = [] # Circle arcs
    list_of_curve_loops = [] # Curve loops
    list_of_surfaces_hole = [] # Surfaces of hole
    list_of_surfaces_adjacent = [] # Region adjacent to the hole    
    #----------------------------------------------------------------
    # Extract centres
    #----------------------------------------------------------------    
    c1 = list_of_centres[0] # Centre of big sphere
    c2 = list_of_centres[1] # Centre of small sphere
    #----------------------------------------------------------------
    # Define centre and bases (tangent ti and bitangent bi) for the
    # tangent plane of the intersecting circle
    #----------------------------------------------------------------
    # Distance between centres
    d = np.sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2+(c1[2]-c2[2])**2)
    # Define a weighting coefficient
    h = 1/2 + ( (R**2-r**2) / (2*d**2))
    # Define a linear combination of our centres which defines where
    # our beloved circle is located
    ci = tuple([c1[i] + h*(c2[i]-c1[i]) for i in range(len(c1))])
    # Define a trigonometric entity which is the radius of the intersecting circle
    ri = np.sqrt(R**2-((h**2)*(d**2)))
    # Define a normal to this plane
    ni = tuple([round(((c2[i]-c1[i])/(d)),4) for i in range(len(list(c1)))])
    # Define the tangent ti as a cross product of one of the unit basis
    # vectors
    if np.cross([1, 0, 0],list(ni))[0]!=0 or np.cross([1, 0, 0],list(ni))[1]!=0 or np.cross([1, 0, 0],list(ni))[2]!=0:
        ti = np.cross([1, 0, 0],list(ni)) # Take cross product
    elif np.cross([0, 1, 0],list(ni))[0]!=0 or np.cross([0, 1, 0],list(ni))[1]!=0 or np.cross([0, 1, 0],list(ni))[2]!=0:
        ti = np.cross([0, 1, 0],list(ni)) # Take cross product
    elif np.cross([0, 0, 1],list(ni))[0]!=0 or np.cross([0, 0, 1],list(ni))[1]!=0 or np.cross([0, 0, 1],list(ni))[2]!=0:
        ti = np.cross([0, 0, 1],list(ni)) # Take cross product        
    # Normalise the same vector
    ti = ti/np.linalg.norm(ti) 
    # Define the bitangent bi
    bi = np.cross(ti,ni)
    #----------------------------------------------------------------
    # Define the hole or the intersecting circle
    #----------------------------------------------------------------
    # Define the mid point of the intersecting circle
    list_of_points.append((ci[0], ci[1], ci[2], Lc1, label_counter)) # Add the point
    label_counter += 1    # Increment label counter
    # Define first point on intersecting circle
    theta = 0 # Angle of the point
    p1 = tuple([ci[j] + ri*(ti[j]*np.cos(theta)+bi[j]*np.sin(theta)) for j in range(len(ci))])
    list_of_points.append((p1[0], p1[1], p1[2], Lc1, label_counter)) # Add the point
    label_counter += 1    # Increment label counter    
    # Define second point on circle
    theta = math.pi # Angle of the point
    p2 = tuple([ci[j] + ri*(ti[j]*np.cos(theta)+bi[j]*np.sin(theta)) for j in range(len(ci))])
    list_of_points.append((p2[0], p2[1], p2[2], Lc1, label_counter)) # Add the point
    label_counter += 1    # Increment label counter
    # Define third point on circle
    theta = math.pi/2 # Angle of the point
    p3 = tuple([ci[j] + ri*(ti[j]*math.cos(theta)+bi[j]*math.sin(theta)) for j in range(len(ci))])
    list_of_points.append((p3[0], p3[1], p3[2], Lc1, label_counter)) # Add the point
    label_counter += 1    # Increment label counter        
    # Define fourth point on circle
    theta = 3*math.pi/2 # Angle of the point
    p4 = tuple([ci[j] + ri*(ti[j]*math.cos(theta)+bi[j]*math.sin(theta)) for j in range(len(ci))])
    list_of_points.append((p4[0], p4[1], p4[2], Lc1, label_counter)) # Add the point
    # CIRCLE ARCS
    # Define four circle arcs between these five points
    list_of_circle_arcs.append((label_counter-3,label_counter-4,label_counter-1,label_counter+1)) # Circle arc 1
    list_of_circle_arcs.append((label_counter-1,label_counter-4,label_counter-2,label_counter+2)) # Circle arc 2
    list_of_circle_arcs.append((label_counter-2,label_counter-4,label_counter,label_counter+3)) # Circle arc 3
    list_of_circle_arcs.append((label_counter,label_counter-4,label_counter-3,label_counter+4)) # Circle arc 4    
    # Curve loop for the circle arcs
    list_of_curve_loops.append([[label_counter+1,label_counter+2,label_counter+3,label_counter+4],label_counter+5])
    # Surface based on circle arc
    list_of_surfaces_hole.append([[label_counter+5],label_counter+6])
    # Increment the label counter
    label_counter += 7
    #----------------------------------------------------------------
    # Define the square around the circle 
    #----------------------------------------------------------------
    # Corner 1 of square
    temp_list =[ci[j] + 1.75*ri*(ti[j]*math.cos(math.pi/4)+bi[j]*math.sin(math.pi/4)) for j in range(len(ci)-1)]
    temp_list.append(np.sqrt(R**2-temp_list[0]**2-temp_list[1]**2))
    p5 = tuple(temp_list)
    # Corner 2 of square    
    temp_list =[ci[j] + 1.75*ri*(ti[j]*math.cos(3*math.pi/4)+bi[j]*math.sin(3*math.pi/4)) for j in range(len(ci)-1)]
    temp_list.append(np.sqrt(R**2-temp_list[0]**2-temp_list[1]**2))    
    p6 = tuple(temp_list)
    # Corner 3 of square 
    temp_list =[ci[j] + 1.75*ri*(ti[j]*math.cos(5*math.pi/4)+bi[j]*math.sin(5*math.pi/4)) for j in range(len(ci)-1)]
    temp_list.append(np.sqrt(R**2-temp_list[0]**2-temp_list[1]**2))
    p7 = tuple(temp_list)
    # Corner 4 of square 
    temp_list =[ci[j] + 1.75*ri*(ti[j]*math.cos(7*math.pi/4)+bi[j]*math.sin(7*math.pi/4)) for j in range(len(ci)-1)]
    temp_list.append(np.sqrt(R**2-temp_list[0]**2-temp_list[1]**2))
    p8 = tuple(temp_list)
    # Create a list of tuples
    tuple_5 = (np.arctan2(p5[1],p5[0]),np.arctan2(np.sqrt(p5[0]**2+p5[1]**2),p5[2]),0)
    tuple_6 = (np.arctan2(p6[1],p6[0]),np.arctan2(np.sqrt(p6[0]**2+p6[1]**2),p6[2]),1)
    tuple_7 = (np.arctan2(p7[1],p7[0]),np.arctan2(np.sqrt(p7[0]**2+p7[1]**2),p7[2]),2)
    tuple_8 = (np.arctan2(p8[1],p8[0]),np.arctan2(np.sqrt(p8[0]**2+p8[1]**2),p8[2]),3)
    list_of_tuples = [tuple_5, tuple_6, tuple_7, tuple_8]
    #tup_theta = sort_tuple(list_of_tuples,0)
    list_phi = sorted(list_of_tuples, key=lambda tup: tup[1])   
    list_theta = sorted(list_of_tuples, key=lambda tup: tup[0])   
    # Extract the minimum and maximum angles
    phi_tuple = (list_phi[0][1], list_phi[len(list_phi)-1][1])
    theta_tuple = (list_theta[0][0], list_theta[len(list_theta)-1][0])
    # Add these four points instead:
    # Point 1
    x = R * np.cos(theta_tuple[0]) * np.sin(phi_tuple[0])
    y = R * np.sin(theta_tuple[0]) * np.sin(phi_tuple[0])
    z = R * np.cos(phi_tuple[0])
    list_of_points.append((x, y, z, Lc1, label_counter)) # Add the point
    label_counter += 1    # Increment label counter
    # Point 2
    x = R * np.cos(theta_tuple[0]) * np.sin(phi_tuple[1])
    y = R * np.sin(theta_tuple[0]) * np.sin(phi_tuple[1])
    z = R * np.cos(phi_tuple[1])
    list_of_points.append((x, y, z, Lc1, label_counter)) # Add the point
    label_counter += 1    # Increment label counter
    # Point 3
    x = R * np.cos(theta_tuple[1]) * np.sin(phi_tuple[1])
    y = R * np.sin(theta_tuple[1]) * np.sin(phi_tuple[1])
    z = R * np.cos(phi_tuple[1])
    list_of_points.append((x, y, z, Lc1, label_counter)) # Add the point
    label_counter += 1    # Increment label counter
    # Point 4
    x = R * np.cos(theta_tuple[1]) * np.sin(phi_tuple[0])
    y = R * np.sin(theta_tuple[1]) * np.sin(phi_tuple[0])
    z = R * np.cos(phi_tuple[0])
    list_of_points.append((x, y, z, Lc1, label_counter)) # Add the point
    label_counter += 1    # Increment label counter
    # Save a magic label for later use
    magic_label = label_counter
    # Print the corner points
    corner_points = [label_counter-4,label_counter-3,label_counter-2,label_counter-1]
    # ASSEMBLE LINES FOR THE SQUARE
    # Calculate increment in both directions
    theta_increment = step_size
    phi_increment = step_size
    # Side 1:
    # Loop over sides
    for side in range(1,5):
        # Allocate a list for all point on
        # one side of the "square" in the
        # (theta,phi)-plane
        side_points = []        
        # Allocate some temporary angles
        phi = 0
        theta = 0
        # Some special cases depending which side we are on
        if side == 1:
            start_point = magic_label - 4
            end_point = magic_label - 3
            increment = phi_increment
            theta = theta_tuple[0]
            start_angle = phi_tuple[0]
        elif side == 2:
            start_point = magic_label - 3
            end_point = magic_label - 2
            increment = theta_increment
            phi = phi_tuple[1]
            start_angle = theta_tuple[0]            
        elif side == 3:
            start_point = magic_label - 2
            end_point = magic_label - 1
            increment = phi_increment
            theta = theta_tuple[1]
            start_angle = phi_tuple[0]            
        elif side == 4:
            start_point = magic_label - 1
            end_point = magic_label - 4
            increment = -theta_increment
            phi = phi_tuple[0]
            start_angle = theta_tuple[1]
        # Add the start point to the side
        side_points.append(start_point)
        # Calculate number of increments
        if side%2 == 0: # Even side
            number_of_points = int(round(((theta_tuple[1]-theta_tuple[0])/(step_size))))
        else: # Odd side
            number_of_points = int(round(((phi_tuple[1]-phi_tuple[0])/(step_size))))            
        # Start to add points along the side
        for increment_index in range(1,number_of_points):            
            # Now, we define the phi or theta angle
            # depending if we work with even or odd sides
            if side %2 == 0: # Even side
                if increment_index == 1:
                    theta = start_angle+increment
                else:
                    theta += increment
                    # Add our beloved point
                    x = R * np.cos(theta) * np.sin(phi)
                    y = R * np.sin(theta) * np.sin(phi)
                    z = R * np.cos(phi)
                    list_of_points.append((x, y, z, Lc1, label_counter)) # Add the point to the mesh
                    side_points.append(label_counter) # Add the points to what defines the side
                    label_counter += 1    # Increment label counter                    
            else: # Odd side
                if increment_index == 1:
                    phi = start_angle+increment
                else:
                    phi += increment
                # Add our beloved point
                x = R * np.cos(theta) * np.sin(phi)
                y = R * np.sin(theta) * np.sin(phi)
                z = R * np.cos(phi)
                list_of_points.append((x, y, z, Lc1, label_counter)) # Add the point to the mesh
                side_points.append(label_counter) # Add the points to what defines the side
                label_counter += 1    # Increment label counter                    
        # Add the end point to the side
        side_points.append(end_point)
        # Loop over all points on the side, draw lines between the points
        # in order to define the sides
        for side_point_index in range(len(side_points)-1):
            # Draw a line
            list_of_lines.append((side_points[side_point_index],side_points[side_point_index+1],label_counter))
            # Increment the label counter
            label_counter += 1
    # Add a curve loop for all our lines
    list_of_curve_loops.append([[lines[2] for lines in list_of_lines],label_counter+1])
    # Define a surface between our circle and the area around
    list_of_surfaces_adjacent.append([[curve[1] for curve in list_of_curve_loops],label_counter+1])        
    # Update the list, we return
    list_of_lists.append(list_of_points) # Points
    list_of_lists.append(list_of_lines) # Lines
    list_of_lists.append(list_of_circle_arcs) # Arcs
    list_of_lists.append(list_of_curve_loops) # Curve loops
    list_of_lists.append(list_of_surfaces_hole) # Surfaces holes
    list_of_lists.append(list_of_surfaces_adjacent) # Surfaces adjacent
    # Update the new label
    new_label = label_counter
    # Return the results
    return new_label, list_of_lists, theta_tuple, phi_tuple
#-----------------------------------------------------------------------------------
# Function 2: "add_parts_of_sphere"
# The function generates the mesh for parts of the sphere
# with centre in the origin (0,0,0) and radius R. 
# The function takes the following four inputs:
# 1. "R": a value of the radius of the sphere with centre in the
# origin (0,0,0),
# 2. "phi_tuple": A tuple of the type phi_tuple=(phi_min,phi_max)
# which determines the extent of the mesh in the phi-direction,
# 3. "theta_tuple": A tuple of the type theta_tuple=(theta_min,theta_max)
# which determines the extent of the mesh in the theta-direction,
# 4. "step_size": determines the density of the mesh,
# 5. "object_counter": counter which keeps track of all points,
# lines, curve loops and surfaces,
# 6. "north_pole": a tuple determining the north pole given by
# (0,0,r,label)
# 6. "south_pole": a tuple determining the north pole given by
# (0,0,-r,label) 
# The function returns a two outputs:
# 1. "new_counter" which is a variable counting the latest object in the mesh,
# 2. "list_of_lists" which is a list containing sublists for the points, the lines
# the curve_loops and the surfaces.
#-----------------------------------------------------------------------------------
def add_parts_of_sphere(R,phi_tuple,theta_tuple,step_size,object_counter,nort_pole,south_pole):
    #-----------------------------------------------------------------------
    # Allocate memory for the ouput
    #-----------------------------------------------------------------------    
    # Allocate memory for the outputs (what we will return).
    list_of_lists = [] # Points, curves and surfaces
    # Allocate memory for the lists that will be included in the    
    # larger list of lists
    list_of_points = [] # Points
    list_of_lines = [] # Lines
    list_of_curve_loops = [] # Curve loops 
    list_of_surfaces = [] # Surfaces of rest if the sphere
    #-----------------------------------------------------------------------
    # Set up the vectors for the angles 
    #-----------------------------------------------------------------------    
    # Define a step size for the phi vector
    phi_step = max(int(round(((phi_tuple[1]-phi_tuple[0])/(step_size)))),2)
    # Define a step size for the theta vector
    theta_step = max(int(round(((theta_tuple[1]-theta_tuple[0])/(step_size)))),2)
    # Define our two linspace vectors which we will loop over
    phi_vec = np.linspace(phi_tuple[0],phi_tuple[1],phi_step) # Phi vector
    theta_vec = np.linspace(theta_tuple[0],theta_tuple[1],theta_step) # Theta vector    # Define a measure of the thickness of each point
    Lc1 = 0.01
    #-----------------------------------------------------------------------
    # Define all points in the list
    #-----------------------------------------------------------------------    
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
                if phi != 0 or phi != math.pi:
                    # Define the 3D point by its cartesian coordinates
                    x = R * np.cos(theta) * np.sin(phi)
                    y = R * np.sin(theta) * np.sin(phi)
                    z = R * np.cos(phi)
                    # Add the point to our list
                    temp_points.append((x,y,z,Lc1,object_counter))
                    # Increment the point counter
                    object_counter += 1
        # Save the points for the line segment that
        # was just computed
        list_of_points.append(temp_points)
    # -------------------------------------------------------------------
    # DEFINING ALL LINES IN THE MESH
    # -------------------------------------------------------------------
    # Divide the lines into two classes,
    # where one has one direction of its lines
    # and the other has the opposite direction
    line_classifier = 0
    # Add all lines
    for line_segments in list_of_points:
        # Allocate an empty list for the lines in this segment
        temp_lines = []
        # Loop over all points in our line segments
        # and combine them to lines
        for i in range(len(line_segments)):
            # Make sure that we loop until the second last point
            if i < len(line_segments)-1:
                # For the even ones we add lines in one direction
                if line_classifier%2 == 0:
                    # Append the new line to our list of lines
                    temp_lines.append((line_segments[i][4],line_segments[i+1][4],object_counter))
                else:
                    # Append the new line to our list of lines
                    temp_lines.append((line_segments[i+1][4],line_segments[i][4],object_counter))                
                # Update the counter for the lines
                object_counter += 1
        # Save all lines in our line list
        if line_classifier%2 == 0: # Even lines
            list_of_lines.append(temp_lines) # Standard direction
        else: # Odd lines
            list_of_lines.append(temp_lines[::-1]) # Reverse direction
        # Lastly, increment the line classifier
        line_classifier += 1
    # -------------------------------------------------------------------
    # Add lines between each arch and define curves
    # -------------------------------------------------------------------
    # Allocate memory for the inbetween lines
    in_between_lines = []
    # Loop through all points and add in between lines
    for arc_index in range(len(list_of_lines)-1):
        # Create a between index
        between_index = -1
        # Allocate a temporary list
        in_between_line = []
        #Loop over all lines in the arc
        for line_index in range(len(list_of_lines[arc_index])):
            # Start with the even lines
            if arc_index %2 == 0:
                # Take the standard surface elements in the middle    
                if phi_vec[line_index] != 0 and phi_vec[line_index] != math.pi: # Normal case
                    # Extract the first point on curve 1
                    point_1_on_curve_1 = list_of_lines[arc_index][line_index][0]
                    # Extract the second point on curve 1
                    point_2_on_curve_1 = list_of_lines[arc_index][line_index][1]            
                    # Extract the first point on curve 2
                    point_1_on_curve_2 = list_of_lines[arc_index+1][len(list_of_lines[arc_index])-(line_index+1)][0]
                    # Extract the second point on curve 2
                    point_2_on_curve_2 = list_of_lines[arc_index+1][len(list_of_lines[arc_index])-(line_index+1)][1]            
                    # Save the first in between line as a 3-tuple
                    in_between_line.append((point_2_on_curve_1,point_1_on_curve_2, object_counter))
                    # Increment the line_counter
                    object_counter += 1
                    # Save the second in between line as a 3-tuple
                    in_between_line.append((point_2_on_curve_2,point_1_on_curve_1, object_counter))
                    # Increment the line_counter
                    object_counter += 1
                    # Define a mini-curve
                    mini_curve = [list_of_lines[arc_index][line_index][2], object_counter-2,list_of_lines[arc_index+1][len(list_of_lines[arc_index])-(line_index+1)][2], object_counter-1]
                    # Append the curve at hand to the list_of_surfaces
                    list_of_curve_loops.append([mini_curve,object_counter])
                    # Increment the line_counter
                    object_counter += 1           
                elif phi_vec[line_index] == 0:# Special case, first surface element
                    # Extract the second point on curve 1
                    point_2_on_curve_1 = list_of_lines[arc_index][line_index][1]            
                    # Extract the first point on curve 2
                    point_1_on_curve_2 = list_of_lines[arc_index+1][len(list_of_lines[arc_index])-(line_index+1)][0]
                    # Save the second in between line as a 3-tuple
                    in_between_line.append((point_2_on_curve_1,point_1_on_curve_2, object_counter))
                    # Increment the object_counter
                    object_counter += 1                
                    # Define a mini-curve
                    mini_curve = [list_of_lines[arc_index][line_index][2], object_counter-1,list_of_lines[arc_index+1][len(list_of_lines[arc_index])-(line_index+1)][2]]
                    # Append the curve at hand to the list_of_curve_loops
                    list_of_curve_loops.append([mini_curve, object_counter])
                    # Increment the object_counter
                    object_counter += 1                                        
                elif phi_vec[line_index] == math.pi:
                    # Extract the second point on curve 1
                    point_1_on_curve_1 = list_of_lines[arc_index][line_index][0]            
                    # Extract the first point on curve 2
                    point_2_on_curve_2 = list_of_lines[arc_index+1][len(list_of_lines[arc_index])-(line_index+1)][1]
                    # Save the second in between line as a 3-tuple
                    in_between_line.append((point_2_on_curve_2,point_1_on_curve_1, object_counter))
                    # Increment the object_counter
                    object_counter += 1                
                    # Define a mini-curve
                    mini_curve = [list_of_lines[arc_index][line_index][2],list_of_lines[arc_index+1][len(list_of_lines[arc_index])-(line_index+1)][2],object_counter-1]
                    # Append the curve at hand to the list_of_curve_loops
                    list_of_curve_loops.append([mini_curve, object_counter])
                    # Increment the object_counter
                    object_counter += 1                                        
                # Increment the between_index
                between_index += 2
            else: # Odd lines
                # Take the standard surface elements in the middle    
                if phi_vec[line_index] != 0 and phi_vec[line_index] != math.pi: # Normal case
                    # Extract the first point on curve 2
                    point_1_on_curve_2 = list_of_lines[arc_index+1][line_index][0]
                    # Extract the second point on curve 2
                    point_2_on_curve_2 = list_of_lines[arc_index+1][line_index][1]
                    # Extract the second point on curve 1
                    point_2_on_curve_1 = list_of_lines[arc_index][len(list_of_lines[arc_index])-(line_index+1)][1]   
                    # Extract the first point on curve 1
                    point_1_on_curve_1 = list_of_lines[arc_index][len(list_of_lines[arc_index])-(line_index+1)][0]
                    # Save the first in between line as a 3-tuple
                    in_between_line.append((point_2_on_curve_1,point_1_on_curve_2, object_counter))
                    # Increment the object_counter
                    object_counter += 1
                    # Save the second in between line as a 3-tuple
                    in_between_line.append((point_2_on_curve_2,point_1_on_curve_1, object_counter))
                    # Increment the object_counter
                    object_counter += 1
                    # Define a mini-curve
                    mini_curve = [list_of_lines[arc_index][len(list_of_lines[arc_index])-(line_index+1)][2], object_counter-1,list_of_lines[arc_index+1][line_index][2], object_counter-2]
                    # Append the curve at hand to the list_of_curve_loops
                    list_of_curve_loops.append([mini_curve, object_counter])
                    # Increment the object_counter
                    object_counter += 1                                        
                elif phi_vec[line_index] == 0:# Special case, first surface element
                    # Extract the second point on curve 1
                    point_1_on_curve_1 = list_of_lines[arc_index][len(list_of_lines[arc_index])-(line_index+1)][0]            
                    # Extract the first point on curve 2
                    point_2_on_curve_2 = list_of_lines[arc_index+1][line_index][1]
                    # Save the second in between line as a 3-tuple
                    in_between_line.append((point_2_on_curve_2,point_1_on_curve_1, object_counter))
                    # Increment the object_counter
                    object_counter += 1                
                    # Define a mini-curve
                    mini_curve = [list_of_lines[arc_index+1][line_index][2], object_counter-1,list_of_lines[arc_index][len(list_of_lines[arc_index])-(line_index+1)][2]]
                    # Append the curve at hand to the list_of_curve_loops
                    list_of_curve_loops.append([mini_curve, object_counter])
                    # Increment the object_counter
                    object_counter += 1                                        
                elif phi_vec[line_index] == math.pi:
                    # Extract the second point on curve 1
                    point_2_on_curve_1 = list_of_lines[arc_index][len(list_of_lines[arc_index])-(line_index+1)][1]            
                    # Extract the first point on curve 2
                    point_1_on_curve_2 = list_of_lines[arc_index+1][line_index][0]
                    # Save the second in between line as a 3-tuple
                    in_between_line.append((point_2_on_curve_1,point_1_on_curve_2, object_counter))
                    # Increment the object_counter
                    object_counter += 1                
                    # Define a mini-curve
                    mini_curve = [list_of_lines[arc_index+1][line_index][2],list_of_lines[arc_index][len(list_of_lines[arc_index])-(line_index+1)][2],object_counter-1]
                    # Append the curve at hand to the list_of_curve_loops
                    list_of_curve_loops.append([mini_curve,object_counter])
                    # Increment the object_counter 
                    object_counter += 1                                        
                # Increment the between_index
                between_index += 2                
            #--------------------------------------------------------------------
            # End of special cases
            #--------------------------------------------------------------------
        # Save all in-between lines
        in_between_lines.append(in_between_line)
    #--------------------------------------------------------------------
    # Final manipulation of lines and surfaces
    #--------------------------------------------------------------------
    # Merge the line we have at hand
    list_of_lines = list_of_lines + in_between_lines
    # Add all surfaces now
    for curve in list_of_curve_loops:
        # Save the surface
        list_of_surfaces.append([[curve[1]], object_counter])
        # Increment the object counter
        object_counter += 1
    #--------------------------------------------------------------------
    # Save all lists and return output
    #--------------------------------------------------------------------
    # Add the lines to the list of lists which lastly is returned
    list_of_lists.append(list_of_points)   # Points
    list_of_lists.append(list_of_lines)    # Lines
    list_of_lists.append(list_of_curve_loops) # Surfaces of rest if the sphere
    list_of_lists.append(list_of_surfaces) # Surfaces of rest if the sphere    
    # Also, make sure to return the new counter as well
    new_counter = object_counter
    return new_counter, list_of_lists
# -------------------------------------------------------------------
# END OF FUNCTIONS
# -------------------------------------------------------------------
#-----------------------------------------------------------------------------------        
# -------------------------------------------------------------------
# Allocating memory for all lists and vectors
# -------------------------------------------------------------------
# Define a centimetre
cm = 1e-02
#Lc1 = 0.01
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
# Define a thickness of the point
Lc1 = 0.01
# Add all centres to a list
list_of_centres = [c1, c2]
# Define a label counter
label_counter = 1
# Define the step size
#step_size = math.pi/20
step_size = math.pi/24
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
new_label, list_of_lists, theta_hole, phi_hole = add_circular_holes_on_sphere(list_of_centres,r,R,label_counter,step_size)
# Update the label counter
label_counter = new_label
#---------------------------------------------------------------
# Add the hole to the mesh
#---------------------------------------------------------------
# The label is currently at
# Add all points
for point in list_of_lists[0]:
    hole_on_sphere.addPoint(point[0], point[1], point[2], point[3], point[4])
# Add all lines
for line in list_of_lists[1]:
    hole_on_sphere.addLine(line[0], line[1], line[2])
# Add all circle arcs    
for circle_arc in list_of_lists[2]:
    hole_on_sphere.addCircleArc(circle_arc[0], circle_arc[1], circle_arc[2], circle_arc[3])
# Add all curve loops
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
# SEGMENT ABOVE HOLE
# Define phi
phi_below = (phi_hole[1],math.pi)
# Define theta
theta_below = theta_hole
# Caclulate the rest of the sphere
new_counter, list_of_lists_sphere_below = add_parts_of_sphere(R,phi_below,theta_below,step_size,label_counter,north_pole,south_pole)
# Update the label_counter
label_counter = new_counter
# SEGMENT BELOW HOLE
# Define phi
phi_above = (0,phi_hole[0])
# Define theta
theta_above = theta_hole
# Caclulate the rest of the sphere
new_counter, list_of_lists_sphere_above = add_parts_of_sphere(R,phi_above,theta_above,step_size,label_counter,north_pole,south_pole)
# Update the label_counter
label_counter = new_counter
# SEGMENT ON THE LEFT
# Define phi
phi_left = (0,math.pi)
# Define theta
theta_left = [-np.pi, theta_hole[0]]
# Print the dimensions of the hole
# Caclulate the rest of the sphere
new_counter, list_of_lists_sphere_left = add_parts_of_sphere(R,phi_left,theta_left,step_size,label_counter,north_pole,south_pole)
# Update the label_counter
label_counter = new_counter
# SEGMENT ON THE RIGHT
# Define phi
phi_right = (0,math.pi)
# Define theta
theta_right = [theta_hole[1],np.pi]
# Caclulate the rest of the sphere
new_counter, list_of_lists_sphere_right = add_parts_of_sphere(R,phi_right,theta_right,step_size,label_counter,north_pole,south_pole)
# Update the label_counter
label_counter = new_counter
# MERGE THE LISTS:
list_of_lists_sphere = []
# Loop over all lengths
for index in range(len(list_of_lists_sphere_above)):
    # Create a temporary list
    temp_list = list_of_lists_sphere_below[index] + list_of_lists_sphere_above[index]+ list_of_lists_sphere_left[index] + list_of_lists_sphere_right[index]
    #temp_list = list_of_lists_sphere_below[index] + list_of_lists_sphere_above[index]+list_of_lists_sphere_left[index]
    # Append the list
    list_of_lists_sphere.append(temp_list)

#---------------------------------------------------------------
# Add the rest of the sphere
#---------------------------------------------------------------
# Add all points
for arc in list_of_lists_sphere[0]:
    for point in arc:
        if point[4]!= 1 and point[4]!= 2:
            hole_on_sphere.addPoint(point[0], point[1], point[2], point[3], point[4])
# Add all lines
for arc in list_of_lists_sphere[1]:
    for line in arc:
        hole_on_sphere.addLine(line[0], line[1], line[2])
# Add all curve loops
for curve_loop in list_of_lists_sphere[2]:
    hole_on_sphere.addCurveLoop(curve_loop[0],curve_loop[1])
# Add all surfaces
for surface in list_of_lists_sphere[3]:
    hole_on_sphere.addPlaneSurface(surface[0],surface[1])    
# -------------------------------------------------------------------
# SYNCHRONIZE THE MESH
# -------------------------------------------------------------------
# Synchronise everything with Gmsh
hole_on_sphere.synchronize()
# -------------------------------------------------------------------
# SYNCHRONIZE THE MESH
# -------------------------------------------------------------------
# PHYSICAL GROUPS FOR CURVES:
circle_arc = gmsh.model.addPhysicalGroup(1,[circle_arc[3] for circle_arc in list_of_lists[2]])
gmsh.model.setPhysicalName(1,circle_arc,"Boundary curve for the hole")
# PHYSICAL GROUPS FOR SURFACES:
# The hole itself
surface_hole = gmsh.model.addPhysicalGroup(2,[hole[1] for hole in list_of_lists[4]])
gmsh.model.setPhysicalName(2,surface_hole,"Surface hole")
# Adjacent hole
surface_adjacent_hole = gmsh.model.addPhysicalGroup(2,[adjacent_hole[1] for adjacent_hole in list_of_lists[5]])
gmsh.model.setPhysicalName(2, surface_adjacent_hole, "Surface adjacent hole")
# The rest of the sphere
surface_rest_of_sphere = gmsh.model.addPhysicalGroup(2,[surface[1] for surface in list_of_lists_sphere[3]])
gmsh.model.setPhysicalName(2, surface_rest_of_sphere, "Surface the rest of the sphere")
# -------------------------------------------------------------------
# DEFINE COLOURS FOR THE MESH
# -------------------------------------------------------------------
# For nice colours, please see (https://colorbrewer2.org/)...
#-----------------------------------------------------------------------------
gmsh.model.setColor([(2,sphere[1]) for sphere in list_of_lists_sphere[3]], 2, 56, 88)  # Dark blue
gmsh.model.setColor([(2,hole[1]) for hole in list_of_lists[4]], 255, 247, 251)  # Lightest blue
gmsh.model.setColor([(1,circle_arc[3]) for circle_arc in list_of_lists[2]], 103, 0, 31)  # Light Magenta
#-----------------------------------------------------------------------------
#gmsh.model.setColor([(2,sphere[1]) for sphere in list_of_lists_sphere[3]], 54, 144, 192)  # Normal blue
#gmsh.model.setColor([(2,sphere[1]) for sphere in list_of_lists_sphere[3]], 2, 56, 88)  # Dark blue
#for arc in list_of_lists_sphere[1]:
    #gmsh.model.setColor([(1,line[2]) for line in arc], 54, 144, 192)  # Normal blue
# Define the colours for the holes
#gmsh.model.setColor([(2,hole[1]) for hole in list_of_lists[4]], 2, 56, 88)  # Dark blue
#gmsh.model.setColor([(2,hole[1]) for hole in list_of_lists[4]], 255, 247, 251)  # Lightest blue
#gmsh.model.setColor([(1,29), (1,30)], 166, 189, 219)  # Light blue
#gmsh.model.setColor([(2,31), (2,32)], 166, 189, 219)  # Light blue
#gmsh.model.setColor([(2,adjacent_hole[1]) for adjacent_hole in list_of_lists[5]], 208, 209, 230)  # Medium blue
#gmsh.model.setColor([(1,line[2]) for line in list_of_lists[1]], 255, 247, 251)  # Lightest blue
#gmsh.model.setColor([(1,circle_arc[3]) for circle_arc in list_of_lists[2]], 2, 56, 88)  # Lightest blue
#gmsh.model.setColor([(1,circle_arc[3]) for circle_arc in list_of_lists[2]], 103, 0, 31)  # Light Magenta
# -------------------------------------------------------------------
# GENERATING THE MESH
# -------------------------------------------------------------------
# Generate the two dimensional mesh (as we work with surfaces) 
gmsh.model.mesh.generate(2)
# Write the mesh to our file
gmsh.write("../Meshes/sphere_with_one_hole.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
