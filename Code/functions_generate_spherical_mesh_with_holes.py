# =================================================================================
# =================================================================================
# Script:"generate_spherical_mesh"
# Date: 2021-08-25
# Implemented by: Johannes Borgqvist
# Description:
# The script generate a spherical FEM-mesh with
# a hole on the surface. It uses the Python interface
# of Gmsh to generate the mesh and the script has been
# written using functions which automates the generation
# of the holes and the adjacent region around the holes
# as well as the rest of the sphere. The sphere that the
# mesh works with is the unit sphere (i.e. r=1). 
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
# Functions for automating the generation of the mesh
# =================================================================================
# =================================================================================
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
def add_circular_holes_on_sphere(list_of_centres,r,R,label_counter,main_centre):
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
    list_of_circle_arcs_hole = [] # Circle arcs for the hole
    list_of_circle_arcs_adjacent = [] # Circle arcs for the adjacent
    list_of_curve_loops = [] # Curve loops hole
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
    list_of_circle_arcs_hole.append((label_counter-3,label_counter-4,label_counter-1,label_counter+1)) # Circle arc 1
    list_of_circle_arcs_hole.append((label_counter-1,label_counter-4,label_counter-2,label_counter+2)) # Circle arc 2
    list_of_circle_arcs_hole.append((label_counter-2,label_counter-4,label_counter,label_counter+3)) # Circle arc 3
    list_of_circle_arcs_hole.append((label_counter,label_counter-4,label_counter-3,label_counter+4)) # Circle arc 4    
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
    temp_list =[ci[j] + 1.45*ri*(ti[j]*math.cos(math.pi/4)+bi[j]*math.sin(math.pi/4)) for j in range(len(ci)-1)]
    temp_list.append(np.sqrt(R**2-temp_list[0]**2-temp_list[1]**2))
    p5 = tuple(temp_list)
    # Corner 2 of square    
    temp_list =[ci[j] + 1.45*ri*(ti[j]*math.cos(3*math.pi/4)+bi[j]*math.sin(3*math.pi/4)) for j in range(len(ci)-1)]
    temp_list.append(np.sqrt(R**2-temp_list[0]**2-temp_list[1]**2))    
    p6 = tuple(temp_list)
    # Corner 3 of square 
    temp_list =[ci[j] + 1.05*ri*(ti[j]*math.cos(5*math.pi/4)+bi[j]*math.sin(5*math.pi/4)) for j in range(len(ci)-1)]
    temp_list.append(np.sqrt(R**2-temp_list[0]**2-temp_list[1]**2))
    p7 = tuple(temp_list)
    # Corner 4 of square 
    temp_list =[ci[j] + 1.05*ri*(ti[j]*math.cos(7*math.pi/4)+bi[j]*math.sin(7*math.pi/4)) for j in range(len(ci)-1)]
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
    #z = np.sqrt(R**2-x**2-y**2)    
    list_of_points.append((x, y, z, Lc1, label_counter)) # Add the point
    label_counter += 1    # Increment label counter
    upper_centre = (0,0,z)
    # Point 2
    x = R * np.cos(theta_tuple[0]) * np.sin(phi_tuple[1])
    y = R * np.sin(theta_tuple[0]) * np.sin(phi_tuple[1])
    z = R * np.cos(phi_tuple[1])
    #z = np.sqrt(R**2-x**2-y**2)
    list_of_points.append((x, y, z, Lc1, label_counter)) # Add the point
    label_counter += 1    # Increment label counter
    lower_centre = (0,0,z)    
    # Point 3
    x = R * np.cos(theta_tuple[1]) * np.sin(phi_tuple[1])
    y = R * np.sin(theta_tuple[1]) * np.sin(phi_tuple[1])
    z = R * np.cos(phi_tuple[1])
    #z = np.sqrt(R**2-x**2-y**2)    
    list_of_points.append((x, y, z, Lc1, label_counter)) # Add the point
    label_counter += 1    # Increment label counter
    # Point 4
    x = R * np.cos(theta_tuple[1]) * np.sin(phi_tuple[0])
    y = R * np.sin(theta_tuple[1]) * np.sin(phi_tuple[0])
    z = R * np.cos(phi_tuple[0])
    #z = np.sqrt(R**2-x**2-y**2)    
    list_of_points.append((x, y, z, Lc1, label_counter)) # Add the point
    label_counter += 1    # Increment label counter
    # Save the corner points
    corner_points = [label_counter-4,label_counter-3,label_counter-2,label_counter-1]
    list_of_points.append((upper_centre[0], upper_centre[1], upper_centre[2], Lc1, label_counter)) # Add the upper centre
    upper_centre_label = label_counter
    label_counter += 1    # Increment label counter
    list_of_points.append((lower_centre[0], lower_centre[1], lower_centre[2], Lc1, label_counter)) # Add the upper centre
    lower_centre_label = label_counter
    label_counter += 1    # Increment label counter    
    # Define circle arcs for the outer surface:
    # ARC 1
    arc1 = [corner_points[0], main_centre[4], corner_points[1], label_counter]
    label_counter +=1
    list_of_circle_arcs_adjacent.append(arc1)
    # ARC 2
    arc2 = [corner_points[1], lower_centre_label, corner_points[2], label_counter]
    label_counter +=1
    list_of_circle_arcs_adjacent.append(arc2)
    # ARC 3
    arc3 = [corner_points[2], main_centre[4], corner_points[3], label_counter]
    label_counter +=1
    list_of_circle_arcs_adjacent.append(arc3)
    # ARC 4
    arc4 = [corner_points[3], upper_centre_label, corner_points[0], label_counter]
    label_counter +=1
    list_of_circle_arcs_adjacent.append(arc4)
    # Save all these arcs in order to construct the outer surface
    outer_boundaries = [arc1, arc2, arc3, arc4]
    # Add a curve loop for all our lines
    list_of_curve_loops.append([[arc[3] for arc in outer_boundaries],label_counter+1])
    # Define a surface between our circle and the area around
    list_of_surfaces_adjacent.append([[curve[1] for curve in list_of_curve_loops],label_counter+1])        
    # Update the list, we return
    list_of_lists.append(list_of_points) # Points
    #list_of_lists.append(list_of_lines) # Lines
    list_of_lists.append(list_of_circle_arcs_hole) # Circle arcs for the hole
    list_of_lists.append(list_of_circle_arcs_adjacent) # Circle arcs for the adjacent region
    list_of_lists.append(list_of_curve_loops) # Curve loops hole
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
def add_parts_of_sphere(R,phi_tuple,theta_tuple,object_counter,north_pole,south_pole,main_centre):
    # Define the thickness of one point
    Lc1 = 0.01
    #-----------------------------------------------------------------------
    # Allocate memory for the ouput
    #-----------------------------------------------------------------------    
    # Allocate memory for the outputs (what we will return).
    list_of_lists = [] # Points, curves and surfaces
    # Allocate memory for the lists that will be included in the    
    # larger list of lists
    list_of_points = [] # Points
    list_of_circle_arcs = [] # Circle arcs
    list_of_curve_loops = [] # Curve loops 
    list_of_surfaces = [] # Surfaces of rest if the sphere
    # OKAY, IT IS NOT MAYBE THE NICEST SOLUTION BUT IT IS CORRECT.
    # SO LET'S DEFINE OUR FOUR CASES AND DO THEM ONE BY ONE.
    # CASE1: Standard surface with four sides
    if phi_tuple[0]!=0 and phi_tuple[1]!=math.pi:
        # Allocate memory for our six arcs
        arc_1 = []
        arc_2 = []
        arc_3 = []
        arc_4 = []
        arc_5 = []
        arc_6 = []
        # Define the mid point for the angle phi
        phi_midpoint = ((phi_tuple[1]+phi_tuple[0])/(2))
        # Define our six points
        # POINT 0:
        x = R * np.cos(theta_tuple[0]) * np.sin(phi_tuple[0])
        y = R * np.sin(theta_tuple[0]) * np.sin(phi_tuple[0])
        z = R * np.cos(phi_tuple[0])
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_1.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter
        # Save the z-coordinate in order to create a lower centre
        z_upper = z        
        # POINT 1:
        x = R * np.cos(theta_tuple[0]) * np.sin(phi_midpoint)
        y = R * np.sin(theta_tuple[0]) * np.sin(phi_midpoint)
        z = R * np.cos(phi_midpoint)
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_1.append(object_counter) # Add the point to our circle arc as well
        arc_2.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter        
        # POINT 2:
        x = R * np.cos(theta_tuple[0]) * np.sin(phi_tuple[1])
        y = R * np.sin(theta_tuple[0]) * np.sin(phi_tuple[1])
        z = R * np.cos(phi_tuple[1])
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_2.append(object_counter) # Add the point to our circle arc as well
        arc_3.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter
        # Save the z-coordinate in order to create a lower centre
        z_lower = z
        # POINT 3:
        x = R * np.cos(theta_tuple[1]) * np.sin(phi_tuple[1])
        y = R * np.sin(theta_tuple[1]) * np.sin(phi_tuple[1])
        z = R * np.cos(phi_tuple[1])
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_3.append(object_counter) # Add the point to our circle arc as well
        arc_4.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter
        # POINT 4:
        x = R * np.cos(theta_tuple[1]) * np.sin(phi_midpoint)
        y = R * np.sin(theta_tuple[1]) * np.sin(phi_midpoint)
        z = R * np.cos(phi_midpoint)
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_4.append(object_counter) # Add the point to our circle arc as well
        arc_5.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter
        # POINT 5:
        x = R * np.cos(theta_tuple[1]) * np.sin(phi_tuple[0])
        y = R * np.sin(theta_tuple[1]) * np.sin(phi_tuple[0])
        z = R * np.cos(phi_tuple[0])
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_5.append(object_counter) # Add the point to our circle arc as well
        arc_6.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter
        # Add the last point to the sixth arc
        arc_6.append(object_counter-6)
        # Add the main centre to arcs 1,2,4 and 5
        # Add smaller centres as well
        list_of_points.append((0, 0, z_upper, Lc1, object_counter)) # Add the point to the mesh
        upper_centre_counter = object_counter
        object_counter+=1
        list_of_points.append((0, 0, z_lower, Lc1, object_counter)) # Add the point to the mesh
        lower_centre_counter = object_counter
        object_counter+=1
        arc_1.append(main_centre[4]) # Add main centre arc 1
        arc_1.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_2.append(main_centre[4]) # Add main centre arc 1
        arc_2.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_3.append(lower_centre_counter) # Add main centre arc 1
        arc_3.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_4.append(main_centre[4]) # Add main centre arc 1
        arc_4.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_5.append(main_centre[4]) # Add main centre arc 1
        arc_5.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_6.append(upper_centre_counter) # Add main centre arc 1
        arc_6.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        # Add all these arcs to the list of circle arcs
        list_of_circle_arcs.append(arc_1)
        list_of_circle_arcs.append(arc_2)
        list_of_circle_arcs.append(arc_3)
        list_of_circle_arcs.append(arc_4)
        list_of_circle_arcs.append(arc_5)
        list_of_circle_arcs.append(arc_6)
        # Add a curve loop for all our circle arcs
        object_counter += 1 # Increment label
        list_of_curve_loops.append([[arcs[3] for arcs in list_of_circle_arcs],object_counter])
        object_counter += 1 # Increment label
        # Define a surface between our circle and the area around
        list_of_surfaces.append([[curve[1] for curve in list_of_curve_loops],object_counter])
        object_counter += 1 # Increment label        
    # CASE2: with only three sides and lower has a pointy corner
    elif phi_tuple[0]==0 and phi_tuple[1]!=math.pi:
        # Allocate memory for our five arcs
        arc_1 = []
        arc_2 = []
        arc_3 = []
        arc_4 = []
        arc_5 = []
        arc_6 = []
        # Define the mid point for the angle phi
        phi_midpoint = ((phi_tuple[1]+phi_tuple[0])/(2))
        # Define our six points
        # POINT 0:
        list_of_points.append((north_pole[0], north_pole[1], north_pole[2], north_pole[3], north_pole[4])) # Add the point to the mesh
        arc_1.append(north_pole[4]) # Add the point to our circle arc as well
        # POINT 1:
        x = R * np.cos(theta_tuple[0]) * np.sin(phi_midpoint)
        y = R * np.sin(theta_tuple[0]) * np.sin(phi_midpoint)
        z = R * np.cos(phi_midpoint)
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_1.append(object_counter) # Add the point to our circle arc as well
        arc_2.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter        
        # POINT 2:
        x = R * np.cos(theta_tuple[0]) * np.sin(phi_tuple[1])
        y = R * np.sin(theta_tuple[0]) * np.sin(phi_tuple[1])
        z = R * np.cos(phi_tuple[1])
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_2.append(object_counter) # Add the point to our circle arc as well
        arc_3.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter
        # Save the z-coordinate in order to create a lower centre
        z_lower = z
        # POINT 3:
        x = R * np.cos(theta_tuple[1]) * np.sin(phi_tuple[1])
        y = R * np.sin(theta_tuple[1]) * np.sin(phi_tuple[1])
        z = R * np.cos(phi_tuple[1])
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_3.append(object_counter) # Add the point to our circle arc as well
        arc_4.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter
        # POINT 4:
        x = R * np.cos(theta_tuple[1]) * np.sin(phi_midpoint)
        y = R * np.sin(theta_tuple[1]) * np.sin(phi_midpoint)
        z = R * np.cos(phi_midpoint)
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_4.append(object_counter) # Add the point to our circle arc as well
        arc_5.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter
        # Add a last point to arc_5 which is the north pole
        arc_5.append(north_pole[4])
        # Add the main centre to arcs 1,2,4 and 5
        list_of_points.append((0, 0, z_lower, Lc1, object_counter)) # Add the point to the mesh
        lower_centre_counter = object_counter
        object_counter+=1
        arc_1.append(main_centre[4]) # Add main centre arc 1
        arc_1.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_2.append(main_centre[4]) # Add main centre arc 1
        arc_2.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_3.append(lower_centre_counter) # Add main centre arc 1
        arc_3.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_4.append(main_centre[4]) # Add main centre arc 1
        arc_4.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_5.append(main_centre[4]) # Add main centre arc 1
        arc_5.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        # Add all these arcs to the list of circle arcs
        list_of_circle_arcs.append(arc_1)
        list_of_circle_arcs.append(arc_2)
        list_of_circle_arcs.append(arc_3)
        list_of_circle_arcs.append(arc_4)
        list_of_circle_arcs.append(arc_5)
        # Add a curve loop for all our circle arcs
        object_counter += 1 # Increment label
        list_of_curve_loops.append([[arcs[3] for arcs in list_of_circle_arcs],object_counter])
        object_counter += 1 # Increment label
        # Define a surface between our circle and the area around
        list_of_surfaces.append([[curve[1] for curve in list_of_curve_loops],object_counter])
        object_counter += 1 # Increment label        
    # CASE3: with only three sides and upper has a pointy corner
    elif phi_tuple[0]!=0 and phi_tuple[1]==math.pi:
        # Allocate memory for our five arcs
        arc_1 = []
        arc_2 = []
        arc_4 = []
        arc_5 = []
        arc_6 = []
        # Define the mid point for the angle phi
        phi_midpoint = ((phi_tuple[1]+phi_tuple[0])/(2))
        # Define our six points
        # POINT 0:
        x = R * np.cos(theta_tuple[0]) * np.sin(phi_tuple[0])
        y = R * np.sin(theta_tuple[0]) * np.sin(phi_tuple[0])
        z = R * np.cos(phi_tuple[0])
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_1.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter
        # Save the z-coordinate in order to create a lower centre
        z_upper = z        
        # POINT 1:
        x = R * np.cos(theta_tuple[0]) * np.sin(phi_midpoint)
        y = R * np.sin(theta_tuple[0]) * np.sin(phi_midpoint)
        z = R * np.cos(phi_midpoint)
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_1.append(object_counter) # Add the point to our circle arc as well
        arc_2.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter        
        # POINT 2:
        list_of_points.append((south_pole[0], south_pole[1], south_pole[2], south_pole[3], south_pole[4])) # Add the point to the mesh
        arc_2.append(south_pole[4]) # Add the point to our circle arc as well
        arc_4.append(south_pole[4]) # Add the point to our circle arc as well        
        # POINT 4:
        x = R * np.cos(theta_tuple[1]) * np.sin(phi_midpoint)
        y = R * np.sin(theta_tuple[1]) * np.sin(phi_midpoint)
        z = R * np.cos(phi_midpoint)
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_4.append(object_counter) # Add the point to our circle arc as well
        arc_5.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter
        # POINT 5:
        x = R * np.cos(theta_tuple[1]) * np.sin(phi_tuple[0])
        y = R * np.sin(theta_tuple[1]) * np.sin(phi_tuple[0])
        z = R * np.cos(phi_tuple[0])
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_5.append(object_counter) # Add the point to our circle arc as well
        arc_6.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter
        # Add the last point to the sixth arc
        arc_6.append(object_counter-4)
        # Add the main centre to arcs 1,2,4 and 5
        # Add smaller centres as well
        list_of_points.append((0, 0, z_upper, Lc1, object_counter)) # Add the point to the mesh
        upper_centre_counter = object_counter
        object_counter+=1
        arc_1.append(main_centre[4]) # Add main centre arc 1
        arc_1.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_2.append(main_centre[4]) # Add main centre arc 1
        arc_2.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_4.append(main_centre[4]) # Add main centre arc 1
        arc_4.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_5.append(main_centre[4]) # Add main centre arc 1
        arc_5.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_6.append(upper_centre_counter) # Add main centre arc 1
        arc_6.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        # Add all these arcs to the list of circle arcs
        list_of_circle_arcs.append(arc_1)
        list_of_circle_arcs.append(arc_2)
        list_of_circle_arcs.append(arc_4)
        list_of_circle_arcs.append(arc_5)
        list_of_circle_arcs.append(arc_6)
        # Add a curve loop for all our circle arcs
        object_counter += 1 # Increment label
        list_of_curve_loops.append([[arcs[3] for arcs in list_of_circle_arcs],object_counter])
        object_counter += 1 # Increment label
        # Define a surface between our circle and the area around
        list_of_surfaces.append([[curve[1] for curve in list_of_curve_loops],object_counter])
        object_counter += 1 # Increment label        
    # CASE4: with two sides since corners are at phi=0 and phi=pi respectively
    elif phi_tuple[0]==0 and phi_tuple[1]==math.pi:
        # Allocate memory for our six arcs
        arc_1 = []
        arc_2 = []
        arc_4 = []
        arc_5 = []
        # Define the mid point for the angle phi
        phi_midpoint = ((phi_tuple[1]+phi_tuple[0])/(2))
        # Define our six points
        # POINT 0:
        list_of_points.append((north_pole[0], north_pole[1], north_pole[2], north_pole[3], north_pole[4])) # Add the point to the mesh
        arc_1.append(north_pole[4]) # Add the point to our circle arc as well
        # POINT 1:
        x = R * np.cos(theta_tuple[0]) * np.sin(phi_midpoint)
        y = R * np.sin(theta_tuple[0]) * np.sin(phi_midpoint)
        z = R * np.cos(phi_midpoint)
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_1.append(object_counter) # Add the point to our circle arc as well
        arc_2.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter        
        # POINT 2:
        list_of_points.append((south_pole[0], south_pole[1], south_pole[2], south_pole[3], south_pole[4])) # Add the point to the mesh
        arc_2.append(south_pole[4]) # Add the point to our circle arc as well
        arc_4.append(south_pole[4]) # Add the point to our circle arc as well
        # POINT 4:
        x = R * np.cos(theta_tuple[1]) * np.sin(phi_midpoint)
        y = R * np.sin(theta_tuple[1]) * np.sin(phi_midpoint)
        z = R * np.cos(phi_midpoint)
        list_of_points.append((x, y, z, Lc1, object_counter)) # Add the point to the mesh
        arc_4.append(object_counter) # Add the point to our circle arc as well
        arc_5.append(object_counter) # Add the point to our circle arc as well
        object_counter += 1 # Increment the label counter
        # POINT 5:
        arc_5.append(north_pole[4]) # Add the point to our circle arc as well
        # Add the main centre to arcs 1,2,4 and 5
        arc_1.append(main_centre[4]) # Add main centre arc 1
        arc_1.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_2.append(main_centre[4]) # Add main centre arc 1
        arc_2.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_4.append(main_centre[4]) # Add main centre arc 1
        arc_4.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        arc_5.append(main_centre[4]) # Add main centre arc 1
        arc_5.append(object_counter) # Add label to this arc 
        object_counter += 1 # Increment label
        # Add all these arcs to the list of circle arcs
        list_of_circle_arcs.append(arc_1)
        list_of_circle_arcs.append(arc_2)
        list_of_circle_arcs.append(arc_4)
        list_of_circle_arcs.append(arc_5)
        # Add a curve loop for all our circle arcs
        object_counter += 1 # Increment label
        list_of_curve_loops.append([[arcs[3] for arcs in list_of_circle_arcs],object_counter])
        object_counter += 1 # Increment label
        # Define a surface between our circle and the area around
        list_of_surfaces.append([[curve[1] for curve in list_of_curve_loops],object_counter])
        object_counter += 1 # Increment label        
     #--------------------------------------------------------------------
    # Save all lists and return output
    #--------------------------------------------------------------------
    # Add the lines to the list of lists which lastly is returned
    list_of_lists.append(list_of_points)   # Points
    list_of_lists.append(list_of_circle_arcs)    # Lines
    list_of_lists.append(list_of_curve_loops) # Surfaces of rest if the sphere
    list_of_lists.append(list_of_surfaces) # Surfaces of rest if the sphere    
    # Also, make sure to return the new counter as well
    new_counter = object_counter
    return new_counter, list_of_lists
#-----------------------------------------------------------------------------------
# Function 3: "add_sphere_without_hole"
# The function generates the mesh for parts of the sphere
# with centre in the origin (0,0,0) and radius R. It calls
# Function 2 "add_parts_of_sphere" repeatedly and it generates
# slices of the sphere for small intervals in the theta direction.
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
def add_sphere_without_hole(R,phi_tuple,theta_tuple,object_counter,north_pole,south_pole,main_centre):
    # Define a step size in the theta direction
    #theta_step = math.pi/5
    theta_step = math.pi/24
    #theta_step = math.pi/75
    #theta_step = math.pi/100
    #theta_step = math.pi/125
    #theta_step = math.pi/150
    #theta_step = math.pi/175
    #theta_step = math.pi/200
    #theta_step = math.pi/300
    #theta_step = math.pi/400
    #theta_step = math.pi/500
    #theta_step = math.pi/350
    #theta_step = math.pi/360
    # Define the number of steps based on this step size
    num_theta_steps = int(round(((theta_tuple[1]-theta_tuple[0])/(theta_step))))
    # Define our theta vector
    theta_vec = np.linspace(theta_tuple[0],theta_tuple[1],num_theta_steps)
    # Allocate memory for the list of lists we will return
    list_of_lists_sphere = []
    list_of_points = []
    list_of_circle_arcs = []
    list_of_curve_loops = []
    list_of_surfaces = []
    # Append all the allocated smaller lists to the main list
    list_of_lists_sphere.append(list_of_points)   # Points
    list_of_lists_sphere.append(list_of_circle_arcs)    # Lines
    list_of_lists_sphere.append(list_of_curve_loops) # Surfaces of rest if the sphere
    list_of_lists_sphere.append(list_of_surfaces) # Surfaces of rest if the sphere        
    # Loop over all angles and calculates the slice at hand
    for index in range(len(theta_vec)-1):
        # Define the test vector for the current slice
        theta_test = [theta_vec[index], theta_vec[index+1]]
        # Call the function for the small interval at hand defined by our theta vector
        new_counter, list_of_lists_temp = add_parts_of_sphere(R,phi_tuple,theta_test,object_counter,north_pole,south_pole,main_centre)
        # Update the object counter
        object_counter = new_counter+1
        # Add all temporary lists to the main list which we return
        list_of_lists_sphere[0].append(list_of_lists_temp[0]) # Add the points for the newest slice
        list_of_lists_sphere[1].append(list_of_lists_temp[1]) # Add the circle arcs for the newest slice
        list_of_lists_sphere[2].append(list_of_lists_temp[2]) # Add the curve loops for the newest slice
        list_of_lists_sphere[3].append(list_of_lists_temp[3]) # Add the surfaces for the newest slice        
    # It is not certain that the step size in the linspace vector matches exactly with the provided
    # interval so if needed we add an extra last slice that is potentially thinner than the provided
    # step size.
    if theta_vec[len(theta_vec)-1] < theta_tuple[1]:
        theta_test = [theta_vec[len(theta_vec)-1], theta_tuple[1]]
        # Call the function for the small interval at hand defined by our theta vector
        new_counter, list_of_lists_temp = add_parts_of_sphere(R,phi_tuple,theta_test,object_counter+1,north_pole,south_pole,main_centre)
        # Update the object counter
        object_counter = new_counter+1
        # Add all temporary lists to the main list which we return
        list_of_lists_sphere[0].append(list_of_lists_temp[0]) # Add the points for the newest slice
        list_of_lists_sphere[1].append(list_of_lists_temp[1]) # Add the circle arcs for the newest slice
        list_of_lists_sphere[2].append(list_of_lists_temp[2]) # Add the curve loops for the newest slice
        list_of_lists_sphere[3].append(list_of_lists_temp[3]) # Add the surfaces for the newest slice               
    #  Lastly, return this massive list and return the updated counter
    new_counter = object_counter
    return new_counter, list_of_lists_sphere   
    
# -------------------------------------------------------------------
# END OF FUNCTIONS
# -------------------------------------------------------------------
