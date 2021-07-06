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
import numpy as np  # Needed to create arrays
from numpy import array # Import the array part to give an initial guess to the Newton solver 
from scipy.optimize import fmin # Import fmin to minimize a function
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
# Function 1: hole_on_sphere_locator
# The function takes in four inputs:
# 1. R the radius of the big sphere,
# 2. r the radius of the small sphere,
# 3. theta_0 the theta-coordinate of the centre of the small sphere,
# 4. phi_0 the phi-coordinate of the center of the small sphere.
# Given this it defines the function that fmin minimises in order
# to see the solutions of (theta, phi) that defines the border of
# the small circle that the intersection between these two spheres
# correspond to.
def hole_on_sphere_locator(t, params):
    # Step 1: Extract the parameters
    R = params[0]
    r = params[1]
    x_0 = params[2]
    y_0 = params[3]
    z_0 = params[4]
    # Extract the provided variables we will
    # optimise
    x, y, z=t
    # huge term
    term1 = (x0**2 + y0**2 + z0**2)+2*(x0*x+y0*y+z0*z) + R**2 - r**2
    term2 = x**2+y**2+z**2-R**2
    # Lastly we return our beloved objective function
    return abs(term1) + abs(term2)
# Function 2: rotate_point
# The function takes two inputs:
# 1. The point p= [x,y,z],
# 2. The angle theta.
# The function returns the rotated point, which has been rotated
# with an angle theta anticlockwise.
def rotate_point(point,theta,centre):
    # Extract the central point
    x0, y0 = centre
    #Extract the original coordinates
    x, y = point
    # Rotate the point 
    x_rot = x*np.cos(theta)-y*np.sin(theta) + x0*np.cos(theta)-y0*np.sin(theta) + x0
    y_rot = x*np.sin(theta)+y*np.cos(theta) + x0*np.sin(theta)+y0*np.cos(theta) + y0
    # Print something
    print("(x,y)\t=\t(%s,%s)"%(str(x),str(y)))
    print("(x_rot,y_rot)\t=\t(%s,%s)"%(str(x_rot),str(y_rot)))    
    # Assemble results in an array
    point_rotated = array([x_rot, y_rot])
    # Return the rotated point
    return point_rotated

def transform_plane(P,centre):
    x0, y0, z0 = centre
    z = x0*P[0] + P[1]*y0 + z0 - ((x0**2+y0**2)/(z0))
    transformed_point = (P[0],P[1],z)
    return transformed_point
# -------------------------------------------------------------------
# Allocating memory for all lists and vectors
# -------------------------------------------------------------------
# Define a centimetre
cm = 1e-02
#Lc1 = 0.01
Lc1 = 0.01
R = 1 # Big radius
r = 0.1*R # Small radius
# Define our centre of the small sphere:
theta_0 = math.pi/4
phi_0 = math.pi/4
# Define a parameter for the start guess
theta_epsilon = theta_0/2
phi_epsilon = phi_0/2 
# where the centre is determined by the spherical
# coordinate (r,theta_0,phi_0).
#---------------------------------------------------------
# CENTRE
#---------------------------------------------------------
#Begin by adding
# this point to the mesh
# Define the 3D point by its cartesian coordinates
#x0 = R * np.cos(theta_0) * np.sin(phi_0)
#y0 = R * np.sin(theta_0) * np.sin(phi_0)
#z0 = R * np.cos(phi_0)
# Define a point counter
#point_counter = 1
# Add the coordinate to the mesh
#hole_on_sphere.addPoint(x0, y0, z0, Lc1, point_counter)        
#print("Centre:")
#print((x0,y0,z0))
#---------------------------------------------------------
# FIRST POINT ON CIRCLE
#---------------------------------------------------------
# Define our parameters for the Newton solver
#params_1 = (R_small,x0,y0)
# Define our start guess for the Newton solver
#centre = array([x0, y0])
# Try to solve the minimisation problem
#x = fmin(hole_on_sphere_locator,centre,args=(params_1,))
# Calculate the x,y and z coordinates
# Define the 3D point by its cartesian coordinates
#z = z0-np.sqrt(R**2-(x[0]-x0)**2-(x[1]-y0)**2)
#z = z0-(1/z0)*((x0*(x[0]-x0))+(y0*(x[1]-y0)))
#print("First point:")
#print((x[0],x[1],z))
# Define a counter for the mesh
#point_counter += 1
# Add the coordinate to the mesh
#hole_on_sphere.addPoint(x[0], x[1], z, Lc1, point_counter)        
#---------------------------------------------------------
# Rotate the point twice
#---------------------------------------------------------
# Rotate the original point
#point_rotated = rotate_point(x,-np.pi/4,centre)
# Calculate the z-coordinate
#z_rotated = z0-np.sqrt(R**2-(point_rotated[0]-x0)**2-(point_rotated[1]-y0)**2)
#z_rotated = z0-(1/z0)*((x0*(point_rotated[0]-x0))+(y0*(point_rotated[1]-y0)))
# Define a counter for the mesh
#point_counter += 1
# Add the coordinate to the mesh
#hole_on_sphere.addPoint(point_rotated[0], point_rotated[1], z_rotated, Lc1, point_counter)
#print("Rotated point")
#print(point_rotated)
# Rotate the original point
#point_rotated = rotate_point(point_rotated,-np.pi/4,centre)
# Calculate the z-coordinate
#z_rotated = z0-np.sqrt(R**2-(point_rotated[0]-x0)**2-(point_rotated[1]-y0)**2)
#z_rotated = z0-(1/z0)*((x0*(point_rotated[0]-x0))+(y0*(point_rotated[1]-y0)))
# Define a counter for the mesh
#point_counter += 1
# Add the coordinate to the mesh
#hole_on_sphere.addPoint(point_rotated[0], point_rotated[1], z_rotated, Lc1, point_counter)
#print("Rotated point second")
#print(point_rotated)

#P1 =
#P1 = (0, 1, 0)
#P2 = (0, 0, 0)
#P3 = (1, 0, 0)

# Add the circle arc
#hole_on_sphere.addPoint(P1[0], P1[1], P1[2], Lc1, 1)
#hole_on_sphere.addPoint(P2[0], P2[1], P2[2], Lc1, 2)
#hole_on_sphere.addPoint(P3[0], P3[1], P3[2], Lc1, 3)
#hole_on_sphere.addCircleArc(1, 2, 3, 55)
# Add line from origin to point 1
#hole_on_sphere.addLine(2,1,4)
#hole_on_sphere.addLine(3,2,5)
# Define a curve loop using 55, 4 and 5
#hole_on_sphere.addCurveLoop([4,55,5],66)
#hole_on_sphere.addPlaneSurface([66],666)



#x0 = 0.75
#y0 = -0.5
#z0 = R**2-(x0**2+y0**2)

# Define our parameters for the Newton solver
#params_1 = (R,R_small,x0,y0,z0)
# Define our start guess for the Newton solver
#centre = array([x0, y0,z0])
# Try to solve the minimisation problem
#x = fmin(hole_on_sphere_locator,centre,args=(params_1,))
#print("Centre")
#print((x0,y0,z0))
#hole_on_sphere.addPoint(x0, y0, z0, Lc1, 1)
#print("Point on intersecting circle")
#print((x[0],x[1],x[2]))
#hole_on_sphere.addPoint(x[0], x[1], x[2], Lc1, 2)

#centre = array([x[0], x[1]])
#original_point = array([x[0], x[1]])
#point_rotated = rotate_point(original_point,-np.pi/4,centre)
#z_temp = R**2 - point_rotated[0]**2-point_rotated[1]**2
#x_1 = array([point_rotated[0], point_rotated[1], z_temp])
#hole_on_sphere.addPoint(x_1[0], x_1[1], x_1[2], Lc1, 3)
#hole_on_sphere.addCircleArc(1, 2, 3, 55)
#centre = array([z0, y0, z0])
#T1 = transform_plane(P1,centre)
#T2 = transform_plane(P2,centre)
#T3 = transform_plane(P3,centre)
#hole_on_sphere.addPoint(T1[0], T1[1], T1[2], Lc1, 6)
#hole_on_sphere.addPoint(T2[0], T2[1], T2[2], Lc1, 7)
#hole_on_sphere.addPoint(T3[0], T3[1], T3[2], Lc1, 8)
# Add line from origin to point 1
#hole_on_sphere.addLine(7,6,9)
#hole_on_sphere.addLine(8,7,10)
# Define circle arc
#hole_on_sphere.addCircleArc(6, 7, 8, 56)
# Define a curve loop using 55, 4 and 5
#hole_on_sphere.addCurveLoop([9,56,10],77)
#hole_on_sphere.addPlaneSurface([77],777)
# -------------------------------------------------------------------
# DEFINE OUR CIRCLE AT HAND
# -------------------------------------------------------------------
# SPHERE 1: Centre
c1 = (0,0,0) # with radius R defined above
# SPHERE 1: Centre
x2 = 0.1
y2 = -0.2
z2 = math.sqrt(R**2-x2**2-y2**2)
c2 = (x2,y2,z2) # with radius R defined above
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
# Define the tangent ti
ti = np.cross([0, 0, 1],list(ni)) # Take cross product
ti = ti/np.linalg.norm(ti) # Normalise the same vector
# Define the bitangent bi
bi = np.cross(ti,ni)
# Print all this
print("Centre sphere 1 with radius r=%s:\tc1\t=\t%s"%(str(R),str(c1)))
print("Centre sphere 2 with radius r=%s\tc2\t=\t%s"%(str(r),str(c2)))
print("Distance between spheres:\td\t=\t%0.3f"%(d))
print("The second centre is on the big sphere as x^2+y^2+z^2=%0.3f and R^2=%0.3f\n"%(x2**2+y2**2+z2**2,R**2))
print("The weighting coefficient:\th\t=\t%0.3f"%(h))
print("The centre of the intersecting circle is:\tci\t=\t%s"%(str(ci)))
print("The radius of the intersecting circle:\tri\t=\t%0.3f"%(ri))
print("The normal of this plane is:\tni\t=\t%s"%(str(ni)))
print("The tangent of this plane is:\tti\t=\t%s"%(str(ti)))
print("The bi-tangent of this plane is:\tbi\t=\t%s"%(str(bi)))
# -------------------------------------------------------------------
# START ADDING POINTS
# -------------------------------------------------------------------
hole_on_sphere.addPoint(c1[0], c1[1], c1[2], Lc1, 1)
hole_on_sphere.addPoint(c2[0], c2[1], c2[2], Lc1, 2)
# Define the mid point of the circle
hole_on_sphere.addPoint(ci[0], ci[1], ci[2], Lc1, 3)
# Define first point on circle
theta = 0
p1 = tuple([ci[j] + ri*(ti[j]*np.cos(theta)+bi[j]*np.cos(theta)) for j in range(len(ci))])
hole_on_sphere.addPoint(p1[0], p1[1], p1[2], Lc1, 4)
# Define second point on circle
theta = math.pi
p2 = tuple([ci[j] + ri*(ti[j]*np.cos(theta)+bi[j]*np.cos(theta)) for j in range(len(ci))])
hole_on_sphere.addPoint(p2[0], p2[1], p2[2], Lc1, 5)
# Define two circle arcs between these three points
hole_on_sphere.addCircleArc(4, 3, 5, 56)
hole_on_sphere.addCircleArc(5, 3, 4, 57)
# Define two lines in order to make a surface
hole_on_sphere.addLine(4, 5, 6)
hole_on_sphere.addLine(5, 4, 7)
# Define two curve loops based on this
hole_on_sphere.addCurveLoop([56,7], 8)
hole_on_sphere.addCurveLoop([57,6], 9)
# Define two surfaces
hole_on_sphere.addPlaneSurface([8], 10)
hole_on_sphere.addPlaneSurface([9], 11)
print("----------------------------------")
print("The centre of the mini circle:\tci\t=\t%s"%(str(ci)))
print("The first point:\tp1\t=\t%s"%(str(p1)))
print("The second point:\tp2\t=\t%s"%(str(p2)))
# -------------------------------------------------------------------
# GENERATING THE MESH
# -------------------------------------------------------------------
# For nice colours, please see (https://colorbrewer2.org/)
# Synchronise everything with Gmsh
hole_on_sphere.synchronize()

v = gmsh.view.add("comments")
gmsh.view.addListDataString(v, [0, 0], ["Created with Gmsh"])
# Add a text string in model coordinates centered at (X,Y,Z) = (0, 0.11, 0),
# with some style attributes:
gmsh.view.addListDataString(v, [c1[0], c1[1], c1[2]], ["Centre of original sphere"],
                            ["Align", "Center", "Font", "Times-BoldItalic", "FontSize","30"])
gmsh.view.addListDataString(v, [c2[0], c2[1], c2[2]], ["Centre of smaller sphere"],
                            ["Align", "Center", "Font", "Times-BoldItalic","FontSize","30"])
# Set the colours of the surfaces
gmsh.model.setColor([(1,56),(1,57),(1,6),(1,7)], 221, 28, 119)  # Magenta
#gmsh.model.setColor([(1,57)], 201, 148, 199)  # Magenta, sligthly lesser
gmsh.model.setColor([(2,10), (2,11)], 103, 0, 31)  # Magenta, sligthly lesser


gmsh.model.mesh.generate(2)
#----------------------------------
# Generate the two dimensional mesh (as we work with surfaces) 
#gmsh.model.mesh.generate(2)
# Write the mesh to our file
#gmsh.write("../Meshes/sphere.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
