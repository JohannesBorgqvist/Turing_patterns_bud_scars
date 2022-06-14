# =================================================================================
# =================================================================================
# Script:"pairwise_comparison"
# Date: 2022-06-13
# Implemented by: Johannes Borgqvist
# Description:
# This script compares the properties of simulations when the hole is placed at
# (0,0,-1) and (1,0,0).
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import meshio # To extract the important parts of the mesh
import numpy as np # Import numpy as well
import Schnakenberg_properties # Home made
from sklearn.cluster import DBSCAN # To calculate the number of poles in the concentration profile
from matplotlib import pyplot as plt # For plotting
from fenics import *
import math as m
# =================================================================================
# =================================================================================
# Functions
# =================================================================================
# =================================================================================
def cart2sph(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r = m.sqrt(XsqPlusYsq + z**2)               # r
    elev = m.atan2(z,m.sqrt(XsqPlusYsq))     # theta
    az = m.atan2(y,x)                           # phi
    return r, elev, az
def plot_LaTeX_2D(t,y,file_str,plot_str,legend_str):
    # Open a file with the append option
    # so that we can write to the same
    # file multiple times
    f = open(file_str, "a")
    # Create a temporary string which
    # is the one that does the plotting.
    # Here we incorporate the input plot_str
    # which contains the color, and the markers
    # of the plot at hand
    if len(legend_str)==0:
        temp_str = "\\addplot[\nforget plot,\n" + plot_str+ "\n]\n"
    else:
        temp_str = "\\addplot[\n" + plot_str+ "\n]\n"
    # Add the coordinates
    temp_str += "coordinates {%\n"
    # Loop over the input files and add
    # them to the file
    for i in range(len(t)):
        temp_str += "(" + str(t[i]) + "," + str(y[i]) + ")\n"
    # The plotting is done, let's close the shop    
    temp_str += "};\n"
    # Add a legend if one is provided
    if len(legend_str) > 0:
        temp_str += "\\addlegendentry{" + legend_str + "}\n"
    # Finally, we write the huge string
    # we have created
    f.write("%s"%(temp_str))
    # Close the file
    f.close()
def read_mesh_Schnakenberg_sphere_with_holes(num_holes,radii_holes):
    # Allocate memory for the mesh and the mesh value collection
    mesh = Mesh()
    mvc_subdomains = MeshValueCollection("size_t", mesh, 2)
    # Define a mesh function taking the subdomains into account
    #mf_subdomains = 0
    # Allocate memory for a list containing all the integration
    # measures involved in the variational formulation
    dx_list = []    
    # Define the string of the mesh that we want to read
    mesh_str = "../Meshes/s_h_" + str(num_holes)
    # Define the string in which we read the mesh
    # depending on the number of holes
    if num_holes > 0:
        # Loop over the radii and add these to the string name of the mesh
        for radius in radii_holes:
            mesh_str += "_r_" + str(round(radius,3)).replace(".","p") + "_"
    # Now, we add the file format to the string as well
    mesh_str += ".xdmf"
    # And tidy the file name up in necessary
    mesh_str = mesh_str.replace("_.xdmf",".xdmf")
    #print(mesh_str)
    # Read in the mesh and the subdomains into these two variables
    with XDMFFile(mesh_str) as infile:
        infile.read(mesh)
        infile.read(mvc_subdomains, "name_to_read")
    # Define a mesh function taking the subdomains into account
    mf_subdomains = cpp.mesh.MeshFunctionSizet(mesh, mvc_subdomains)
    # Add our integration measure to the list of integration measure
    dx_list.append(Measure("dx", domain=mesh, subdomain_data=mf_subdomains, subdomain_id=1))
    # Lastly, return our mesh, the mesh value collection, the mesh
    # function and the integration measures.
    return mesh, mvc_subdomains, mf_subdomains, dx_list
    
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Define the parameter pairs
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The parameters in the Schnakenberg model
a = 0.2
b = 1.0
# The wavenumber k^2
n = 2
k_squared = n*(n+1)
# Calculate the steady states and the critical parameters
u_0, v_0, d_c, gamma_c = Schnakenberg_properties.calculate_steady_states_and_critical_parameters_Schnakenberg(a,b,k_squared)
# Save the steady states in a list
steady_states = [u_0,v_0]
# Set the value of the relative diffusion
d = 18.0
# Set the value of the reaction strength to its critical value
gamma = gamma_c
# Define the number of holes
num_holes = 0
# Define the radius
radii_holes = []
# Define that we have the ICs around the steady states
ICs_around_steady_states = True
# Define the perturbation in the initial conditions
sigma = 1e-4
# Define the end time for the simulations
T = 50
# Let's start with the zeroth repitition
repitition_index = 0
# Define a parameter epsilon for the clustering
epsilon = 1
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Mesh at (0,0,-1)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Read the msh file
conc_profile = meshio.read("../Output/north_pole/iteration_0/u000101.vtu")
# Extract the concentration profile
u = np.asarray(list(conc_profile.point_data.values())[0])
## Extract the spatial coordinates
spatial_coordinates = conc_profile.points
# Define a threshold concentration
threshold_concentration = 0.8*max(u)
# Maximum conc
max_conc_north = max(u)
# Save the spatial coordinates for which the concentration profile is above our threshold value
pole_coordinates = np.asarray([spatial_coordinate for index,spatial_coordinate in enumerate(spatial_coordinates) if u[index]>=threshold_concentration])
# CLUSTERING: conduct the density based scan using DBSCAN
db = DBSCAN(eps=epsilon, min_samples=1).fit(pole_coordinates)
cluster_labels = db.labels_
# Get the number of poles
num_poles_north = len(set(cluster_labels))
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Calculate the relative pole area by loading the old mesh and integrating
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Load the old mesh
mesh_old = Mesh("../Output/north_pole/iteration_0/final_timestep_mesh.xml")
gdim = mesh_old.geometry().dim()
# Define a function space on the old mesh
H_old = FunctionSpace(mesh_old, "P", 1)
# Load the old initial conditions on the old function space
u_old = Function(H_old, "../Output/north_pole/iteration_0/final_timestep_u.xml")
v_old = Function(H_old, "../Output/north_pole/iteration_0/final_timestep_v.xml")
# Find out the coordinates in the mesh
# Create a dofmap
dofmap = H_old.dofmap()
sphere_dofs = dofmap.dofs()
# Find all coordinates for the nodes in our mesh
coordinates = H_old.tabulate_dof_coordinates()
# Calculate the coordinates in the pole
pole_1_coordinates = []
pole_1_dofs = []
pole_2_coordinates = []
pole_2_dofs = []
non_pole_dofs = []
for index,coordinate in enumerate(coordinates):
    for sub_index,temp_coordinate in enumerate(pole_coordinates):
        if (coordinate[0]==temp_coordinate[0]) and ((coordinate[1]==temp_coordinate[1])) and ((coordinate[2]==temp_coordinate[2])):
            if cluster_labels[sub_index] == 0:
                pole_1_coordinates.append(coordinate)
                pole_1_dofs.append(sphere_dofs[index])
            elif cluster_labels[sub_index] == 1:
                pole_2_coordinates.append(coordinate)
                pole_2_dofs.append(sphere_dofs[index])                
            else:
                non_pole_dofs.append(sphere_dofs[index])                
# Create a function as well 
pole_1_function = Function(H_old)
pole_2_function = Function(H_old)
# Set this pole function to 1 in the pole and 0 outside
pole_1_function.vector()[non_pole_dofs] = 0
pole_1_function.vector()[pole_1_dofs] = 1
pole_2_function.vector()[non_pole_dofs] = 0
pole_2_function.vector()[pole_2_dofs] = 1
# Define an integration measure
dx = Measure("dx", domain=mesh_old)
# Integrate over domain
pole_1_area_integration_north = 100*assemble(pole_1_function*dx)/assemble(Constant(1.0)*dx)
pole_2_area_integration_north = 100*assemble(pole_2_function*dx)/assemble(Constant(1.0)*dx)    
relative_pole_area_north = pole_1_area_integration_north + pole_2_area_integration_north
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Calculate the relative pole area by loading the old mesh and integrating
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Calculate the centroids as just the centres of mass of the poles
centroid_1 = [np.mean(np.asarray([point[0] for index,point in enumerate(pole_1_coordinates)])),np.mean(np.asarray([point[1] for index,point in enumerate(pole_1_coordinates)])),np.mean(np.asarray([point[2] for index,point in enumerate(pole_1_coordinates)]))]
centroid_2 = [np.mean(np.asarray([point[0] for index,point in enumerate(pole_2_coordinates)])),np.mean(np.asarray([point[1] for index,point in enumerate(pole_2_coordinates)])),np.mean(np.asarray([point[2] for index,point in enumerate(pole_2_coordinates)]))]
# Calculate these points in sphercial coordinates
r_1, theta_1, phi_1 = cart2sph(centroid_1[0],centroid_1[1],centroid_1[2])
r_2, theta_2, phi_2 = cart2sph(centroid_2[0],centroid_2[1],centroid_2[2])
# Now, we re-define our centroids based on these coordinates so that we are sure that we end up on the sphere
centroid_1 = [np.cos(phi_1)*np.cos(theta_1), np.cos(phi_1)*np.sin(theta_1), np.sin(phi_1)]
centroid_2 = [np.cos(phi_2)*np.cos(theta_2), np.cos(phi_2)*np.sin(theta_2), np.sin(phi_2)]
# Calculate the minimal distance between these points
dist_temp_north = m.acos(centroid_1[0]*centroid_2[0]+centroid_1[1]*centroid_2[1]+centroid_1[2]*centroid_2[2])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Calculate the angle between the poles and the holes
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Define the location of the hole in the mesh
hole = [0, 0, -1]
# Define two vectors between the hole and the centroids
#v1 = [centroid_1[0]-hole[0], centroid_1[1]-hole[1], centroid_1[2]-hole[2]]
#v2 = [centroid_2[0]-hole[0], centroid_2[1]-hole[1], centroid_2[2]-hole[2]]
v1 = hole
v2 = centroid_1
# Normalise both vectors
unit_vector_1 = v1 / np.linalg.norm(v1)
unit_vector_2 = v2 / np.linalg.norm(v2)
# Calculate the dot product between these vectors
dot_product = np.dot(unit_vector_1, unit_vector_2)
# Calculate the angle between these two vectors
angle_north = np.arccos(dot_product)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Mesh at (1,0,0)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Read the msh file
conc_profile = meshio.read("../Output/equator/iteration_0/u000101.vtu")
# Extract the concentration profile
u = np.asarray(list(conc_profile.point_data.values())[0])
## Extract the spatial coordinates
spatial_coordinates = conc_profile.points
# Define a threshold concentration
threshold_concentration = 0.8*max(u)
# Maximum conc
max_conc_equator = max(u)
# Save the spatial coordinates for which the concentration profile is above our threshold value
pole_coordinates = np.asarray([spatial_coordinate for index,spatial_coordinate in enumerate(spatial_coordinates) if u[index]>=threshold_concentration])
# CLUSTERING: conduct the density based scan using DBSCAN
db = DBSCAN(eps=epsilon, min_samples=1).fit(pole_coordinates)
cluster_labels = db.labels_
# Get the number of poles
num_poles_equator = len(set(cluster_labels))
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Calculate the relative pole area by loading the old mesh and integrating
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Load the old mesh
mesh_old = Mesh("../Output/equator/iteration_0/final_timestep_mesh.xml")
gdim = mesh_old.geometry().dim()
# Define a function space on the old mesh
H_old = FunctionSpace(mesh_old, "P", 1)
# Load the old initial conditions on the old function space
u_old = Function(H_old, "../Output/equator/iteration_0/final_timestep_u.xml")
v_old = Function(H_old, "../Output/equator/iteration_0/final_timestep_v.xml")
# Find out the coordinates in the mesh
# Create a dofmap
dofmap = H_old.dofmap()
sphere_dofs = dofmap.dofs()
# Find all coordinates for the nodes in our mesh
coordinates = H_old.tabulate_dof_coordinates()
# Calculate the coordinates in the pole
pole_1_coordinates = []
pole_1_dofs = []
pole_2_coordinates = []
pole_2_dofs = []
non_pole_dofs = []
for index,coordinate in enumerate(coordinates):
    for sub_index,temp_coordinate in enumerate(pole_coordinates):
        if (coordinate[0]==temp_coordinate[0]) and ((coordinate[1]==temp_coordinate[1])) and ((coordinate[2]==temp_coordinate[2])):
            if cluster_labels[sub_index] == 0:
                pole_1_coordinates.append(coordinate)
                pole_1_dofs.append(sphere_dofs[index])
            elif cluster_labels[sub_index] == 1:
                pole_2_coordinates.append(coordinate)
                pole_2_dofs.append(sphere_dofs[index])                
            else:
                non_pole_dofs.append(sphere_dofs[index])                
# Create a function as well 
pole_1_function = Function(H_old)
pole_2_function = Function(H_old)
# Set this pole function to 1 in the pole and 0 outside
pole_1_function.vector()[non_pole_dofs] = 0
pole_1_function.vector()[pole_1_dofs] = 1
pole_2_function.vector()[non_pole_dofs] = 0
pole_2_function.vector()[pole_2_dofs] = 1
# Define an integration measure
dx = Measure("dx", domain=mesh_old)
# Integrate over domain
pole_1_area_integration_equator = 100*assemble(pole_1_function*dx)/assemble(Constant(1.0)*dx)
pole_2_area_integration_equator = 100*assemble(pole_2_function*dx)/assemble(Constant(1.0)*dx)    
relative_pole_area_equator = pole_1_area_integration_equator + pole_2_area_integration_equator
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Calculate the relative pole area by loading the old mesh and integrating
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Calculate the centroids as just the centres of mass of the poles
centroid_1 = [np.mean(np.asarray([point[0] for index,point in enumerate(pole_1_coordinates)])),np.mean(np.asarray([point[1] for index,point in enumerate(pole_1_coordinates)])),np.mean(np.asarray([point[2] for index,point in enumerate(pole_1_coordinates)]))]
centroid_2 = [np.mean(np.asarray([point[0] for index,point in enumerate(pole_2_coordinates)])),np.mean(np.asarray([point[1] for index,point in enumerate(pole_2_coordinates)])),np.mean(np.asarray([point[2] for index,point in enumerate(pole_2_coordinates)]))]
# Calculate these points in sphercial coordinates
r_1, theta_1, phi_1 = cart2sph(centroid_1[0],centroid_1[1],centroid_1[2])
r_2, theta_2, phi_2 = cart2sph(centroid_2[0],centroid_2[1],centroid_2[2])
# Now, we re-define our centroids based on these coordinates so that we are sure that we end up on the sphere
centroid_1 = [np.cos(phi_1)*np.cos(theta_1), np.cos(phi_1)*np.sin(theta_1), np.sin(phi_1)]
centroid_2 = [np.cos(phi_2)*np.cos(theta_2), np.cos(phi_2)*np.sin(theta_2), np.sin(phi_2)]
# Calculate the minimal distance between these points
dist_temp_equator = m.acos(centroid_1[0]*centroid_2[0]+centroid_1[1]*centroid_2[1]+centroid_1[2]*centroid_2[2])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Calculate the angle between the poles and the holes
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Define the location of the hole in the mesh
hole = [1, 0, 0]
# Define two vectors between the hole and the centroids
#v1 = [centroid_1[0]-hole[0], centroid_1[1]-hole[1], centroid_1[2]-hole[2]]
#v2 = [centroid_2[0]-hole[0], centroid_2[1]-hole[1], centroid_2[2]-hole[2]]
v1 = hole
v2 = centroid_1
# Normalise both vectors
unit_vector_1 = v1 / np.linalg.norm(v1)
unit_vector_2 = v2 / np.linalg.norm(v2)
# Calculate the dot product between these vectors
dot_product = np.dot(unit_vector_1, unit_vector_2)
# Calculate the angle between these two vectors
angle_equator = np.arccos(dot_product)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Print comparison
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
print("-------------------------------------------------------------------------------------")
print("Property\t\t\t(0,0,-1)\t(1,0,0)")
print("-------------------------------------------------------------------------------------")
print("Max.conc\t\t\t%0.3f\t\t%0.3f"%(max_conc_north,max_conc_equator))
print("Num. poles\t\t\t%d\t\t%d"%(num_poles_north,num_poles_equator))
print("Pole 1\t\t\t\t%0.3f\t\t%0.3f"%(pole_1_area_integration_north,pole_1_area_integration_equator))
print("Pole 2\t\t\t\t%0.3f\t\t%0.3f"%(pole_2_area_integration_north,pole_2_area_integration_equator))
print("Tot area\t\t\t%0.3f\t\t%0.3f"%(relative_pole_area_north,relative_pole_area_equator))
print("Distance\t\t\t%0.3f\t\t%0.3f"%(dist_temp_north,dist_temp_equator))
print("Angle\t\t\t\t%0.3f\t\t%0.3f"%(angle_north,angle_equator))
print("-------------------------------------------------------------------------------------")
