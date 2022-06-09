# =================================================================================
# =================================================================================
# Script:"analyse_data"
# Date: 2022-06-08
# Implemented by: Johannes Borgqvist
# Description:
# A temporary script for calculating the pole area using integration and finding
# the distance to the poles.
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import meshio # To extract the important parts of the mesh
import pandas as pd # Pandas for dataframes
import numpy as np # Import numpy as well
import Schnakenberg_properties # Home made
from matplotlib import pyplot as plt # Do the plotting using matplotlib
from mpl_toolkits.mplot3d import Axes3D # <--- This is important for 3d plotting 
from sklearn.cluster import DBSCAN # To calculate the number of poles in the concentration profile
import math as m # Mathematical functions
import matplotlib.font_manager as font_manager # For the font manager
from shapely.geometry import MultiPoint
# Most good things are contained in FEniCS
from fenics import *
import os
from dolfin import *
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
#------------------------------------------------------------------
# Function 1: "read_mesh_Schnakenberg_sphere_with_holes"
# The function takes the number of holes as inputs and
# returns the following four outputs:
# 1. The number of holes in a variable called "num_holes",
# 2. A list of all the radii of the holes called "radii_holes",
# 3. A mesh_function "mf_subdomain",
# 4. A list of integration measures "dx_list" used in the
# variational formulation.
#------------------------------------------------------------------

def read_mesh_Schnakenberg_sphere_with_holes(num_holes,radii_holes):
    # Allocate memory for the mesh and the mesh value collection
    mesh = Mesh()
    mvc_subdomains = MeshValueCollection("size_t", mesh, 2)
    # Define a mesh function taking the subdomains into account
    mf_subdomains = 0
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
    print(mesh_str)
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
def get_centermost_point(cluster):
    centroid = (MultiPoint(cluster).centroid.x, MultiPoint(cluster).centroid.y)
    centermost_point = min(cluster, key=lambda point: great_circle(point, centroid).m)
    return tuple(centermost_point)

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
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# DEFINE THE FOLDERS WE LOOK THROUGH
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
folder_str = "../Output/"
a_str = "a_" + str(round(a,3)).replace(".","p") + "_"
b_str = "b_" + str(round(b,3)).replace(".","p") + "_"
d_str = "d_" + str(round(d,3)).replace(".","p") + "_"
gamma_str = "gamma_" + str(round(gamma,3)).replace(".","p") + "_"
sigma_str = "sigma_" + str(round(sigma,5)).replace(".","p") + "_"
T_str = "T_" + str(round(T,3)).replace(".","p") + "_"
if ICs_around_steady_states:
    IC_str = "ICs_around_steady_states/"
else:
    IC_str = "ICs_at_zero/"
# Define the radius and the hole string
hole_str = "h_0_"
# Radius string
radius_str = ""
# Gather all these substrings into one giant string where we will save the output files
output_folder = folder_str + hole_str + radius_str + a_str + b_str + d_str + gamma_str + sigma_str + T_str + IC_str
mesh_name =  output_folder + "iteration_" + str(repitition_index) + "/u000101"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Read the concentration profile and spatial coordinates
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Read the msh file
conc_profile = meshio.read(mesh_name + ".vtu")
# Extract the concentration profile
u = np.asarray(list(conc_profile.point_data.values())[0])
# Extract the spatial coordinates
spatial_coordinates = conc_profile.points
# Save each of the coordinates as numpy arrays (i.e. we get three arrays named x, y and z)
x = np.asarray([spatial_coordinates[index][0] for index in range(len(spatial_coordinates))])
y = np.asarray([spatial_coordinates[index][1] for index in range(len(spatial_coordinates))])
z = np.asarray([spatial_coordinates[index][2] for index in range(len(spatial_coordinates))])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the concentration profile using matplotlib
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
img = ax.scatter(x, y, z, c=u, cmap='viridis', alpha=1)
cb = plt.colorbar(img, orientation="horizontal")
cb.set_label(label='Concentration profile $u(\mathbf{x},t=50)$',fontname="Times New Roman", size=15)
# changing the fontsize of yticks
plt.xticks(fontname="Times New Roman", size=10)
plt.yticks(fontname="Times New Roman", size=10)
ax.set_xlabel('$x$',fontname="Times New Roman", size=15)
ax.set_ylabel('$y$',fontname="Times New Roman", size=15)
ax.set_zlabel('$z$',fontname="Times New Roman", size=15)
ax.set_title('Concentration profile at $t=50$ with Chaplain\'s reference parameters',fontname="Times New Roman", size=30,fontweight="bold")

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Lazy way of calculating pole area
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Define a threshold concentration
threshold_concentration = 0.8*max(u)
# Save the spatial coordinates for which the concentration profile is above our threshold value
spatial_coordinates = np.asarray([[x[index], y[index], z[index]] for index in range(len(u)) if u[index]>=threshold_concentration])
# CLUSTERING: conduct the density based scan using DBSCAN
epsilon = 1
db = DBSCAN(eps=epsilon, min_samples=1).fit(spatial_coordinates)
cluster_labels = db.labels_
# Get the number of poles
num_poles = len(set(cluster_labels))
# Get the clusters alright
#clusters = pd.Series([coords[cluster_labels == n] for n in range(num_poles)])
print("The number of poles is:\t\t\t%d,"%(num_poles))
print("Now, we are calculating the pole area the lazy way:")
relative_pole_area = 100*len(spatial_coordinates)/len(u)
print("The relative pole area:\t\t\t%2.2f%%,"%(relative_pole_area))
pole_area_1 = 100*(len([label for label in cluster_labels if label==0])/len(u))
pole_area_2 = 100*(len([label for label in cluster_labels if label==1])/len(u))
print("The relative pole area of pole 1:\t%2.2f%%,"%(pole_area_1))
print("The relative pole area of pole 2:\t%2.2f%%."%(pole_area_2))
# Extract the coordinates belong to pole 1
pole_1_coordinates = np.asarray([spatial_coordinates[index] for index,label in enumerate(cluster_labels) if label==0])
# Covert to spherical coordinates
pole_1_spherical = np.asarray([list(cart2sph(coordinate[0],coordinate[1],coordinate[2])) for coordinate in pole_1_coordinates])
# Extract the coordinates belong to pole 2
pole_2_coordinates = np.asarray([spatial_coordinates[index] for index,label in enumerate(cluster_labels) if label==1])
# Covert to spherical coordinates
pole_2_spherical = np.asarray([list(cart2sph(coordinate[0],coordinate[1],coordinate[2])) for coordinate in pole_2_coordinates])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the poles in spherical coordinates
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
fig = plt.figure()
ax = fig.add_subplot(111)
# Pole 1
plt.scatter(pole_1_spherical[0][2], pole_1_spherical[0][1], marker='o',color=(77/256,0/256,75/256),label='Pole 1')
for index in range(1,len(pole_1_spherical)):
    plt.scatter(pole_1_spherical[index][2], pole_1_spherical[index][1], marker='o',color=(77/256,0/256,75/256))
# Pole 2
plt.scatter(pole_2_spherical[0][2], pole_2_spherical[0][1], marker='o',color=(102/256,37/256,6/256),label='Pole 2')
for index in range(1,len(pole_2_spherical)):
    plt.scatter(pole_2_spherical[index][2], pole_2_spherical[index][1], marker='o',color=(102/256,37/256,6/256))
font = font_manager.FontProperties(family="Times New Roman",
                                   weight="bold",
                                   style="normal", size=16)
plt.legend(prop=font)    
# changing the fontsize of yticks
plt.xticks(fontname="Times New Roman", size=10)
plt.yticks(fontname="Times New Roman", size=10)
plt.xlim(-np.pi, np.pi)
plt.ylim(-np.pi/2, np.pi/2)
ax.set_xlabel('Angle, $\\theta$',fontname="Times New Roman", size=15)
ax.set_ylabel('Angle, $\\phi$',fontname="Times New Roman", size=15)
ax.set_title('Poles projected on the $(\\theta,\\phi)$-plane',fontname="Times New Roman", size=30,fontweight="bold")
plt.show()
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Try to calculate the pole area the non-lazy way
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Load the old mesh
mesh_old = Mesh(output_folder + "final_timestep_mesh.xml")
gdim = mesh_old.geometry().dim()
# Define a function space on the old mesh
H_old = FunctionSpace(mesh_old, "P", 1)
# Load the old initial conditions on the old function space
u_old = Function(H_old, output_folder + "final_timestep_u.xml")
v_old = Function(H_old, output_folder + "final_timestep_v.xml")
# Convert these two into vectors
u_vec = u_old.vector()
# Calculate the maximum value
u_max_vec = u_vec.max()
print("FEniCS maximum,\tu_max\t=\t%0.3f\n"%(u_max_vec))
print("Data file maximum,\tu_max\t=\t%0.3f\n"%(max(u)))
print(u_vec[0])
# Find out the coordinates in the mesh
# Create a dofmap
dofmap = H_old.dofmap()
sphere_dofs = dofmap.dofs()



coordinates = H_old.tabulate_dof_coordinates()
# Calculate the coordinates in the pole
pole_1_coordinates = []
pole_1_dofs = []
pole_2_coordinates = []
pole_2_dofs = []
non_pole_dofs = []
for index,coordinate in enumerate(coordinates):
    for sub_index,temp_coordinate in enumerate(spatial_coordinates):
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
pole_1_area_integration = assemble(pole_1_function*dx)/assemble(Constant(1.0)*dx)
pole_2_area_integration = assemble(pole_2_function*dx)/assemble(Constant(1.0)*dx)
print("The relative pole area of pole 1:\t\t\t%0.5f%%,"%(100*pole_1_area_integration))
print("The relative pole area of pole 2:\t\t\t%0.5f%%,"%(100*pole_2_area_integration))
print("The total relative pole area:\t\t\t\t%0.5f%%,"%(100*(pole_1_area_integration+pole_2_area_integration)))
# Get all the coordinates
# Get coordinates as len(dofs) x gdim array

#print(dofmap.dofs())

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Try to calculate the center point
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
poles = []
for index,coordinate in enumerate(pole_1_coordinates):
    poles.append(coordinate)
for index,coordinate in enumerate(pole_2_coordinates):
    poles.append(coordinate)    
cluster_1 = pd.Series(pole_1_coordinates)
cluster_2 = pd.Series(pole_2_coordinates)
clusters = pd.Series(poles)
#cluster_1 = [label for label in cluster_labels if label==0]
#cluster_2 = [label for label in cluster_labels if label==1]
print(MultiPoint(cluster_1).centroid)
print(MultiPoint(cluster_1).centroid.x)
print(MultiPoint(cluster_1).centroid.y)
#centroid_1 = (MultiPoint(cluster_1).centroid.x, MultiPoint(cluster_1).centroid.y,MultiPoint(cluster_1).centroid.z)
