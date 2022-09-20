# =================================================================================
# =================================================================================
# Script:"data_analysis_of_spatial_patterns"
# Date: 2022-06-08
# Implemented by: Johannes Borgqvist
# Description:
# This script analyses the spatial patterns by reading in the concentration profile
# from the vtu-files generated by the FEM simulations. The main scripts here are
# numpy for handling arrays etc, meshio for reading the vtu file as well as
# extracting the data and lastly 
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
n = 1
k_squared = n*(n+1)
# Calculate the steady states and the critical parameters
u_0, v_0, d_c, gamma_c = Schnakenberg_properties.calculate_steady_states_and_critical_parameters_Schnakenberg(a,b,k_squared)
# Save the steady states in a list
steady_states = [u_0,v_0]
# Set the value of the relative diffusion
d = 20.0 # n=1
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
# Let's start with the zeroth repitition
repitition_index = 0
# Define the meshes we want to loop over
hole_radius_array = np.arange(0,0.75,0.05)
#hole_radius_array = np.arange(0, 0.35, 0.05)
# Define a parameter epsilon for the clustering
epsilon = 1
# Define the location of the hole in the mesh
hole = [0, 0, -1]
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Read the concentration profile and spatial coordinates
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Allocate memory for the three metric we will plot
individual_pole_area = []
rel_pol_area = []
num_poles_vec = []
max_conc = []
min_dist = []
# Read ALL the data by looping over all hole radii
for hole_index in range(len(hole_radius_array)):
    if hole_index == 0:
        hole_str = "h_0_"
        radii_holes = []
    else:
        hole_str = "h_1_"
        radii_holes = [hole_radius_array[hole_index]]
    # Define an empty hole radius string    
    radius_str = ""
    # Loop over all raddi and read them one by one
    for radius in radii_holes:
        radius_str += "r_" + str(round(radius,3)).replace(".","p") + "_"
    # Allocate the memory for all the outputs
    num_poles = []
    u_max = []
    individual_pole_area_integration = []
    relative_pole_area = []
    dist_temp = []
    # Loop over all repititions
    for repitition_index in range(20):
        # Gather all these substrings into one giant string where we will save the output files
        mesh_name = folder_str + hole_str + radius_str + a_str + b_str + d_str + gamma_str + sigma_str + T_str + IC_str + "iteration_" + str(repitition_index) + "/u000001"
        # Read the msh file
        conc_profile = meshio.read(mesh_name + ".vtu")
        # Extract the concentration profile
        u = np.asarray(list(conc_profile.point_data.values())[0])
        # Extract the spatial coordinates
        spatial_coordinates = conc_profile.points
        # Define a threshold concentration
        threshold_concentration = 0.95*max(u)
        # Save the spatial coordinates for which the concentration profile is above our threshold value
        pole_coordinates = np.asarray([spatial_coordinate for index,spatial_coordinate in enumerate(spatial_coordinates) if u[index]>=threshold_concentration])
        # CLUSTERING: conduct the density based scan using DBSCAN
        #db = DBSCAN(eps=epsilon, min_samples=1).fit(pole_coordinates)
        db = DBSCAN(eps=epsilon, min_samples=1).fit(pole_coordinates)
        cluster_labels = db.labels_
        # Get the number of poles
        num_poles.append(len(set(cluster_labels)))
        u_max.append(max(u))
        #----------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------
        # Calculate the relative pole area by loading the old mesh and integrating
        #----------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------
        # Load the old mesh
        mesh_old = Mesh(folder_str + hole_str + radius_str + a_str + b_str + d_str + gamma_str + sigma_str + T_str + IC_str + "iteration_" + str(repitition_index) +"/final_timestep_mesh.xml")
        gdim = mesh_old.geometry().dim()
        # Define a function space on the old mesh
        H_old = FunctionSpace(mesh_old, "P", 1)
        # Load the old initial conditions on the old function space
        u_old = Function(H_old, folder_str + hole_str + radius_str + a_str + b_str + d_str + gamma_str + sigma_str + T_str + IC_str + "iteration_" + str(repitition_index) +"/final_timestep_u.xml")
        v_old = Function(H_old, folder_str + hole_str + radius_str + a_str + b_str + d_str + gamma_str + sigma_str + T_str + IC_str  + "iteration_" + str(repitition_index) + "/final_timestep_v.xml")
        # Find out the coordinates in the mesh
        # Create a dofmap
        dofmap = H_old.dofmap()
        sphere_dofs = dofmap.dofs()
        # Find all coordinates for the nodes in our mesh
        coordinates = H_old.tabulate_dof_coordinates()
        # Calculate the coordinates in the pole
        # Allocate memory
        specific_pole_coordinates = []
        specific_pole_dofs = []
        # For each pole we add a new list
        for label in set(cluster_labels):
            specific_pole_coordinates.append([])
            specific_pole_dofs.append([])
            individual_pole_area_integration.append([])            
        # All the dofs that are not in the pole
        non_pole_dofs = []
        # Loop over coordinates and save the ones in the respective poles
        for index,coordinate in enumerate(coordinates): # Loop over all coordinates
            # We assume that the coordinate at hand is not in the pole
            non_pole_indicator = True
            # Loop over all pole coordinates
            for sub_index,temp_coordinate in enumerate(pole_coordinates): 
                # We have a coordinate in the pole
                if (coordinate[0]==temp_coordinate[0]) and ((coordinate[1]==temp_coordinate[1])) and ((coordinate[2]==temp_coordinate[2])):
                    # Save our coordinate and the dof of the specific pole
                    specific_pole_coordinates[cluster_labels[sub_index]].append(coordinate)
                    specific_pole_dofs[cluster_labels[sub_index]].append(sphere_dofs[index])
                    # Set the non_pole_indicator to false
                    non_pole_indicator = False
                    # Stop looping
                    break
            # The coordinate is not in the pole, hey?
            if non_pole_indicator:
                # Append the dof
                non_pole_dofs.append(sphere_dofs[index])
        # Define an integration measure
        dx = Measure("dx", domain=mesh_old)
        # Do a temporary sum of the total area
        total_area = 0
        # Loop over all poles and calculate their area
        for index,label in enumerate(set(cluster_labels)):
            # Create a function
            exec("pole_%d_function=Function(H_old)"%(index+1))
            # Set this pole function to 1 in the pole and 0 outside
            exec("pole_%d_function.vector()[non_pole_dofs] = 0"%(index+1))
            exec("pole_%d_function.vector()[specific_pole_dofs[%d]] = 1"%(index+1,label))
            # Calculate the individual area of the pole by integrating
            exec("individual_pole_area_integration[%d].append(100*assemble(pole_%d_function*dx)/assemble(Constant(1.0)*dx))"%(label,index+1))
            # Add this to the sum as well
            exec("total_area += 100*assemble(pole_%d_function*dx)/assemble(Constant(1.0)*dx)"%(index+1))
        # Append the toal area
        relative_pole_area.append(total_area)

        #----------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------
        # Calculate the minimial distance between our points
        #----------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------
        # Calculate the distances between the hole in the mesh and the centre points of the hole
        distances = []
        # Loop over all poles and calculate their centroid being the centre of mass
        for index,label in enumerate(set(cluster_labels)):
            # Calculate the centroids as just the centres of mass of the poles
            centroid = [np.mean(np.asarray([point[0] for index,point in enumerate(specific_pole_coordinates[label])])),np.mean(np.asarray([point[1] for index,point in enumerate(specific_pole_coordinates[label])])),np.mean(np.asarray([point[2] for index,point in enumerate(specific_pole_coordinates[label])]))]
            # Calculate these points in spherical coordinates
            r, theta, phi = cart2sph(centroid[0],centroid[1],centroid[2])
            # Now, we re-define our centroids based on these coordinates so that we are sure that we end up on the sphere
            centroid = [np.cos(phi)*np.cos(theta), np.cos(phi)*np.sin(theta), np.sin(phi)]
            # Now, calculate the distance between the centroid and the hole
            distances.append(m.acos(centroid[0]*hole[0]+centroid[1]*hole[1]+centroid[2]*hole[2]))
        # Lastly, take the minimum of all distances, hey?
        dist_temp.append(np.amin(np.array(distances)))
    #----------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------
    # Append all the values we return in the end
    #----------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------
    # Append all values to our output vectors
    rel_pol_area.append(relative_pole_area)
    individual_pole_area.append(individual_pole_area_integration)
    num_poles_vec.append(num_poles)
    max_conc.append(u_max)
    min_dist.append(dist_temp)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot our metrics as a function of the hole area
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Set all parameters to tex
plt.rcParams['text.usetex'] = True
# Figures
fig, axes = plt.subplots(2,2,figsize=(30,10))
plt.rc('axes', labelsize=35)    # fontsize of the x and y label
plt.rc('legend', fontsize=15)    # legend fontsize
plt.rc('xtick', labelsize=30)    # fontsize of the tick labels
plt.rc('ytick', labelsize=30)    # fontsize of the tick labels
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# Subplot 1 of 4
axes[0][0].plot(hole_radius_array,np.asarray([np.percentile(rel_pol_area[index],50) for index in range(len(rel_pol_area))]),'-',color=(0/256,68/256,27/256),label="Total pole area in \\%")
axes[0][0].plot(hole_radius_array,np.asarray([np.percentile(rel_pol_area[index],95) for index in range(len(rel_pol_area))]),'-',color=(0/256,68/256,27/256))
axes[0][0].plot(hole_radius_array,np.asarray([np.percentile(rel_pol_area[index],5) for index in range(len(rel_pol_area))]),'-',color=(0/256,68/256,27/256))
axes[0][0].fill_between(hole_radius_array,np.asarray([np.percentile(rel_pol_area[index],5) for index in range(len(rel_pol_area))]),np.asarray([np.percentile(rel_pol_area[index],95) for index in range(len(rel_pol_area))]), facecolor=(0/256,68/256,27/256),alpha=0.5,interpolate=True)
axes[0][0].legend(loc='best')
axes[0][0].set_ylim([0,100])
axes[0][0].yaxis.set_ticks(np.arange(0,110,10))
axes[0][0].set_xlim([0,hole_radius_array[-1]])
# We change the fontsize of minor ticks label 
axes[0][0].tick_params(axis='both', which='major', labelsize=15)
axes[0][0].tick_params(axis='both', which='minor', labelsize=15)

# Subplot 2 of 4
axes[0][1].plot(hole_radius_array,np.asarray([np.percentile(num_poles_vec[index],50) for index in range(len(num_poles_vec))]),'-',color=(77/256,0/256,75/256),label="Number of poles")
axes[0][1].plot(hole_radius_array,np.asarray([np.percentile(num_poles_vec[index],95) for index in range(len(num_poles_vec))]),'-',color=(77/256,0/256,75/256))
axes[0][1].plot(hole_radius_array,np.asarray([np.percentile(num_poles_vec[index],5) for index in range(len(num_poles_vec))]),'-',color=(77/256,0/256,75/256))
axes[0][1].fill_between(hole_radius_array,np.asarray([np.percentile(num_poles_vec[index],5) for index in range(len(num_poles_vec))]),np.array([np.percentile(num_poles_vec[index],95) for index in range(len(num_poles_vec))]), facecolor=(77/256,0/256,75/256),alpha=0.5,interpolate=True)
axes[0][1].legend(loc='best')
axes[0][1].set_ylim([0,5])
axes[0][1].set_xlim([0,hole_radius_array[-1]])
axes[0][1].tick_params(axis='both', which='major', labelsize=15)
axes[0][1].tick_params(axis='both', which='minor', labelsize=15)
# Subplot 3 of 4
axes[1][0].plot(hole_radius_array,np.asarray([np.percentile(max_conc[index],50) for index in range(len(max_conc))]),'-',color=(129/256,15/256,124/256),label="Max. conc.")
axes[1][0].plot(hole_radius_array,np.asarray([np.percentile(max_conc[index],95) for index in range(len(max_conc))]),'-',color=(129/256,15/256,124/256))
axes[1][0].plot(hole_radius_array,np.asarray([np.percentile(max_conc[index],5) for index in range(len(max_conc))]),'-',color=(129/256,15/256,124/256))
axes[1][0].fill_between(hole_radius_array,np.asarray([np.percentile(max_conc[index],5) for index in range(len(max_conc))]),np.array([np.percentile(max_conc[index],95) for index in range(len(max_conc))]), facecolor=(129/256,15/256,124/256),alpha=0.5,interpolate=True)
axes[1][0].legend(loc='best')
axes[1][0].set_ylim([0,2])
axes[1][0].set_xlim([0,hole_radius_array[-1]])
axes[1][0].tick_params(axis='both', which='major', labelsize=15)
axes[1][0].tick_params(axis='both', which='minor', labelsize=15)
# Subplot 4 of 4
axes[1][1].plot(hole_radius_array,np.asarray([np.percentile(min_dist[index],50) for index in range(len(min_dist))]),'-',color=(64/256,0/256,125/256),label="Min. dist.")
axes[1][1].plot(hole_radius_array,np.asarray([np.percentile(min_dist[index],95) for index in range(len(min_dist))]),'-',color=(64/256,0/256,125/256))
axes[1][1].plot(hole_radius_array,np.asarray([np.percentile(min_dist[index],5) for index in range(len(min_dist))]),'-',color=(64/256,0/256,125/256))
axes[1][1].fill_between(hole_radius_array,np.asarray([np.percentile(min_dist[index],5) for index in range(len(min_dist))]),np.array([np.percentile(min_dist[index],95) for index in range(len(min_dist))]), facecolor=(64/256,0/256,125/256),alpha=0.5,interpolate=True)
axes[1][1].legend(loc='best')
axes[1][1].set_ylim([0,4])
axes[1][1].yaxis.set_ticks(np.arange(0,4.5,0.5))
axes[1][1].set_xlim([0,hole_radius_array[-1]])
axes[1][1].tick_params(axis='both', which='major', labelsize=15)
axes[1][1].tick_params(axis='both', which='minor', labelsize=15)
# displaying the title
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Cylindrical hole radius, $\\varepsilon$")
plt.ylabel("Quantitative metric")
# displaying the title
plt.title("Metrics of $u(\\mathbf{x},t=50)$ as functions of the hole radius $\\varepsilon$",fontsize=30, fontweight='bold')
# Save fig and show it
plt.savefig("../Figures/patterns_are_preserved_growing_hole_radii_n_" + str(n) + ".png")
plt.show()





