# =================================================================================
# =================================================================================
# Script:"launch_simulations_Schnakenberg_sphere_with_holes"
# Date: 2022-06-09
# Implemented by: Johannes Borgqvist and Carl Lundholm
# Description:
# This is the script which launches the FEM simulations, and
# it uses all functions that are stored in the script
# "toolbox_FEM_simulations_Schnakenberg_sphere_with_holes.py".
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import numpy as np
import toolbox_FEM_simulations_Schnakenberg_sphere_with_holes as FEM_toolbox  # Home-made
import Schnakenberg_properties # Home-made
# =================================================================================
# =================================================================================
# The experiments
# =================================================================================
# =================================================================================
# The parameters in the Schnakenberg model
a = 0.2
b = 1
# Prompt to the user
print("---------------------------------------------------------------------------------------------------------\n")
print("\tREPLICATING CHAPLAINS SIMULATIONS")
print("---------------------------------------------------------------------------------------------------------\n")
# Set the value of the relative diffusion
d = 18
# We have no holes, so no radius necessary
radii_holes = []
# Define the perturbation in the initial conditions
sigma = 1e-4
# Define the end time for the simulations
T = 50
# Collect these latter two parameters in a list as well
numerical_parameters = [sigma, T]
# Looping over the varius radii and run all simulations there!
# Define the experimental design of holes with increasing radii
experimental_design = []
# Define the meshes we want to loop over
#hole_radius_array = np.arange(0,0.75,0.05)
hole_radius_array = np.asarray([0.2, 0.2])
# Define the eigenvalues we want to consider
#n_vec = [1, 2, 3, 4]
n_vec = [2]
# Loop over the eigenvalues
for n in n_vec:
    k_squared = n*(n+1)
    # Calculate the steady states and the critical parameters
    u_0, v_0, d_c, gamma_c = Schnakenberg_properties.calculate_steady_states_and_critical_parameters_Schnakenberg(a,b,k_squared)
    # Save the steady states in a list
    steady_states = [u_0,v_0]
    # Set the value of the reaction strength to its critical value
    gamma = gamma_c
    # Collect all parameters in a list
    parameters = [a, b, d, gamma]
    # Print the results
    print("\n\t\tThe steady states:\t\t\t(u_0,v_0)\t=\t(%0.4f,%0.4f)"%(u_0,v_0))
    print("\t\tThe critical parameters:\t\t(d_c,gamma_c)\t=\t(%0.4f,%0.4f)"%(d_c,gamma_c))    
    # Loop over the hole_radii and add the experiments
    for hole_index,hole_radius in enumerate(hole_radius_array):
        # Special case for the mesh with no hole
        if hole_radius == 0:
            experimental_design.append((0,parameters,steady_states,numerical_parameters,[],True,False))
        else:
            #experimental_design.append((1,parameters,steady_states,numerical_parameters,[hole_radius],True,False))
            if hole_index == 0:
                mesh_name = "../Meshes/s_h_1_r_0p2_north_pole.xdmf"
            else:
                mesh_name = "../Meshes/s_h_1_r_0p2_equator.xdmf"                
            experimental_design.append((1,parameters,steady_states,numerical_parameters,[hole_radius],True,True,mesh_name,hole_index))
# We repeat the experiments a certain number of times due to the stochasticity in the intial conditions
number_of_repititions = 1        
# Loop over the experiments in the experimental design and run them all (with the appropriate number of repititions)
for experiment in experimental_design:
    # Prompt to the user
    print("---------------------------------------------------------------------------------------------------------\n")
    print("\tNUM_HOLES\t=\t%d,\tRADII\t=\t%s"%(int(experiment[0]),str(experiment[3])))
    print("---------------------------------------------------------------------------------------------------------\n")    
    # Solve the FEM system with the given parameters
    #FEM_toolbox.FEMFD_simulation_Schnakenberg_sphere_with_holes(experiment[0],experiment[1],experiment[2],experiment[3],experiment[4],experiment[5],number_of_repititions,experiment[6])
    FEM_toolbox.FEMFD_simulation_Schnakenberg_sphere_with_holes(experiment[0],experiment[1],experiment[2],experiment[3],experiment[4],experiment[5],number_of_repititions,experiment[6],experiment[7],experiment[8])    

