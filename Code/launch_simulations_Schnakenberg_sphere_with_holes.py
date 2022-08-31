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
print("\tRUNNING SCHNAKENBERG'S MODEL ON A SPHERICAL MESH WITH A HOLE IN IT")
print("---------------------------------------------------------------------------------------------------------\n")
# Set the value of the relative diffusion
# BEFORE d-interval-calibration/maxing d (for n > 2)
#d = 30 # For n=1
#d = 18 # For n=2 
#d = 22 # For n=3
# AFTER d-interval-calibration/maxing d (for n > 2)
#d = 30 # For n=1
#d = 18 # For n=2 
#d = 21.75 # For n=3
#d = 19.75 # For n=4
#d = 18.75 # For n=5
#d = 18.23 # For n=6
# Calibrated values of d based on simulations (for n=1,2,3)
#d_vec = [20, 18, 19] # n=1,2,3
#d_vec = [18] # n=4
d_vec = [18] # New attempt with d=18 for n=3. 
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
#hole_radius_array = np.arange(0,0.75,0.05) # The full experimental design
hole_radius_array = np.array([0]) # Calibrate d-value on the mesh without hole
# Define the eigenvalues we want to consider
#n_vec = [1, 3]
#n_vec = [4]
n_vec = [5]
# Loop over the eigenvalues
for n_index,n in enumerate(n_vec):
    k_squared = n*(n+1)
    # Calculate the steady states and the critical parameters
    u_0, v_0, d_c, gamma_c = Schnakenberg_properties.calculate_steady_states_and_critical_parameters_Schnakenberg(a,b,k_squared)
    # Save the steady states in a list
    steady_states = [u_0,v_0]
    # Set the value of the relative diffusion to its given value
    d = d_vec[n_index]
    # Set the value of the reaction strength to its critical value
    gamma = gamma_c
    # Collect all parameters in a list
    parameters = [a, b, d, gamma, n]
    # Check isolated spectral modes for choices of a, b, d, and gamma.
    Turing_conditions,L,M = Schnakenberg_properties.check_Turing_conditions_Scnakenberg(a,b,d)
    print(gamma*L, gamma*M)
    ninterval = [1, 10]
    Schnakenberg_properties.compute_isolated_spectral_modes(ninterval, gamma, L, M)
    # Print the results
    print("\n\t\tThe steady states:\t\t\t(u_0,v_0)\t=\t(%0.4f,%0.4f)"%(u_0,v_0))
    print("\t\tThe critical parameters:\t\t(d_c,gamma_c)\t=\t(%0.4f,%0.4f)"%(d_c,gamma_c))    
    # Loop over the hole_radii and add the experiments
    for hole_index,hole_radius in enumerate(hole_radius_array):
        # Special case for the mesh with no hole
        if hole_radius == 0:
            experimental_design.append((0,parameters,steady_states,numerical_parameters,[],True,False))
        else:
            experimental_design.append((1,parameters,steady_states,numerical_parameters,[hole_radius],True,False))
# Here, we define the start repititions and the number of repititions
#number_of_repititions = 20 # For the full experimental design
number_of_repititions = 1 # For a single repitition when calibrating the d-value
start_repitition = 0 # This value we can tweek if we want to add extra simulations afterwards
# Loop over the experiments in the experimental design and run them all (with the appropriate number of repititions)
for experiment in experimental_design:
    # Prompt to the user
    print("---------------------------------------------------------------------------------------------------------\n")
    print("\tNUM_HOLES\t=\t%d,\tRADII\t=\t%s"%(int(experiment[0]),str(experiment[4])))
    print("---------------------------------------------------------------------------------------------------------\n")    
    # Solve the FEM system with the given parameters
    FEM_toolbox.FEMFD_simulation_Schnakenberg_sphere_with_holes(experiment[0],experiment[1],experiment[2],experiment[3],experiment[4],experiment[5],number_of_repititions,experiment[6],start_repitition)

