# =================================================================================
# =================================================================================
# Script:"launch_simulations_Schnakenberg_sphere_with_holes"
# Date: 2022-03-30
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
# The wavenumber k^2
n = 2
k_squared = n*(n+1)
# Calculate the steady states and the critical parameters
u_0, v_0, d_c, gamma_c = Schnakenberg_properties.calculate_steady_states_and_critical_parameters_Schnakenberg(a,b,k_squared)
# Save the steady states in a list
steady_states = [u_0,v_0]
# Prompt to the user
print("---------------------------------------------------------------------------------------------------------\n")
print("\tREPLICATING CHAPLAINS SIMULATIONS")
print("---------------------------------------------------------------------------------------------------------\n")
# Print the results
print("\n\t\tThe steady states:\t\t\t(u_0,v_0)\t=\t(%0.4f,%0.4f)"%(u_0,v_0))
print("\t\tThe critical parameters:\t\t(d_c,gamma_c)\t=\t(%0.4f,%0.4f)"%(d_c,gamma_c))
# Set the value of the relative diffusion
#d = d_c + 1.5
#d = d_c + 10
d = 18
# Set the value of the reaction strength to its critical value
gamma = gamma_c
# Compute minimal critical hole radius for pattern disturbance
n_largest = 1
eps_tuple, n_tuple, m_tuple = Schnakenberg_properties.compute_minimal_holeradius_for_pattern_disturbance(a,b,d,gamma,n,n_largest)
hole_cylinder_radius = np.sin(eps_tuple[1])
# Print the results
print("\n\t\tTheoretical parameter values for pattern disturbances") 
print("\t\tThe critical radii:\t\t(geodesic, cylindric)\t=\t(%0.4f,%0.4f)"%(eps_tuple[1], hole_cylinder_radius))
print("\t\tThe critical spectral parameters:\t\t(n,m)\t=\t(%0.4f,%0.4f)"%(n_tuple[1], m_tuple[1]))
# Collect all parameters in a list
parameters = [a, b, d, gamma]
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
hole_radius_array = np.arange(0,0.75,0.05)
# Loop over the hole_radii and add the experiments
for hole_radius in hole_radius_array:
    # Special case for the mesh with no hole
    if hole_radius == 0:
        experimental_design.append((0,parameters,numerical_parameters,[],True))
    else:
        experimental_design.append((1,parameters,numerical_parameters,[hole_radius],True))
# We repeat the experiments a number of time due to the stochasticity in the intial conditions
number_of_repititions = 10        
# Loop over the experiments in the experimental design and run them all (with the appropriate number of repititions)
for experiment in experimental_design:
    # Prompt to the user
    print("---------------------------------------------------------------------------------------------------------\n")
    print("\tNUM_HOLES\t=\t%d,\tRADII\t=\t%s"%(int(experiment[0]),str(experiment[3])))
    print("---------------------------------------------------------------------------------------------------------\n")    
    # Solve the FEM system with the given parameters
    FEM_toolbox.FEMFD_simulation_Schnakenberg_sphere_with_holes(experiment[0],experiment[1],steady_states,experiment[2],experiment[3],experiment[4],number_of_repititions)    

