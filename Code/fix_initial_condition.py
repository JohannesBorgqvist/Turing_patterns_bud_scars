# =================================================================================
# =================================================================================
# Script:"fix_initial_condition"
# Date: 2022-04-29
# Implemented by: Johannes Borgqvist
# Description:
# This script generates an initial condition that is perturbed from the steady-state
# and then it saves this IC so that it can be read from a file.
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
# Set the number of holes
num_holes = 0
# We have no holes, so no radius necessary
radii_holes = []
# Define the perturbation in the initial conditions
sigma = 1e-4
# Define the end time for the simulations
T = 50
# Collect these latter two parameters in a list as well
numerical_parameters = [sigma, T]
# We set the initial conditions around the steady states
ICs_around_steady_states = True
# Save the initial conditio
FEM_toolbox.save_IC(num_holes,steady_states,numerical_parameters,radii_holes,ICs_around_steady_states)
