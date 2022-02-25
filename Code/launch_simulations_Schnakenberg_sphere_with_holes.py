# =================================================================================
# =================================================================================
# Script:"launch_simulations_Schnakenberg_sphere_with_holes"
# Date: 2022-02-22
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
# =================================================================================
# =================================================================================
# The experiments
# =================================================================================
# =================================================================================
#----------------------------------------------------------------------------------
# Experiment 1 in Chaplain:
# No cell growth, no holes and no local activation
#----------------------------------------------------------------------------------
# The parameters in the Schnakenberg model
a = 0.2
b = 1
# The wavenumber k^2
n = 2
k_squared = n*(n+1)
# Calculate the steady states and the critical parameters
u_0, v_0, d_c, gamma_c = FEM_toolbox.calculate_steady_states_and_critical_parameters_Schnakenberg(a,b,k_squared)
# Save the steady states in a list
steady_states = [u_0,v_0]
# Print the results
print("\n\t\tThe steady states:\t\t\t(u_0,v_0)\t=\t(%0.4f,%0.4f)"%(u_0,v_0))
print("\t\tThe critical parameters:\t\t(d_c,gamma_c)\t=\t(%0.4f,%0.4f)"%(d_c,gamma_c))
# Set the value of the relative diffusion
#d=5
#d=10
#d = 100000000000000
d = 18
# Set the value of the reaction strength to its critical value
gamma = gamma_c
# Compute minimal critical hole radius for pattern disturbance
#n_largest = 731 # Largest eigenvalue parameter to test for.
#n_largest = 6
n_largest = 6
eps_min, n_min, m_min = FEM_toolbox.compute_minimal_holeradius_for_pattern_disturbance(a,b,d,gamma,n,n_largest)
hole_cylinder_radius = np.sin(eps_min)
# Print the results
print("\n\t\tTheoretical parameter values for pattern disturbances") 
print("\t\tThe critical radii:\t\t(geodesic, cylindric)\t=\t(%0.4f,%0.4f)"%(eps_min, hole_cylinder_radius))
print("\t\tThe critical spectral parameters:\t\t(n,m)\t=\t(%0.4f,%0.4f)"%(n_min, m_min))
# Collect all parameters in a list
parameters = [a, b, d, gamma]
# For this experiment we have no hole in the mesh
num_holes = 0
#num_holes = 1
#num_holes = 2
#num_holes = 5
# Define the perturbation in the initial conditions
#sigma = 0.01
sigma = 1e-4
# Define the end time for the simulations
T = 50
#T = 0.05
#T = 1.5
#T = 0.0
#T = 2.0
# Collect these latter two parameters in a list as well
numerical_parameters = [sigma, T]
# Since, we have no hole we have no adjacent region and therefore it is not meaningful to talk about local activation in this region. So we set the activation parameters to 1 meaning that we have no extra activation in the region adjacent to the hole
activation_parameters = [1.0, 1.0]
# We do not have any cell growth here so we set this to false
cell_growth = False
# Solve the FEM system with the given parameters
FEM_toolbox.FEMFD_simulation_Schnakenberg_sphere_with_holes(num_holes,parameters,steady_states,numerical_parameters,activation_parameters,cell_growth)
#FEM_toolbox.FEM_FD_simulation_pure_diffusion_sphere_with_holes()

