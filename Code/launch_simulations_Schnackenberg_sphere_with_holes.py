# =================================================================================
# =================================================================================
# Script:"launch_simulations_Schnackenberg_sphere_with_holes"
# Date: 2021-09-02
# Implemented by: Johannes Borgqvist
# Description:
# This is the script which launches the FEM simulations, and
# it uses all functions that are stored in the script
# "toolbox_FEM_simulations_Schnackenberg_sphere_with_holes.py".
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import toolbox_FEM_simulations_Schnackenberg_sphere_with_holes as FEM_toolbox  # Home-made
# =================================================================================
# =================================================================================
# The experiments
# =================================================================================
# =================================================================================
#----------------------------------------------------------------------------------
# Experiment 1 in Chaplain:
# No cell growth, no holes and no local activation
#----------------------------------------------------------------------------------
# The parameters in the Schnackenberg model
a = 0.2
b = 1
# The wavenumber
k_squared = 6
# Calculate the steady states and the critical parameters
u_0, v_0, d_c, gamma_c = FEM_toolbox.calculate_steady_states_and_critical_parameters_Schnackenberg(a,b,k_squared)
# Save the steady states in a list
steady_states = [u_0,v_0]
# Print the results
print("\n\t\tThe steady states:\t\t\t(u_0,v_0)\t=\t(%0.4f,%0.4f)"%(u_0,v_0))
print("\n\t\tThe critical parameters:\t\t(d_c,gamma_c)\t=\t(%0.4f,%0.4f)"%(d_c,gamma_c))
# Set the value of the relative diffusion
d = 15
# Set the value of the reaction strength to its critical value
gamma = gamma_c
# Collect all parameters in a list
parameters = [a, b, d, gamma]
# For this experiment we have no hole in the mesh
num_holes = 0
# Define the perturbation in the initial conditions
sigma = 0.05
# Define the end time for the simulations
T = 50
# Collect these latter two parameters in a list as well
numerical_parameters = [sigma, T]
# Since, we have no hole we have no adjacent region and therefore it is not meaningful to talk about local activation in this region. So we set the activation parameters to 1 meaning that we have no extra activation in the region adjacent to the hole
activation_parameters = [1, 1]
# We do not have any cell growth here so we set this to false
cell_growth = False
# Solve the FEM system with the given parameters
FEM_toolbox.FEMFD_simulation_Schnackenberg_sphere_with_holes(num_holes,parameters,steady_states,numerical_parameters,activation_parameters,cell_growth)
#FEM_toolbox.FEM_FD_simulation_pure_diffusion_sphere_with_holes()

