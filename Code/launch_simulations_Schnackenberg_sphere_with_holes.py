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
import toolbox_FEM_simulations_Schnackenberg_sphere_with_holes  # Home-made
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
u_0, v_0, d_c, gamma_c = toolbox_FEM_simulations_Schnackenberg_sphere_with_holes.calculate_steady_states_and_critical_parameters_Schnackenberg(a,b,k_squared)
# Print the results
print("\n\t\tThe steady states:\t\t\t(u_0,v_0)\t=\t(%0.4f,%0.4f)"%(u_0,v_0))
print("\n\t\tThe critical parameters:\t\t(d_c,gamma_c)\t=\t(%0.4f,%0.4f)"%(d_c,gamma_c))
# Define the number of holes
num_holes = 1
# Test the mesh function just for fun
mesh, mvc_subdomains, mf_subdomains, dx_list = toolbox_FEM_simulations_Schnackenberg_sphere_with_holes.read_mesh_Schnackenberg_sphere_with_holes(num_holes)

