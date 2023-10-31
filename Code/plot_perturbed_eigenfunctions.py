# =================================================================================
# =================================================================================
# Script:"plot_perturbed_eigenfunctions"
# Date: 2022-04-29
# Implemented by: Johannes Borgqvist
# Description:
# This script reads the eigenfunctions from the final time step in the launched
# simulations
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import Schnakenberg_properties # Home made
import numpy as np # For numerical calculations
from matplotlib import pyplot as plt # For plotting
import pandas as pd # Import pandas for reading the csv files
from pandas import read_csv # Pandas to read the data
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
# Define the meshes we want to loop over
hole_radius_array = np.arange(0,0.75,0.05)
# Allocate a list of all the basis functions
basis_functions = []
# Let's add 21 lists for each basis function
# corresponding to n=1,2,3,4,5
for index in range(21):
    basis_functions.append([])   
# Now for each hole radius add an empty list as well with all iterations
for index in range(21):
    for sub_index in range(len(hole_radius_array)):
        basis_functions[index].append([])
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# DEFINE THE FOLDERS WE LOOK THROUGH
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
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Allocate memory for the legend_strings
legend_strings = []
# Add a help variable
help_variable = 0
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
    # Now, the hole and the radii string is finished, so let's loop through all of the iterations
    for repitition_index in range(20):
        # Gather all these substrings into one giant string where we will save the output files
        output_folder_str = folder_str + hole_str + radius_str + a_str + b_str + d_str + gamma_str + sigma_str + T_str + IC_str + "iteration_" + str(repitition_index) + "/"
        # Read the csv file
        dataframe = read_csv(output_folder_str + "spectral_coefficients.csv", header=None)
        # Save the legends one time
        if help_variable == 0:
            # Save the legend strings
            legend_strings = list(dataframe.values[:,1])
            # Remove the very first element hey?
            del legend_strings[0]
            # Rename the legends
            for index,legend in enumerate(legend_strings):
                legend_strings[index] = legend_strings[index].replace("\\gamma","U").replace(",","}^{")
            # Increment the legend string so that we do not save it anymore
            help_variable += 1
        # Loop through our data frame and save each value (we need to cast it as a double first)
        for index in range(len(dataframe.values[:,3])):
            if index > 0:
                # Append the casted value to our data frame
                basis_functions[index-1][hole_index].append(float(dataframe.values[index,3]))            
# Colours for plotting
colour_list_for_plotting = [(115/256,115/256,115/256),(77/256,0/256,75/256), (129/256,15/256,124/256), (0/256,68/256,27/256),(0/256,109/256,44/256),(35/256,139/256,69/256), (4/256,90/256,141/256),(5/256,112/256,176/256),(54/256,144/256,192/256),(116/256,169/256,207/256),(102/256,37/256,6/256),(153/256,52/256,4/256),(204/256,76/256,2/256),(236/256,112/256,20/256),(254/256,153/256,41/256),(1/256,70/256,54/276),(1/256,108/256,89/256),(2/256,219/256,138/256),(54/256,144/256,192/256),(103/256,169/256,207/256),(166/256,189/256,219/256)]





#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the perturbed eigenfunctions
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Set all parameters to tex
plt.rcParams['text.usetex'] = True
# Figures
fig, axes = plt.subplots(2,3,figsize=(30,10))
plt.rc('axes', labelsize=35)    # fontsize of the x and y label
plt.rc('legend', fontsize=15)    # legend fontsize
plt.rc('xtick', labelsize=30)    # fontsize of the tick labels
plt.rc('ytick', labelsize=30)    # fontsize of the tick labels
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# Loop over all basis functions and plot them
for index in range(21):
    # Define where we are putting the figure
    row = 0
    col = 0
    if index ==0:
        row=0
        col=0
    elif index<3:
        row=0
        col=1
    elif index<6:
        row=0
        col=2        
    elif index<10:
        row=1
        col=0
    elif index<15:
        row=1
        col=1
    else:
        row=1
        col=2                
    # Plot the 90 percentile
    axes[row][col].plot(hole_radius_array,np.array([np.percentile(basis_functions[index][sub_index],95) for sub_index in range(len(hole_radius_array))]),'--',color=colour_list_for_plotting[index])
    # Plot the 50 percentile
    axes[row][col].plot(hole_radius_array,np.array([np.percentile(basis_functions[index][sub_index],50) for sub_index in range(len(hole_radius_array))]),'-',color=colour_list_for_plotting[index],label=legend_strings[index])
    # Plot the 5 percentile
    axes[row][col].plot(hole_radius_array,np.array([np.percentile(basis_functions[index][sub_index],5) for sub_index in range(len(hole_radius_array))]),'--',color=colour_list_for_plotting[index])
    # See if we can fill between
    axes[row][col].fill_between(hole_radius_array,np.array([np.percentile(basis_functions[index][sub_index],5) for sub_index in range(len(hole_radius_array))]), np.array([np.percentile(basis_functions[index][sub_index],95) for sub_index in range(len(hole_radius_array))]), facecolor=colour_list_for_plotting[index],alpha=0.5,interpolate=True)
    # Add the legend in the end
    axes[row][col].legend(bbox_to_anchor=(0.50,0.5), loc="center left", borderaxespad=0)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Geodesic hole radius, $\\varepsilon$")
plt.ylabel('Spectral coefficients $U_{n}^{m}$')
# displaying the title
plt.title("Spectral decomposition of $u(\\mathbf{x},t=50)$ as a function of the hole radius $\\varepsilon$",fontsize=30, fontweight='bold')
# Save fig and show it
plt.savefig("../Figures/eigenfunctions_vs_hole_radius_n_" + str(n) + ".png")
#plt.show() # Uncomment if you want a figure to pop-up

