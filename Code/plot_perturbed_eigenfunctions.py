# =================================================================================
# =================================================================================
# Script:"plot_perturbed_eigenfunctions"
# Date: 2021-04-29
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
# =================================================================================
# =================================================================================
# Functions
# =================================================================================
# =================================================================================
def plot_LaTeX_2D(t,y,file_str,plot_str,legend_str):
    # Open a file with the append option
    # so that we can write to the same
    # file multiple times
    f = open(file_str, "a")
    # Create a temporary string which
    # is the one that does the plotting.
    # Here we incorporate the input plot_str
    # which contains the color, and the markers
    # of the plot at hand
    if len(legend_str)==0:
        temp_str = "\\addplot[\nforget plot,\n" + plot_str+ "\n]\n"
    else:
        temp_str = "\\addplot[\n" + plot_str+ "\n]\n"
    # Add the coordinates
    temp_str += "coordinates {%\n"
    # Loop over the input files and add
    # them to the file
    for i in range(len(t)):
        temp_str += "(" + str(t[i]) + "," + str(y[i]) + ")\n"
    # The plotting is done, let's close the shop    
    temp_str += "};\n"
    # Add a legend if one is provided
    if len(legend_str) > 0:
        temp_str += "\\addlegendentry{" + legend_str + "}\n"
    # Finally, we write the huge string
    # we have created
    f.write("%s"%(temp_str))
    # Close the file
    f.close()
def plot_LaTeX_3D(data,file_str,plot_str,legend_str,surfaceNotCurve):
    # Open a file with the append option
    # so that we can write to the same
    # file multiple times
    f = open(file_str, "a")
    # Create a temporary string which
    # is the one that does the plotting.
    # Here we incorporate the input plot_str
    # which contains the color, and the markers
    # of the plot at hand
    if surfaceNotCurve:
        if len(legend_str)==0:
            temp_str = "\\addplot3[forget plot," + plot_str+ "]\n"
        else:
            temp_str = "\\addplot3[" + plot_str+ "]\n"
    else:
        if len(legend_str)==0:
            temp_str = "\\addplot3+[forget plot," + plot_str+ "]\n"
        else:
            temp_str = "\\addplot3+[" + plot_str+ "]\n"        
    # Add the coordinates
    temp_str += "coordinates {%\n"
    # Loop over the input files and add
    # them to the file
    for index in range(len(data)):
        temp_str += "("+str(data[index][0]) + "," + str(data[index][1]) + "," + str(data[index][2]) + ")"
        if index>0:
            if index < len(data)-1 and data[index][1] < data[index+1][1]:
                temp_str += "\n"
            elif index == len(data)-1:
                temp_str += "\n"
            else:
                temp_str += "  "
    # The plotting is done, let's close the shop    
    temp_str += "};\n"
    # Add a legend if one is provided
    if len(legend_str) > 0:
        temp_str += "\\addlegendentry{" + legend_str + "}\n"
    # Finally, we write the huge string
    # we have created
    f.write("%s"%(temp_str))
    # Close the file
    f.close()    
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
# Define the meshes we want to loop over
hole_radius_array = np.arange(0,0.75,0.05)
# Allocate a list of all the basis functions
basis_functions = []
# Let's add 14 lists for each basis function
for index in range(14):
    basis_functions.append([])   
# Now for each hole radius add an empty list as well with all iterations
for index in range(14):
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
sigma_str = "sigma_" + str(round(sigma,3)).replace(".","p") + "_"
T_str = "T_" + str(round(T,3)).replace(".","p") + "_"
if ICs_around_steady_states:
    IC_str = "ICs_around_steady_states/"
else:
    IC_str = "ICs_at_zero/"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Allocate memory for the legend_strings
legen_strings = []
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
    for repitition_index in range(10):
        # Gather all these substrings into one giant string where we will save the output files
        output_folder_str = folder_str + hole_str + radius_str + a_str + b_str + d_str + gamma_str + sigma_str + T_str + IC_str + "iteration_" + str(repitition_index) + "/"
        # Read the csv file
        dataframe = read_csv(output_folder_str + "spectral_coefficients.csv", header=None)
        # Save the legends one time
        if help_variable == 0:
            # Save the legend strings
            legend_strings = list(dataframe.values[1:-1,1])
            # Increment the legend string so that we do not save it anymore
            help_variable += 1
        # Loop through our data frame and save each value (we need to cast it as a double first)
        for index in range(len(dataframe.values[1:-1,3])):
            # Append the casted value to our data frame
            basis_functions[index][hole_index].append(float(dataframe.values[index+1,3]))            

colour_list_for_plotting = [(115/256,115/256,115/256),(77/256,0/256,75/256), (129/256,15/256,124/256), (0/256,68/256,27/256),(0/256,109/256,44/256),(35/256,139/256,69/256), (4/256,90/256,141/256),(5/256,112/256,176/256),(54/256,144/256,192/256),(116/256,169/256,207/256),(102/256,37/256,6/256),(153/256,52/256,4/256),(204/256,76/256,2/256),(236/256,112/256,20/256),(254/256,153/256,41/256)]
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the perturbed eigenfunctions
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
fig, axes = plt.subplots(1,1,figsize=(15,5))
plt.rc('axes', labelsize=25)    # fontsize of the x and y label
plt.rc('legend', fontsize=20)    # legend fontsize
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# Loop over all basis functions and plot them
for index in range(14):
    # Plot the 90 percentile
    axes.plot(hole_radius_array,np.array([np.percentile(basis_functions[index][sub_index],95) for sub_index in range(len(hole_radius_array))]),'--',color=colour_list_for_plotting[index])
    # Plot the 50 percentile
    axes.plot(hole_radius_array,np.array([np.percentile(basis_functions[index][sub_index],50) for sub_index in range(len(hole_radius_array))]),'-',color=colour_list_for_plotting[index],label=legend_strings[index])
    # Plot the 5 percentile
    axes.plot(hole_radius_array,np.array([np.percentile(basis_functions[index][sub_index],5) for sub_index in range(len(hole_radius_array))]),'--',color=colour_list_for_plotting[index])
    # See if we can fill between
    plt.fill_between(hole_radius_array,np.array([np.percentile(basis_functions[index][sub_index],5) for sub_index in range(len(hole_radius_array))]), np.array([np.percentile(basis_functions[index][sub_index],95) for sub_index in range(len(hole_radius_array))]), facecolor=colour_list_for_plotting[index],alpha=0.5,interpolate=True)
# Add the legend in the end
axes.legend(bbox_to_anchor=(1.00,0.5), loc="center left", borderaxespad=0)
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
#plt.xlabel("Geodesic hole radius, $\\varepsilon$")
plt.xlabel("Cylindrical hole radius, $\\varepsilon$")
plt.ylabel("Eigenfunctions, $\\lambda_{n,m}(\\varepsilon)$")
# displaying the title
plt.title("Perturbed eigenfunctions $\\gamma_{n,m}(\\varepsilon)$ as a function of the hole radius $\\varepsilon$",fontsize=30, fontweight='bold')
plt.show()
plt.savefig("../Figures/perturbed_eigenfunctions.png")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# ILLUSTRATE THE EIGENVALUES IN LATEX
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Create a list of all tuples of m and n values
eigen_value_list = [(n,m) for n in [0,1,2,3,4] for m in range(n+1)]
# Plot the eigenvalues
# Loop over all basis functions and plot them
for index in range(14):
#for index in range(3,6):
    plot_LaTeX_2D(hole_radius_array,np.array([np.percentile(basis_functions[index][sub_index],95) for sub_index in range(len(hole_radius_array))]),"../Figures/illustrate_eigenfunctions/Input/perturbed_eigenfunctions.tex","forget plot, densely dashed, thin, color=eigen_" + str(eigen_value_list[index][0]) + "_" + str(eigen_value_list[index][1]) + ",line width=0.2pt,name path=up_n_" + str(eigen_value_list[index][0]) + "_m_" + str(eigen_value_list[index][1]) + ",",[])
    plot_LaTeX_2D(hole_radius_array,np.array([np.percentile(basis_functions[index][sub_index],50) for sub_index in range(len(hole_radius_array))]),"../Figures/illustrate_eigenfunctions/Input/perturbed_eigenfunctions.tex","forget plot, densely dashed, thin, color=eigen_" + str(eigen_value_list[index][0]) + "_" + str(eigen_value_list[index][1]) + ",line width=1pt,",[])
    plot_LaTeX_2D(hole_radius_array,np.array([np.percentile(basis_functions[index][sub_index],5) for sub_index in range(len(hole_radius_array))]),"../Figures/illustrate_eigenfunctions/Input/perturbed_eigenfunctions.tex","color=eigen_" + str(eigen_value_list[index][0]) + "_" + str(eigen_value_list[index][1]) + ",line width=0.2pt,name path=down_n_" + str(eigen_value_list[index][0]) + "_m_" + str(eigen_value_list[index][1]) + ",",legend_strings[index])    
