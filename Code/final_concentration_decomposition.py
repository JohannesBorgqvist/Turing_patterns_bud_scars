# =================================================================================
# =================================================================================
# Script:"final_concentration_decomposition"
# Date: 2022-08-15
# Implemented by: Johannes Borgqvist
# Description:
# This script reads the final concentration profile and decomposes it in terms
# of the eigenfunctions of the Laplace Beltrami operator
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
#n = 5
n = 4
k_squared = n*(n+1)
# Calculate the steady states and the critical parameters
u_0, v_0, d_c, gamma_c = Schnakenberg_properties.calculate_steady_states_and_critical_parameters_Schnakenberg(a,b,k_squared)
# Save the steady states in a list
steady_states = [u_0,v_0]
# Set the value of the relative diffusion
d = 18.0 # n=2,4,5
#d = 20.0 # n=1
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
#hole_radius_array = np.arange(0,0.75,0.05)
hole_radius_array = np.array([0])
# Allocate a list of all the basis functions
basis_functions = []
# Let's add 20 lists for each basis function corresponding to
# n = 1,2,3,4,5...
for index in range(21):
    basis_functions.append([])   
# Now for each hole radius add an empty list as well with all iterations
#for index in range(14):
    #for sub_index in range(len(hole_radius_array)):
        #basis_functions[index].append([])
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
    #for repitition_index in range(10):
    for repitition_index in range(1):
        # Gather all these substrings into one giant string where we will save the output files
        output_folder_str = folder_str + hole_str + radius_str + a_str + b_str + d_str + gamma_str + sigma_str + T_str + IC_str + "iteration_" + str(repitition_index) + "/"
        print(output_folder_str)
        # Read the csv file
        dataframe = read_csv(output_folder_str + "spectral_coefficients.csv", header=None)
        # Save the legends one time
        if help_variable == 0:
            # Save the legend strings
            legend_strings = list(dataframe.values[:,1])
            # Remove the first value
            del legend_strings[0]
            # Loop over all legends and change their name
            for index,legend in enumerate(legend_strings):
                legend_strings[index] = legend_strings[index].replace("\\gamma","U")
            # Increment the legend string so that we do not save it anymore
            help_variable += 1
        # Loop through our data frame and save each value (we need to cast it as a double first)
        for index in range(len(dataframe.values[:,3])):
            if index>0:
                # Append the casted value to our data frame
                basis_functions[index-1].append(float(dataframe.values[index,3]))

# Save indices            
basis_functions = [value[0] for value in basis_functions]
basis_function_index = [index for index,value in enumerate(basis_functions)]
print(basis_function_index)



plot_LaTeX_2D(basis_function_index,basis_functions,"../Figures/n_" + str(n) + "_final_concentration_decomposition/Input/spectral_analysis.tex","only marks,mark size=1.5pt,color=clr_1,",[])

            
