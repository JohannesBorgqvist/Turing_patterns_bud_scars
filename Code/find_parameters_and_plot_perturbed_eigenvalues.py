# =================================================================================
# =================================================================================
# Script:"find_parameters_and_plot_perturbed_eigenvalues"
# Date: 2021-03-08
# Implemented by: Johannes Borgqvist
# Description:
# This script checks whether certain parameters satisfy the Turing conditions and
# then it plots the eigenvalues as functions of the hole radius epsilon
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import Schnakenberg_properties # Home made
import numpy as np # For numerical calculations
from matplotlib import pyplot as plt # For plotting
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
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Define the parameter pairs
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Experiment 1 in Chaplain constitute our reference parameters:
# a = 0.2
# b = 1
# d = 18
# gamma = 20.62
#----------------------------------------------------------------------------------
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
# Set the value of the relative diffusion
d = 18
# Set the value of the reaction strength to its critical value
gamma = gamma_c
# Prompt to the user
print("\n\n==============================================================================================================================\n")
print("\t Testing parameters and plotting perturbed eigenvalues\n")
print("==============================================================================================================================\n")
# Print the results
print("\n\t\tThe parameters are:\t\t(a,b,gamma,d)\t=\t(%0.3f,%0.3f,%0.3f,%0.3f)"%(a,b,gamma,d))
print("\t\tThe steady states:\t\t(u_0,v_0)\t=\t(%0.4f,%0.4f)"%(u_0,v_0))
print("\t\tThe critical parameters:\t(d_c,gamma_c)\t=\t(%0.4f,%0.4f)"%(d_c,gamma_c))
# Calculate the Turing conditions
Turing_conditions,L,M = Schnakenberg_properties.check_Turing_conditions_Scnakenberg(a,b,d)
# Promt to the user
if Turing_conditions:
    print("\n\t\tThe Turing conditions were satisfied! The lower and upper bounds are:")
    print("\t\t\t Upper:\t\t\tgamma M\t\t=\t%0.3f"%(gamma*M))    
    print("\t\t\t Lower:\t\t\tgamma L\t\t=\t%0.3f\n\n"%(gamma*L))
else:
    print("\n\t\tThe Turing conditions were not satisfied.\n\n")
# Create a list of all tuples of m and n values
eigen_value_list = [(n,m) for n in [2,3,4] for m in range(n+1)]
print(eigen_value_list)
# Find some really nice colours for the curves based on the colour maps on https://colorbrewer2.org/#type=sequential&scheme=PuBu&n=9
colour_list_for_plotting = [(77/256,0/256,75/256), (129/256,15/256,124/256), (136/256,65/256,157/256),(0/256,68/256,27/256),(0/256,109/256,44/256),(35/256,139/256,69/256),(65/256,174/256,118/256),(2/256,56/256,18/256),(4/256,90/256,141/256),(5/256,112/256,176/256),(54/256,144/256,192/256),(116/256,169/256,207/256)]
# Create a list of all labels as well
label_strings = ["$\\lambda_{" + str(eigen_value_list[index][0]) + "," + str(eigen_value_list[index][1]) + "}$" for index in range(len(eigen_value_list))]
# Create an np array called epsilon vector with hole radii
epsilon_vector = np.linspace(0,0.50,100,endpoint=True)
# Create a list of np arrays with the corresponding eigenvalues
lambda_vec = [np.array([Schnakenberg_properties.perturbed_eigenvalue_Schnakenberg(nm_tuple[0],nm_tuple[1],epsilon) for epsilon in list(epsilon_vector)]) for nm_tuple in eigen_value_list]
# Create vector for the lower and upper bounds as well
upper_bound = np.ones(len(epsilon_vector))*gamma*M
lower_bound = np.ones(len(epsilon_vector))*gamma*L
# Create an index list which we will loop over subsequently
index_list = list(range(len(eigen_value_list)))
# Reverse this list so that the legend looks nice
index_list.sort(reverse=True)

print("\n\n==============================================================================================================================")

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the eigenvalues as a function 
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
fig, axes = plt.subplots(1,1,figsize=(15,5))
plt.rc('axes', labelsize=25)    # fontsize of the x and y label
plt.rc('legend', fontsize=20)    # legend fontsize
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# Plot the eigenvalues
for index in index_list:
    axes.plot(epsilon_vector,lambda_vec[index],'-',color=colour_list_for_plotting[index],label=label_strings[index])
# Plot the upper bound
axes.plot(epsilon_vector,upper_bound,'--sk',label="$\\gamma\\;M$")
axes.plot(epsilon_vector,lower_bound,'--pk',label="$\\gamma\\;L$")
axes.legend()
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Geodesic hole radius, $\\varepsilon$")
plt.ylabel("Eigenvalues, $\\lambda_{n,m}(\\varepsilon)$")
# displaying the title
plt.title("Perturbed eigenvalues $\\lambda_{n,m}(\\varepsilon)$ as a function of the hole radius $\\varepsilon$",fontsize=30, fontweight='bold')
plt.show()
#plt.savefig("../Figures/perturbed_eigenvalues.png")



#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# ILLUSTRATE THE EIGENVALUES IN LATEX
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the eigenvalues
for index in index_list:
    plot_LaTeX_2D(epsilon_vector,lambda_vec[index],"../Figures/illustrate_eigenvalues/Input/perturbed_eigenvalues.tex","color=eigen_" + str(eigen_value_list[index][0]) + "_" + str(eigen_value_list[index][1]) + ",line width=2pt,",label_strings[index])

