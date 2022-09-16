# =================================================================================
# =================================================================================
# Script:"find_parameters_and_plot_perturbed_eigenvalues"
# Date: 2021-03-30
# Implemented by: Johannes Borgqvist
# Description:
# This script checks whether certain parameters satisfy the Turing conditions and
# then it plots the eigenvalues as functions of the hole radius epsilon. We can refer
# to this script as Philip's spaghetti plots. 
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
n = 1
#n = 2
#n = 4 
k_squared = n*(n+1)
# Calculate the steady states and the critical parameters
u_0, v_0, d_c, gamma_c = Schnakenberg_properties.calculate_steady_states_and_critical_parameters_Schnakenberg(a,b,k_squared)
# Save the steady states in a list
steady_states = [u_0,v_0]
# Set the value of the relative diffusion
#d = d_c + 1.5
d = 20# n=1
#d=18 # n=2 och n=4
# Set the value of the reaction strength to its critical value
gamma = gamma_c
# Calculate the critical hole radius for all eigenvalues between n=2 and n=4
#n_largest = 6
n_largest = 2
eps_tuple, n_tuple, m_tuple = Schnakenberg_properties.compute_minimal_holeradius_for_pattern_disturbance(a,b,d,gamma,n,n_largest)
# Extract the minimal radii
eps_min = eps_tuple[0]
n_min = n_tuple[0]
m_min = m_tuple[0]
# Extract the maximal radii
eps_max = eps_tuple[1]
n_max = n_tuple[1]
m_max = m_tuple[1]
# Calculate the critical radius
#hole_cylinder_radius = np.sin(eps_min)
hole_cylinder_radius = np.sin(eps_max)
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
eigen_value_list = [(n,m) for n in [1,2,3,4] for m in range(n+1)]
print(eigen_value_list)
# Find some really nice colours for the curves based on the colour maps on https://colorbrewer2.org/#type=sequential&scheme=PuBu&n=9
#colour_list_for_plotting = [(77/256,0/256,75/256), (129/256,15/256,124/256), (136/256,65/256,157/256),(0/256,68/256,27/256),(0/256,109/256,44/256),(35/256,139/256,69/256),(65/256,174/256,118/256),(2/256,56/256,88/256),(4/256,90/256,141/256),(5/256,112/256,176/256),(54/256,144/256,192/256),(116/256,169/256,207/256)]
#colour_list_for_plotting = [(77/256,0/256,75/256), (129/256,15/256,124/256), (0/256,68/256,27/256),(0/256,109/256,44/256),(35/256,139/256,69/256), (4/256,90/256,141/256),(5/256,112/256,176/256),(54/256,144/256,192/256),(116/256,169/256,207/256)]
colour_list_for_plotting = [(77/256,0/256,75/256), (129/256,15/256,124/256), (0/256,68/256,27/256),(0/256,109/256,44/256),(35/256,139/256,69/256), (4/256,90/256,141/256),(5/256,112/256,176/256),(54/256,144/256,192/256),(116/256,169/256,207/256),(102/256,37/256,6/256),(153/256,52/256,4/256),(204/256,76/256,2/256),(236/256,112/256,20/256),(254/256,153/256,41/256)]
# Create a list of all labels as well
label_strings = ["$\\lambda_{" + str(eigen_value_list[index][0]) + "," + str(eigen_value_list[index][1]) + "}$" for index in range(len(eigen_value_list))]
# Create an np array called epsilon vector with hole radii
#epsilon_vector = np.linspace(0,np.sin(1.0),100,endpoint=True)
epsilon_vector = np.linspace(0,0.7,100,endpoint=True)
# Create a list of np arrays with the corresponding eigenvalues
lambda_vec = [np.array([Schnakenberg_properties.perturbed_eigenvalue_Schnakenberg(nm_tuple[0],nm_tuple[1],epsilon) for epsilon in list(epsilon_vector)]) for nm_tuple in eigen_value_list]
# Create vector for the lower and upper bounds as well
upper_bound = np.ones(len(epsilon_vector))*gamma*M
lower_bound = np.ones(len(epsilon_vector))*gamma*L
# Create an index list which we will loop over subsequently
index_list = list(range(len(eigen_value_list)))
# Reverse this list so that the legend looks nice
index_list.sort(reverse=True)
# Generate a plottable vertical line for the eigenvalues as well
#print("\t\tThe critical radius:\t(n,m,r)\t=\t(%d,%d,%0.3f)"%(n_min,m_min,hole_cylinder_radius))
print("\t\tThe critical radius:\t(n,m,r)\t=\t(%d,%d,%0.3f)"%(n_max,m_max,hole_cylinder_radius))
vertical_line = np.linspace(0,21,50)
crit_radius = hole_cylinder_radius*np.ones(len(vertical_line))
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
axes.plot(epsilon_vector,upper_bound,'--k',label="$\\gamma\\;M$")
axes.plot(epsilon_vector,lower_bound,'--k',label="$\\gamma\\;L$")
axes.plot(crit_radius,vertical_line,'+k',label="$\\varepsilon_{\\mathrm{crit}}$")
axes.legend()
#hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
#plt.xlabel("Geodesic hole radius, $\\varepsilon$")
plt.xlabel("Cylindrical hole radius, $\\varepsilon$")
plt.ylabel("Eigenvalues, $\\lambda_{n,m}(\\varepsilon)$")
# displaying the title
plt.title("Perturbed eigenvalues $\\lambda_{n,m}(\\varepsilon)$ as a function of the hole radius $\\varepsilon$",fontsize=30, fontweight='bold')
plt.show()
plt.savefig("../Figures/perturbed_eigenvalues.png")



#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# ILLUSTRATE THE EIGENVALUES IN LATEX
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Plot the eigenvalues
for index in index_list:
    # We need marks to distinguish between the various cases
    mark_str = " every mark/.append style={solid, fill=eigen_" + str(eigen_value_list[index][0]) + "_" + str(eigen_value_list[index][1])
    # Decide the type of mark based on the value of m
    if eigen_value_list[index][1]==0:
        mark_str += "}, mark=*, "
    elif eigen_value_list[index][1]==1:
        mark_str += "}, mark=square*, "
    elif eigen_value_list[index][1]==2:
        mark_str += "}, mark=otimes*, "
    elif eigen_value_list[index][1]==3:
        mark_str += "}, mark=triangle*, "
    elif eigen_value_list[index][1]==4:
        mark_str +=  "}, mark=diamond*, "
    elif eigen_value_list[index][1]==5:
        mark_str +=  "}, mark=star, "
    # Plot the perturbed eigenvalue    
    plot_LaTeX_2D(epsilon_vector,lambda_vec[index],"../Figures/illustrate_eigenvalues/Input/n_" + str(n) + "_d_" + str(d) + ".tex","color=eigen_" + str(eigen_value_list[index][0]) + "_" + str(eigen_value_list[index][1]) + ",line width=0.5pt,"+mark_str+"mark size=0.75pt",label_strings[index])

# Plot the thresholds
plot_LaTeX_2D(epsilon_vector,upper_bound,"../Figures/illustrate_eigenvalues/Input/n_" + str(n) + "_d_" + str(d) + ".tex","only marks, mark=halfcircle*,mark size=1.0pt,color=black,","$\gamma M(a,b,d)$")
plot_LaTeX_2D(epsilon_vector,lower_bound,"../Figures/illustrate_eigenvalues/Input/n_" + str(n) + "_d_" + str(d) + ".tex","only marks, mark=diamond*,mark size=1.0pt,color=black,","$\gamma L(a,b,d)$")
# Plot the vertical line as well
#plot_LaTeX_2D(crit_radius,vertical_line,"../Figures/illustrate_eigenvalues/Input/perturbed_eigenvalues.tex","only marks, mark=diamond*,mark size=0.5pt,color=black,","$\\varepsilon_{\\mathrm{crit}}$")
print("\n\n==============================================================================================================================\n")
print("\t Looking at the parameter space\n")
print("==============================================================================================================================\n")


num_cols = 77
a_vec = np.linspace(0.01,0.25,num_cols,endpoint=True)
d_vec = np.linspace(5,25,num_cols,endpoint=True)

dist_vec = [(a_temp,d_temp,Schnakenberg_properties.calculate_distance_between_bounds(a_temp,b,d_temp)) for d_temp in list(d_vec) for a_temp in list(a_vec)]



plot_LaTeX_3D(dist_vec,"../Figures/parameter_plot_Schnakenberg/Input/parameter_surface.tex","surf","Turing region",True)

d_vec = np.array([Schnakenberg_properties.calculate_steady_states_and_critical_parameters_Schnakenberg(a_temp,b,k_squared)[2] for a_temp in list(a_vec)])


dist_vec_crit = [(a_vec[index],d_vec[index],0) for index in range(len(a_vec))]


plot_LaTeX_3D(dist_vec_crit,"../Figures/parameter_plot_Schnakenberg/Input/critical_diffusion.tex","color=black,line width=2pt,mark=none,","Critical diffusion $d_c$",False)
