# =================================================================================
# =================================================================================
# Script:"slope_test"
# Date: 2023-10-31
# Implemented by: Philip Gerlee
# Description:
# This script analyses the variance in the various simulations that are presented in the manuscript. All data files corresponding to these figures are stored in xlsx files that can be found in the folder "../Output/data_analysis_variance". Since we have stochasticity in the initial conditions when we run the PDE simulations, we run multiple repititions. A relevant question then is to determine whether the spread in the outputs is due solely to variations in the initial conditions or if there is actually some effect of changing a certain parameter. To this end, we do linear regression of the various simulations, and if the slope is zero that means that the observed variations are solely caused by stochasticity in the initial conditions.
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import math
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import pandas as pd
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Read all the data and analyse the variance
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Define the data folder
data_folder = "../Output/data_analysis_variance/"
#----------------------------------------------------------------------------------
# Read data: spectral (n,d)=(1,20)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_5_spectral_n_1_d_20.xlsx', sheet_name='$U_{1}^{0}$')
df['mean'] = df.mean(axis=1)
#----------------------------------------------------------------------------------
# Plot data: spectral (n,d)=(1,20)
#----------------------------------------------------------------------------------
# Set all parameters to tex
plt.rcParams['text.usetex'] = True
# Define the figure and define all label sizes etc.
fig, axes = plt.subplots(figsize=(30,10))
plt.rc('axes', labelsize=35)    # fontsize of the x and y label
plt.rc('legend', fontsize=15)    # legend fontsize
plt.rc('xtick', labelsize=30)    # fontsize of the tick labels
plt.rc('ytick', labelsize=30)    # fontsize of the tick labels
axes.plot(df['Unnamed: 0'],df['mean'],'.',color=(0,0,0))
plt.xlabel("Geodesic hole radius, $\\varepsilon$",fontsize=20)
plt.ylabel('Mean value of spectral coefficients $U_{1}^{0}$',fontsize=20)
# displaying the title
plt.title("Spectral decomposition of $u(\\mathbf{x},t=50)$ as a function of the hole radius $\\varepsilon$",fontsize=30, fontweight='bold')
plt.show()
#----------------------------------------------------------------------------------
# Coefficient U_{1}^{0} for (n,d)=(1,20)
#----------------------------------------------------------------------------------
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for $U_{1}^{0}$: ",p_value)
#----------------------------------------------------------------------------------
# Coefficient U_{1}^{1} for (n,d)=(1,20)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_5_spectral_n_1_d_20.xlsx', sheet_name='$U_{1}^{1}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for $U_{1}^{1}$: ",p_value)
#----------------------------------------------------------------------------------
# Coefficient U_{2}^{0} for (n,d)=(2,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_6_spectral_n_2_d_18.xlsx', sheet_name='$U_{2}^{0}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for $U_{2}^{0}$: ",p_value)
#----------------------------------------------------------------------------------
# Coefficient U_{2}^{1} for (n,d)=(2,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_6_spectral_n_2_d_18.xlsx', sheet_name='$U_{2}^{1}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for $U_{2}^{1}$: ",p_value)
#----------------------------------------------------------------------------------
# Coefficient U_{2}^{2} for (n,d)=(2,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_6_spectral_n_2_d_18.xlsx', sheet_name='$U_{2}^{2}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for $U_{2}^{2}$: ",p_value)
#----------------------------------------------------------------------------------
# Metrics for (n,d)=(1,20)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_7_metrics_n_1_d_20.xlsx', sheet_name='Minimal distance')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for 'Minimal distance (n,d)=(1,20)': ",p_value)
#----------------------------------------------------------------------------------
# Metrics for (n,d)=(2,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_8_metrics_n_2_d_18.xlsx', sheet_name='Minimal distance')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for 'Minimal distance (n,d)=(2,18)': ",p_value)
#----------------------------------------------------------------------------------
# U_{3}^{0} for (n,d)=(3,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_S7_spectral_n_3_d_18.xlsx', sheet_name='$U_{3}^{0}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{3}^{0}$': ",p_value)
#----------------------------------------------------------------------------------
# U_{3}^{1} for (n,d)=(3,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_S7_spectral_n_3_d_18.xlsx', sheet_name='$U_{3}^{1}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{3}^{1}$': ",p_value)
#----------------------------------------------------------------------------------
# U_{3}^{2} for (n,d)=(3,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_S7_spectral_n_3_d_18.xlsx', sheet_name='$U_{3}^{2}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{3}^{2}$': ",p_value)
#----------------------------------------------------------------------------------
# U_{3}^{3} for (n,d)=(3,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_S7_spectral_n_3_d_18.xlsx', sheet_name='$U_{3}^{3}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{3}^{3}$': ",p_value)
#----------------------------------------------------------------------------------
# U_{4}^{0} for (n,d)=(4,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_S10_spectral_n_4_d_18.xlsx', sheet_name='$U_{4}^{0}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{4}^{0}$': ",p_value)
#----------------------------------------------------------------------------------
# U_{4}^{1} for (n,d)=(4,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_S10_spectral_n_4_d_18.xlsx', sheet_name='$U_{4}^{1}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{4}^{1}$': ",p_value)
#----------------------------------------------------------------------------------
# U_{4}^{2} for (n,d)=(4,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_S10_spectral_n_4_d_18.xlsx', sheet_name='$U_{4}^{2}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{4}^{2}$': ",p_value)
#----------------------------------------------------------------------------------
# U_{4}^{3} for (n,d)=(4,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_S10_spectral_n_4_d_18.xlsx', sheet_name='$U_{4}^{3}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{4}^{3}$': ",p_value)
#----------------------------------------------------------------------------------
# U_{4}^{4} for (n,d)=(4,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_S10_spectral_n_4_d_18.xlsx', sheet_name='$U_{4}^{4}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{4}^{4}$': ",p_value)
#----------------------------------------------------------------------------------
# Metrics for (n,d)=(3,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_S8_metrics_n_3_d_18.xlsx', sheet_name='Minimal distance')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for 'Minimal distance (n,d)=(3,18)': ",p_value)
#----------------------------------------------------------------------------------
# Metrics for (n,d)=(4,18)
#----------------------------------------------------------------------------------
df = pd.read_excel(data_folder + 'Fig_S11_metrics_n_4_d_18.xlsx', sheet_name='Minimal distance')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for 'Minimal distance (n,d)=(4,18)': ",p_value)






