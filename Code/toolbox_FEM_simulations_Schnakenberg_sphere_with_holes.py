# =================================================================================
# =================================================================================
# Script:"toolbox_FEM_simulations_Schnakenberg_sphere_with_holes"
# Date: 2022-03-21
# Implemented by: Johannes Borgqvist and Carl Lundholm
# Description:
# This is the main script containing all important functions needed in order to
# run the FEM simulations of the Schnakenberg model on the sphere with holes.
# ### NOTE ###
# This script uses the finite difference scheme 1-SBEM in time, resulting in two
# (n x n) systems of equations. The code has been hard coded to be more efficient
# for precisely the scheme 1-SBEM. To play around with and test other types of
# finite difference schemes it's recommended to use the script
# "toolbox_FEM_simulations_Schnakenberg_sphere_with_holes_MixedSystems".
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
# Most good things are contained in FEniCS
from fenics import *
# Interpolate non_matching_mesh
#from fenicstools import interpolate_nonmatching_mesh
# Import numpy as well
import numpy as np
# Write less trick
np.fac = np.math.factorial
# Import pandas as well
import pandas as pd 
# Import scipy to load the initial conditions nicely
#from scipy.io import loadmat, savemat
# =================================================================================
# =================================================================================
# Functions for conducting the FEM simulations
# =================================================================================
# =================================================================================
#------------------------------------------------------------------
# Function 1: "read_mesh_Schnakenberg_sphere_with_holes"
# The function takes the number of holes as inputs and
# returns the following four outputs:
# 1. The number of holes in a variable called "num_holes",
# 2. A list of all the radii of the holes called "radii_holes",
# 3. A mesh_function "mf_subdomain",
# 4. A list of integration measures "dx_list" used in the
# variational formulation.
#------------------------------------------------------------------

def read_mesh_Schnakenberg_sphere_with_holes(num_holes,radii_holes):
    # Allocate memory for the mesh and the mesh value collection
    mesh = Mesh()
    mvc_subdomains = MeshValueCollection("size_t", mesh, 2)
    # Define a mesh function taking the subdomains into account
    mf_subdomains = 0
    # Allocate memory for a list containing all the integration
    # measures involved in the variational formulation
    dx_list = []    
    # Define the string of the mesh that we want to read
    mesh_str = "../Meshes/s_h_" + str(num_holes)
    # Define the string in which we read the mesh
    # depending on the number of holes
    if num_holes > 0:
        # Loop over the radii and add these to the string name of the mesh
        for radius in radii_holes:
            mesh_str += "_r_" + str(round(radius,3)).replace(".","p") + "_"
    # Now, we add the file format to the string as well
    mesh_str += ".xdmf"
    # And tidy the file name up in necessary
    mesh_str = mesh_str.replace("_.xdmf",".xdmf")
    print(mesh_str)
    # Read in the mesh and the subdomains into these two variables
    with XDMFFile(mesh_str) as infile:
        infile.read(mesh)
        infile.read(mvc_subdomains, "name_to_read")
    # Define a mesh function taking the subdomains into account
    mf_subdomains = cpp.mesh.MeshFunctionSizet(mesh, mvc_subdomains)
    # Add our integration measure to the list of integration measure
    dx_list.append(Measure("dx", domain=mesh, subdomain_data=mf_subdomains, subdomain_id=1))
    # Lastly, return our mesh, the mesh value collection, the mesh
    # function and the integration measures.
    return mesh, mvc_subdomains, mf_subdomains, dx_list

#------------------------------------------------------------------
# Function 2: "set_up_FEM_Schnakenberg_sphere_with_holes"
# The function sets up the FEM framework for the Schnakenberg
# model on the sphere with holes. It takes a mesh as input
# and then it returns the function space H.
#------------------------------------------------------------------

def define_Hilbert_space_Schnakenberg_sphere_with_holes(mesh):

    # Define the basis finite elements and their order. We pick linear basis functions. 
    P1 = FiniteElement("P", mesh.ufl_cell(), 1)

    # Define a mixed element since we have a coupled system of two PDEs
    element = MixedElement([P1, P1])

    # Collect everything in the function space (name it H for a Hilbert space)
    H = FunctionSpace(mesh, element)

    # Return all of these elements
    return H

#------------------------------------------------------------------
# Function 3: "perturbations_ic__Schnakenberg_sphere_with_holes"
# The function is a help function for the Function 5 called "initial_conditions_Schnakenberg_sphere_with_holes" which sets the initial conditions of the system. It generates a vector where each element corresponds to a node in the mesh of interest, and each element in this vector corresponds to a perturbation error drawn from a standard normal distribution. The inputs are the vector itself called epsilon and the parameter sigma determining the variance of the standard normal distribution from which the perturbations of the steady states in the initial conditions are drawn.
#------------------------------------------------------------------

def perturbations_ic__Schnakenberg_sphere_with_holes(epsilon,sigma):
    #epsilon.vector().set_local(np.random.normal(0,sigma,epsilon.vector().local_size()))
    # We use a uniform random distribution between -sigma and sigma
    epsilon.vector().set_local(np.random.uniform(-sigma,sigma,epsilon.vector().local_size()))
    return epsilon

#------------------------------------------------------------------
# Function 4: "initial_conditions_Schnakenberg_sphere_with_holes"
# The function sets the initial conditions for the system and writes these two the function u_prev which is given as an input. The initial conditions for the two states in the Schnakenberg are set to a small perturbation around the steady state. The perturbation is given by a normal distribution with zero mean and where the variance is determined by the input sigma. Hence, no output is returned and the following input must be provided:
# 1. The Hilbert space H,
# 2. The mesh stored in the variable mesh,
# 3. The mesh function mf_subdomains indicating the subdomains in the mesh,
# 4. The number of holes on the sphere stored in num_holes,
# 5. The steady states of the Schnakenberg model stored in steady_states,
# 6. The variance of the perturbation determined by sigma,
# 7. The FEM functions on the function space H stored in u and v in which we set the initial conditions,
# 8. A logical variable called ICs_around_steady_states which determines whether the initial conditions are set to the steady states values or not. 
#------------------------------------------------------------------

def initial_conditions_Schnakenberg_sphere_with_holes(H,mesh,mf_subdomains,num_holes,steady_states,sigma,u,v,ICs_around_steady_states):
    #----------------------------------------------------------------
    # DETAILED DESCRIPTION OF WHAT THE FUNCTION DOES
    # We assign the following initial conditions to the sphere:
    #1. u(t=0,*)=u0+epsilon
    #2. v(t=0,*)=v0+epsilon
    # where u0 and v0 are the provided steady states stored in the list steady_states and epsilon is a uniformly distributed variable that is drawn from a uniform distribution with variance determined by the input sigma. 
    #----------------------------------------------------------------
    # Create a dofmap
    dofmap = H.dofmap()
    # Classify the dofs as either in the hole or on the rest of the sphere
    sphere_dofs = []
    # Loop over the cells in the mesh and add the dofs in the respective list
    for cell in cells(mesh):
        # Add all the dofs in the sphere
        sphere_dofs.extend(dofmap.cell_dofs(cell.index()))            
    # Find the unique dofs
    sphere_dofs = list(set(sphere_dofs))
    # Access the component steady states as well
    u0 = steady_states[0]
    v0 = steady_states[1]    
    # Define an expression for the initial conditions based on the steady states
    if ICs_around_steady_states:
        ic_u = Expression("ic_u",ic_u=u0,degree=1)
        ic_v = Expression("ic_v",ic_v=v0,degree=1)
    else:
        ic_u = Expression("ic_u",ic_u=0,degree=1)
        ic_v = Expression("ic_v",ic_v=0,degree=1)        
    # Interpolate onto the function space
    ss_u = interpolate(ic_u,H)
    ss_v = interpolate(ic_v,H)
    # Define a function for the random perturbations
    epsilon_u = Function(H)
    epsilon_v = Function(H)
    # Fill this vector with a bunch of errors drawn
    # from a standard normal distribution with variance sigma
    epsilon_u = perturbations_ic__Schnakenberg_sphere_with_holes(epsilon_u,sigma)
    epsilon_v = perturbations_ic__Schnakenberg_sphere_with_holes(epsilon_v,sigma)
    # Assign the value to our input vector U which will correspond to the initial conditions
    u.vector()[sphere_dofs]=ss_u.vector()[sphere_dofs]+epsilon_u.vector()[sphere_dofs]
    v.vector()[sphere_dofs]=ss_v.vector()[sphere_dofs]+epsilon_v.vector()[sphere_dofs]
#------------------------------------------------------------------
# Function 5: "VF_and_FEM_Schnakenberg_sphere_with_holes"
# The function calculates the matrices in the FEM formulation by
# defining the variational formulation and projecting it onto the
# space of continuous picewise linear bases functions over the mesh.
# The function takes the following inputs:
# 1. The list parameters containing the parameters of the Schnakenberg model,
# 2. u and v being the states of the Schnakenberg model,
# 3. phi_1 and phi_2 being the trial functions,
# 4. u_prev, v_prev being the solution at the previous time step,
# 5. dx_list containing the integration measures needed in order to define the variational formulation,
# 6. mesh containing the FEM-mesh of the sphere with a hole (potentially).
# The function returns the following output:
# 1. The mass_matrix in the LHS of the matrix equation resulting from FEM,
# 2.The stiffness_matrix in the LHS of the matrix equation resulting from FEM,
# 3.The reaction_matrix in the LHS of the matrix equation resulting from FEM,
# 4. The mass_load_vector_rhs in the RHS of the matrix equation resulting from FEM,
# 5.The stiffness_load_vector_rhs in the RHS of the matrix equation resulting from FEM,
# 6.The reaction_load_vector_rhs in the RHS of the matrix equation resulting from FEM.
#------------------------------------------------------------------

def VF_and_FEM_Schnakenberg_sphere_with_holes(parameters, u, v, phi_1, phi_2, u_prev, v_prev, dx_list, mesh):
    # Extract the parameters of the Schnakenberg model
    a = parameters[0]
    b = parameters[1]
    d = parameters[2]    
    gamma = parameters[3]
    # DEFINE THE REACTION TERMS IN THE SCHNAKENBERG MODEL IN A MIXED IMPLICIT EXPLICIT FASHION
    # (here we exclude the gamma factor as we include it later in the time stepping)
    # DISCRETISATION BASED ON THE METHOD CALLED 1-SBEM
    f = a - u + u*u_prev*v_prev
    g = b - u_prev*u_prev*v
    # DEFINE OUR THREE TYPE OF TERMS IN THE VARIATIONAL FORMULATION (VF):
    # NOTE THAT ALL TERMS SHOULD BE CONSIDERED AS BEING ON THE LHS, i.e., LHS=0
    # 1. Mass_form: originating from the time derivatives in the PDEs
    # 2. Stiffness form: originating from the Laplacian in the PDEs
    # 3. Reaction form: originating from the reaction terms in the PDEs
    #--------------------------------------------------------
    # THE REST OF THE SPHERE
    # Extract our first measure for the rest of the sphere
    dx = dx_list[0]
    # Initiate our three forms
    # Time evolution
    mass_form_u = (u-u_prev)*phi_1*dx
    mass_form_v = (v-v_prev)*phi_2*dx
    # Diffusion
    stiffness_form_u = dot(grad(u), grad(phi_1))*dx
    stiffness_form_v = d*dot(grad(v), grad(phi_2))*dx
    # Reactions
    reaction_form_u = -f*phi_1*dx
    reaction_form_v = -g*phi_2*dx
    # Return the output
    return mass_form_u, mass_form_v, stiffness_form_u, stiffness_form_v, reaction_form_u, reaction_form_v 

#------------------------------------------------------------------
# Function 6: "residual_Schnakenberg_sphere_with_holes"
#------------------------------------------------------------------
# The function calculates the residual forms being a measure of
# how much the states or concentration profiles change over the course
# of one time step. This is used in order to guide the choice of the
# next time step in the adaptive time stepping procedure where a large
# change in the concentration profile will lead to a small time step
# being chosen and vice versa a small change in the concentration
# profile will lead to a large time step being chosen.
# The input of the function is the following:
# 1. u and v being the states or concentration profiles at the current time step,
# 2. phi_1 and phi_2 being the test functions,
# 3. U_prev being the concentration profiles at the previous time step,
# 4. dx_list being the list of all integration measures.
# The function returns the residual forms.
#------------------------------------------------------------------

def residual_Schnakenberg_sphere_with_holes(parameters,phi_1,phi_2,u_prev,v_prev,u_curr,v_curr,dx_list,mesh,k):
    # Extract the parameters of the Schackenberg model
    a = parameters[0]
    b = parameters[1]
    d = parameters[2]    
    gamma = parameters[3]
    # DEFINE THE REACTION TERMS IN THE SCHNAKENBERG MODEL IN A MIXED IMPLICIT EXPLICIT FASHION
    # (here we exclude the gamma factor as we include it later in the time stepping)
    f = a - u_curr + u_curr*u_prev*v_prev
    g = b - u_prev*u_prev*v_curr
    # DEFINE OUR THREE TYPE OF TERMS IN THE VARIATIONAL FORMULATION (VF):
    # NOTE THAT ALL TERMS SHOULD BE CONSIDERED AS BEING ON THE LHS, i.e., LHS=0
    # 1. Mass_form: originating from the time derivatives in the PDEs
    # 2. Stiffness form: originating from the Laplacian in the PDEs
    # 3. Reaction form: originating from the reaction terms in the PDEs
    # THE REST OF THE SPHERE
    # Extract our first measure for the rest of the sphere
    dx = dx_list[0]
    # Initiate our three forms
    # Time evolution
    mass_form_u = (u_curr-u_prev)*phi_1*dx
    mass_form_v = (v_curr-v_prev)*phi_2*dx
    # Diffusion
    stiffness_form_u = dot(grad(u_curr), grad(phi_1))*dx
    stiffness_form_v = d*dot(grad(v_curr), grad(phi_2))*dx
    # Reactions
    reaction_form_u = -f*phi_1*dx
    reaction_form_v = -g*phi_2*dx
    # Lastly define the residual as well
    residual_form_u = mass_form_u + k*(stiffness_form_u + gamma*reaction_form_u)
    residual_form_v = mass_form_v + k*(stiffness_form_v + gamma*reaction_form_v)    
    # Return the residual form
    return residual_form_u + residual_form_v
#------------------------------------------------------------------
# Function 7: "compute_spectral_coefficients_nohole"
#------------------------------------------------------------------
def compute_spectral_coefficients_nohole(u, dx, hole_radius,folder_str):
    # Define the radius and the degree of the local approximation
    r = 1.0     # Radius of the sphere
    degree = 1  # Degree of local approximation
    # DEFINE THE BASIS FUNCTIONS FOR THE SPHERICAL HARMONICS
    # The real spherical harmonics in Cartesian coordinates 
    # n=0
    Y_00 = Expression("sqrt(1/(4*pi))", degree=degree)
    # n=1
    Y_10 = Expression("sqrt(3/(4*pi))*x[2]/r", r=r, degree=degree)
    Y_1p1 = Expression("sqrt(3/(4*pi))*x[0]/r", r=r, degree=degree)
    # n=2
    Y_20 = Expression("sqrt(5/(16*pi))*(3*pow(x[2], 2) - pow(r, 2))/pow(r, 2)", r=r, degree=degree)
    Y_2p1 = Expression("sqrt(15/(4*pi))*x[0]*x[2]/pow(r, 2)", r=r, degree=degree)
    Y_2p2 = Expression("sqrt(15/(16*pi))*(pow(x[0], 2) - pow(x[1], 2))/pow(r, 2)", r=r, degree=degree)
    # n=3
    Y_30 = Expression("sqrt(7/(16*pi))*(5*pow(x[2], 3) - 3*x[2]*pow(r, 2))/pow(r, 3)", r=r, degree=degree)
    Y_3p1 = Expression("sqrt(21/(32*pi))*x[0]*(5*pow(x[2], 2) - pow(r, 2))/pow(r, 3)", r=r, degree=degree)
    Y_3p2 = Expression("sqrt(105/(16*pi))*x[2]*(pow(x[0], 2) - pow(x[1], 2))/pow(r, 3)", r=r, degree=degree)
    Y_3p3 = Expression("sqrt(35/(16*pi))*x[0]*(pow(x[0], 2) - 3*pow(x[1], 2))/pow(r, 3)", r=r, degree=degree)
    # n=4
    Y_40 = Expression("sqrt(9/(256*pi))*(35*pow(x[2], 4) - 30*pow(x[2], 2)*pow(r, 2) + 3*pow(r, 4))/pow(r, 4)", r=r, degree=degree)
    Y_4p1 = Expression("sqrt(45/(32*pi))*x[0]*(7*pow(x[2], 3) - 3*x[2]*pow(r, 2))/pow(r, 4)", r=r, degree=degree)
    Y_4p2 = Expression("sqrt(45/(64*pi))*(pow(x[0], 2) - pow(x[1], 2))*(7*pow(x[2], 2) - pow(r, 2))/pow(r, 4)", r=r, degree=degree)
    Y_4p3 = Expression("sqrt(315/(32*pi))*x[0]*x[2]*(pow(x[0], 2) - 3*pow(x[1], 2))/pow(r, 4)", r=r, degree=degree)
    Y_4p4 = Expression("sqrt(315/(256*pi))*(pow(x[0], 2)*(pow(x[0], 2) - 3*pow(x[1], 2)) - pow(x[1], 2)*(3*pow(x[0], 2) - pow(x[1], 2)))/pow(r, 4)", r=r, degree=degree)
    # ASSEMBLE THE SOLUTION AS A FUNCTION OF THESE BASIS FUNCTIONS
    # The real part of the complex spectral coefficients of
    # the supplied function according to Chaplain,
    # i.e., <u, Re(Y_n^m)> = 1/sqrt(2)<u, Y_(n|m|)> for m non-zero
    # and   <u, Re(Y_n^0)> = <u, Y_(n0)> for m = 0, where
    # Y_n^m is a complex spherical harmonic, Y_(n|m|) is a real one
    # and the Condon-Shortley factor has been neglected  
    # n=0
    U_00 = assemble(u*Y_00*dx)
    # n=1
    U_10 = assemble(u*Y_10*dx)
    U_1p1 = 1/sqrt(2)*assemble(u*Y_1p1*dx)
    # n=2
    U_20 = assemble(u*Y_20*dx)
    U_2p1 = 1/sqrt(2)*assemble(u*Y_2p1*dx)
    U_2p2 = 1/sqrt(2)*assemble(u*Y_2p2*dx)
    # n=3
    U_30 = assemble(u*Y_30*dx)
    U_3p1 = 1/sqrt(2)*assemble(u*Y_3p1*dx)
    U_3p2 = 1/sqrt(2)*assemble(u*Y_3p2*dx)
    U_3p3 = 1/sqrt(2)*assemble(u*Y_3p3*dx)
    # n=4
    U_40 = assemble(u*Y_40*dx)
    U_4p1 = 1/sqrt(2)*assemble(u*Y_4p1*dx)
    U_4p2 = 1/sqrt(2)*assemble(u*Y_4p2*dx)
    U_4p3 = 1/sqrt(2)*assemble(u*Y_4p3*dx)
    U_4p4 = 1/sqrt(2)*assemble(u*Y_4p4*dx)
    # Save the results
    # Create a list of the coefficients of the basis functions
    coeff_list = [U_00, U_10, U_1p1, U_20, U_2p1, U_2p2, U_30, U_3p1, U_3p2, U_3p3, U_40, U_4p1, U_4p2, U_4p3, U_4p4]
    # Create a list of the corresponding hole radius
    hole_radii_list = [hole_radius, hole_radius, hole_radius, hole_radius, hole_radius, hole_radius, hole_radius, hole_radius, hole_radius, hole_radius, hole_radius, hole_radius, hole_radius, hole_radius, hole_radius]
    # Create an array with the temporary data
    temp_data = np.array([hole_radii_list,coeff_list])
    # Create a list of the column names
    col_names = ["Hole radius in the mesh", "Spectral coefficients of the eigenfunctions at final time t=T"]
    # Create a list of the names of the eigenfunctions
    row_names = ["$\\gamma_{0,0}$", "$\\gamma_{1,0}$", "$\\gamma_{1,1}$", "$\\gamma_{2,0}$", "$\\gamma_{2,1}$", "$\\gamma_{2,2}$", "$\\gamma_{3,0}$", "$\\gamma_{3,1}$", "$\\gamma_{3,2}$", "$\\gamma_{3,3}$", "$\\gamma_{4,0}$", "$\\gamma_{4,1}$", "$\\gamma_{4,2}$", "$\\gamma_{4,3}$", "$\\gamma_{4,4}$"]
    # Create a datafrane
    data = temp_data.T
    # Create a pandas data frame which we want to save
    df = pd.DataFrame(data,columns=col_names)
    # Add the row_name afterwards
    df.insert(0,'Eigenfunctions', row_names)
    # Finally, save the data frame to the given file name
    df.to_csv(folder_str + "/spectral_coefficients.csv")    
#------------------------------------------------------------------
# Function 8: "save_IC"
#------------------------------------------------------------------
# The function saves the initial condition for the current mesh and parameters. It takes the same input as
# the subsequent function, so read the comments for that function if you wish to know the name of this input.
def save_IC(num_holes,steady_states,numerical_parameters,radii_holes,ICs_around_steady_states):
    #--------------------------------------------------------------
    # STEP 1 OUT OF 7: EXTRACT PARAMETERS
    #--------------------------------------------------------------    
    # Extract the numerical parameters
    sigma = numerical_parameters[0] # Determining the perturbation in the initial conditions
    T = numerical_parameters[1] # Determining the end time for the FD time stepping scheme
    #--------------------------------------------------------------
    # STEP 2 OUT OF 7: DEFINE MESHES, INTEGRATION MEASURES AND
    # DEFINE THE HILBERT SPACE FOR THE FEM FORMULATION
    #--------------------------------------------------------------    
    # Read in the mesh depending on the number of holes on the sphere
    mesh, mvc_subdomains, mf_subdomains, dx_list = read_mesh_Schnakenberg_sphere_with_holes(num_holes,radii_holes)
    # Define the finite element space for the Schackenberg model using the mesh
    H = FunctionSpace(mesh, "P", 1)
    #--------------------------------------------------------------
    # STEP 3 OUT OF 7: DEFINE TEST FUNCTIONS AND TRIAL FUNCTIONS (I.E.
    # THE SOLUTION TO THE SCHNAKENBERG MODEL)
    #--------------------------------------------------------------    
    # Define the previous time step as a function of the function space H
    u_prev = Function(H) # Previous time step in the FD time stepping scheme
    v_prev = Function(H) # Previous time step in the FD time stepping scheme
    #--------------------------------------------------------------
    # STEP 4 OUT OF 7: SET THE INITIAL CONDITIONS OF THE SYSTEM,
    # CREATE AN OUTPUT FOLDER WHERE THE RESULTS ARE STORED
    # AND SAVE THE INITIAL CONDITIONS.
    #--------------------------------------------------------------    
    # Calculate the initial conditions
    initial_conditions_Schnakenberg_sphere_with_holes(H,mesh,mf_subdomains,num_holes,steady_states,sigma,u_prev,v_prev,ICs_around_steady_states)
    # Save the initial conditions to files
    File('../Output/fixed_IC_u.xml') << u_prev
    File('../Output/fixed_IC_v.xml') << v_prev
    # Save the mesh to a file as well
    File('../Output/fixed_IC_mesh.xml') << mesh
#------------------------------------------------------------------
# Function 9: "FEMFD_simulation_Schnakenberg_sphere_with_holes"
#------------------------------------------------------------------
# The functions solves the Schnakenberg RD model on the sphere with potential holes and potentially the parameters are altered in the regions adjacent to the holes. The function does not return any output but it writes the concentration profiles of u and v respectively to vtk files which are stored in an appropriately named sub folder of the folder named "../Output". The function takes the following inputs:
# 1. The parameter num_holes determining which mesh that is read as they are classified according to the number of holes that are added on the sphere,
# 2. The list parameters=[a,b,d,gamma] containing all the parameters of the Schnakenberg model,
# 3. The list steady_states=[u0,v0] containing the two states of the Schnakenberg model,
# 4. The list numerical_parameters=[sigma,T] where sigma determines the perturbation in the initial condition and T determines the end time for the FD time stepping scheme,
# 5. The list radii_holes containing the list of the radii of the holes in the mesh,
# 6. A logical variable called ICs_around_steady_states which determines whether the initial conditions are set to the steady states values or not,
# 7. A logical variable called load_IC which determines whether we generate new ICs or if load some fixed ICs,

def FEMFD_simulation_Schnakenberg_sphere_with_holes(num_holes,parameters,steady_states,numerical_parameters,radii_holes,ICs_around_steady_states,number_of_repititions,load_IC):
    #--------------------------------------------------------------
    # STEP 1 OUT OF 7: EXTRACT PARAMETERS
    #--------------------------------------------------------------    
    # Extract the parameters of the Schnakenberg model
    a = parameters[0]
    b = parameters[1]
    d = parameters[2]    
    gamma = parameters[3]
    # Re-define parameters as constants
    a = Constant(a)
    b = Constant(b)
    d = Constant(d)
    gamma = gamma
    gamma_const = Constant(gamma)
    # Save these in a new list
    parameters_as_constants = [a, b, d, gamma]
    # Extract the numerical parameters
    sigma = numerical_parameters[0] # Determining the perturbation in the initial conditions
    T = numerical_parameters[1] # Determining the end time for the FD time stepping scheme
    #--------------------------------------------------------------
    # STEP 2 OUT OF 7: DEFINE MESHES, INTEGRATION MEASURES AND
    # DEFINE THE HILBERT SPACE FOR THE FEM FORMULATION
    #--------------------------------------------------------------    
    # Read in the mesh depending on the number of holes on the sphere
    mesh, mvc_subdomains, mf_subdomains, dx_list = read_mesh_Schnakenberg_sphere_with_holes(num_holes,radii_holes)
    # Define the finite element space for the Schackenberg model using the mesh
    H = FunctionSpace(mesh, "P", 1)
    #--------------------------------------------------------------
    # STEP 3 OUT OF 7: DEFINE TEST FUNCTIONS AND TRIAL FUNCTIONS (I.E.
    # THE SOLUTION TO THE SCHNAKENBERG MODEL)
    #--------------------------------------------------------------    
    # Define test functions for the variational formulation (VF) 
    phi_1 = TestFunction(H)
    phi_2 = TestFunction(H)
    # Define the trial functions (i.e. solutions of the Schnakenberg model) for the VF.
    # Think of them as the "analytical solution" which we want to approximate.
    u = TrialFunction(H)
    v = TrialFunction(H)
    # Define the previous time step as a function of the function space H
    u_prev = Function(H) # Previous time step in the FD time stepping scheme
    u_curr = Function(H) # Current time step in the FD time stepping scheme
    v_prev = Function(H) # Previous time step in the FD time stepping scheme
    v_curr = Function(H) # Current time step in the FD time stepping scheme
    #--------------------------------------------------------------
    # STEP 4 OUT OF 7: SET THE INITIAL CONDITIONS OF THE SYSTEM,
    # CREATE AN OUTPUT FOLDER WHERE THE RESULTS ARE STORED
    # AND SAVE THE INITIAL CONDITIONS.
    #--------------------------------------------------------------    
    # We wish to have the most explanatory name of our output folder explaining the chosen parameters. So we create a series of strings which gives us all the information we need in the very name of the folder
    folder_str = "../Output/"
    hole_str = "h_" + str(num_holes) + "_"    
    radius_str = ""
    for radius in radii_holes:
        radius_str += "r_" + str(round(radius,3)).replace(".","p") + "_"
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
    # Gather all these substrings into one giant string where we will save the output files
    output_folder_str = folder_str + hole_str + radius_str + a_str + b_str + d_str + gamma_str + sigma_str + T_str + IC_str
    for repitition_index in range(number_of_repititions):
        # Define two output files based on this giant result folder where we have one output file for each of the two states
        vtkfile_u = File(output_folder_str+ "iteration_" + str(repitition_index) + "/" + "u.pvd")
        vtkfile_v = File(output_folder_str+"iteration_" + str(repitition_index) + "/" + "v.pvd")        
        # Set the time to zero as we are looking at the initial conditions
        t = 0.0
        if load_IC:
            # Load the old mesh
            mesh_old = Mesh('../Output/fixed_IC_mesh.xml')
            # Define a function space on the old mesh
            H_old = FunctionSpace(mesh_old, "P", 1)
            # Load the old initial conditions on the old function space
            u_old = Function(H_old, '../Output/fixed_IC_u.xml')
            v_old = Function(H_old, '../Output/fixed_IC_v.xml')
            # Project the old initial conditions onto the new function space
            u_prev = interpolate(u_old, H)
            #u_prev.set_allow_extrapolation(True)
            v_prev = interpolate(v_old, H)
            #v_prev.set_allow_extrapolation(True)            
        else:
            # Calculate the initial conditions
            initial_conditions_Schnakenberg_sphere_with_holes(H,mesh,mf_subdomains,num_holes,steady_states,sigma,u_prev,v_prev,ICs_around_steady_states)
        # Save the two initial conditions in the output folder
        u_prev.rename("Concentration profile, $u(\mathbf{x},t)$","u")
        vtkfile_u << (u_prev, t)
        v_prev.rename("Concentration profile, $v(\mathbf{x},t)$","v")    
        vtkfile_v << (v_prev, t)
        #--------------------------------------------------------------
        # STEP 5 OUT OF 7: DEFINE THE MATRICES RESULTING FROM THE
        # VF AND THE FEM
        #--------------------------------------------------------------
        # Compute the forms from the VF
        mass_form_u, mass_form_v, stiffness_form_u, stiffness_form_v, reaction_form_u, reaction_form_v = VF_and_FEM_Schnakenberg_sphere_with_holes(parameters, u, v, phi_1, phi_2, u_prev, v_prev, dx_list, mesh)
        # Assemble time-independent matrices for LHS
        mass_matrix_u = assemble(lhs(mass_form_u), keep_diagonal=True)
        mass_matrix_v = assemble(lhs(mass_form_v), keep_diagonal=True)
        stiffness_matrix_u = assemble(lhs(stiffness_form_u), keep_diagonal=True)
        stiffness_matrix_v = assemble(lhs(stiffness_form_v), keep_diagonal=True)
        # Assemble time-independent vectors for RHS
        reaction_vector_u = assemble(rhs(reaction_form_u))
        reaction_vector_v = assemble(rhs(reaction_form_v))
        # Compute time-dependent forms for LHS
        reaction_form_u_lhs = lhs(reaction_form_u)
        reaction_form_v_lhs = lhs(reaction_form_v)
        # Compute time-dependent forms for RHS
        mass_form_u_rhs = rhs(mass_form_u)
        mass_form_v_rhs = rhs(mass_form_v)
        #--------------------------------------------------------------
        # STEP 6 OUT OF 7: TIME STEPPING USING FD IN TIME AND FEM IN SPACE
        #--------------------------------------------------------------
        # Define the constant time step
        dt = 1e-2
        k = Constant(dt) # For the fem solver as well
        # Define an iterator for the time stepping keeping track of
        # how many iterations that has passed
        t_it = 0
        # Also, define the current time outside the loop, so that we can save
        # the final concentration profile after the looping is done.
        t = 0
        # Previous time step
        t_prev = 0
        # Save time every time with value 0.5
        save_iteration = 0.5
        #--------------------------------------------------------------
        # STEP 7 OUT OF 7: CALCULATE THE RESIDUAL FORMS NEEDED FOR THE
        #ADAPTIVE TIME STEPPING
        #--------------------------------------------------------------
        residual_form = residual_Schnakenberg_sphere_with_holes(parameters_as_constants, phi_1, phi_2, u_prev, v_prev, u_curr, v_curr, dx_list, mesh, k)
        #----------------------------------------------------------------------------------
        # We solve the time stepping adaptively until the end time is reached 
        while t < T:
            # Save the previous time step
            t_prev = t
            # Update current time and iteration number 
            t_it += 1
            t += dt
            k = Constant(dt)
            # Assemble time-dependent matrices for LHS
            reaction_matrix_u = assemble(reaction_form_u_lhs, keep_diagonal=True)
            reaction_matrix_v = assemble(reaction_form_v_lhs, keep_diagonal=True)
            # Assemble time-dependent vectors for RHS
            mass_vector_u = assemble(mass_form_u_rhs)
            mass_vector_v = assemble(mass_form_v_rhs)            
            # Assemble system matrices and rhs vector
            # SYSTEM WITH REACTIONS
            A_u = mass_matrix_u + k*(stiffness_matrix_u + gamma*reaction_matrix_u)
            A_v = mass_matrix_v + k*(stiffness_matrix_v + gamma*reaction_matrix_v)
            b_u = mass_vector_u + k*gamma*reaction_vector_u
            b_v = mass_vector_v + k*gamma*reaction_vector_v
            # Solve linear variational problems for time step
            solve(A_u, u_curr.vector(), b_u)
            solve(A_v, v_curr.vector(), b_v)
            # Save and check the solution (every whatever iteration)
            if  (t_prev < save_iteration) and (t > save_iteration):
                # Increase the iterations
                save_iteration += 0.5
                # Save the components in the data files
                u_curr.rename("Concentration profile, $u(\mathbf{x},t)$","u")
                vtkfile_u << (u_curr, t)
                v_curr.rename("Concentration profile, $v(\mathbf{x},t)$","v")
                vtkfile_v << (v_curr, t)
                # Iteration health check
                R = assemble(residual_form)
                l2_norm_R = norm(R, 'l2')
                print("\t\tIteration %d, t\t=\t%0.15f out of %0.3f"%(t_it,t,T))
                print("\t\tl2_norm_of_R = ", l2_norm_R)
                print("\t\tdt = ", dt)
                # In case we have negative concentrations, we break
                if u_curr.vector().min() < 0.0 or v_curr.vector().min() < 0.0:
                    print("\n\t\tINSTABILITY DETECTED! ### TERMINATE SIMULATION ###\n")
                    break
            # Update old solution
            u_prev.assign(u_curr)
            v_prev.assign(v_curr)
        # WE ALSO SAVE THE VERY LAST ITERATION WHEN ALL THE TIME STEPPING IS DONE.
        u_curr.rename("Concentration profile, $u(\mathbf{x},t)$","u")
        vtkfile_u << (u_curr, t)
        v_curr.rename("Concentration profile, $v(\mathbf{x},t)$","v")
        vtkfile_v << (v_curr, t)        
        print("\n\n\t\tALL IS FINE AND DANDY HERE!\n\n")
        print("Iterations are finished!")   
        # Compute and save the spectral coefficients as well
        if len(radii_holes)>0:
            compute_spectral_coefficients_nohole(u_curr, dx_list[0], radii_holes[0],output_folder_str+ "iteration_" + str(repitition_index) + "/")
        else:
            compute_spectral_coefficients_nohole(u_curr, dx_list[0], 0,output_folder_str+ "iteration_" + str(repitition_index) + "/")
#------------------------------------------------------------------
# Function 7: "project_eigen_functions_onto_mesh"
#------------------------------------------------------------------
def project_eigen_functions_onto_mesh(num_holes,radii_holes):
    # Read in the mesh depending on the number of holes on the sphere
    mesh, mvc_subdomains, mf_subdomains, dx_list = read_mesh_Schnakenberg_sphere_with_holes(num_holes,radii_holes)
    # Define the finite element space for the Schackenberg model using the mesh
    H = FunctionSpace(mesh, "P", 1)
    # Define the radius and the degree of the local approximation
    r = 1.0     # Radius of the sphere
    degree = 1  # Degree of local approximation
    # DEFINE THE OUTPUT STRING
    folder_str = "../Output/eigenfunctions/"
    hole_str = "h_" + str(num_holes)    
    if num_holes == 0:
        radius_str = "/"
    else:
        radius_str = "_"
        for radius in radii_holes:
            radius_str += "r_" + str(round(radius,3)).replace(".","p") + "_/"
    # Gather all these substrings into one giant string where we will save the output files
    output_folder_str = folder_str + hole_str + radius_str
    # DEFINE THE BASIS FUNCTIONS FOR THE SPHERICAL HARMONICS
    # The real spherical harmonics in Cartesian coordinates 
    # n=0
    Y_00 = Expression("sqrt(1/(4*pi))", degree=degree)
    P_Y_00 = project(Y_00, H)
    vtkfile_00 = File(output_folder_str+ "Y_00.pvd")
    P_Y_00.rename("$Y_{00}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_00")
    vtkfile_00 << (P_Y_00, 0)    
    # n=1
    Y_10 = Expression("sqrt(3/(4*pi))*x[2]/r", r=r, degree=degree)
    P_Y_10 = project(Y_10, H)
    vtkfile_10 = File(output_folder_str+ "Y_10.pvd")
    P_Y_10.rename("$Y_{10}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_10")
    vtkfile_10 << (P_Y_10, 0)        
    Y_1p1 = Expression("sqrt(3/(4*pi))*x[0]/r", r=r, degree=degree)
    P_Y_11 = project(Y_1p1, H)
    vtkfile_11 = File(output_folder_str+ "Y_11.pvd")
    P_Y_11.rename("$Y_{11}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_1p1")
    vtkfile_11 << (P_Y_11, 0)            
    # n=2
    Y_20 = Expression("sqrt(5/(16*pi))*(3*pow(x[2], 2) - pow(r, 2))/pow(r, 2)", r=r, degree=degree)
    P_Y_20 = project(Y_20, H)
    vtkfile_20 = File(output_folder_str+ "Y_20.pvd")
    P_Y_20.rename("$Y_{20}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_20")
    vtkfile_20 << (P_Y_20, 0)                
    Y_2p1 = Expression("sqrt(15/(4*pi))*x[0]*x[2]/pow(r, 2)", r=r, degree=degree)
    P_Y_21 = project(Y_2p1, H)
    vtkfile_21 = File(output_folder_str+ "Y_21.pvd")
    P_Y_21.rename("$Y_{21}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_2p1")
    vtkfile_21 << (P_Y_21, 0)                
    Y_2p2 = Expression("sqrt(15/(16*pi))*(pow(x[0], 2) - pow(x[1], 2))/pow(r, 2)", r=r, degree=degree)
    P_Y_22 = project(Y_2p2, H)
    vtkfile_22 = File(output_folder_str+ "Y_22.pvd")
    P_Y_22.rename("$Y_{22}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_2p2")
    vtkfile_22 << (P_Y_22, 0)                    
    # n=3
    Y_30 = Expression("sqrt(7/(16*pi))*(5*pow(x[2], 3) - 3*x[2]*pow(r, 2))/pow(r, 3)", r=r, degree=degree)
    P_Y_30 = project(Y_30, H)
    vtkfile_30 = File(output_folder_str+ "Y_30.pvd")
    P_Y_30.rename("$Y_{30}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_30")
    vtkfile_30 << (P_Y_30, 0)                
    Y_3p1 = Expression("sqrt(21/(32*pi))*x[0]*(5*pow(x[2], 2) - pow(r, 2))/pow(r, 3)", r=r, degree=degree)
    P_Y_31 = project(Y_3p1, H)
    vtkfile_31 = File(output_folder_str+ "Y_31.pvd")
    P_Y_31.rename("$Y_{31}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_3p1")
    vtkfile_31 << (P_Y_31, 0)                
    Y_3p2 = Expression("sqrt(105/(16*pi))*x[2]*(pow(x[0], 2) - pow(x[1], 2))/pow(r, 3)", r=r, degree=degree)
    P_Y_32 = project(Y_3p2, H)
    vtkfile_32 = File(output_folder_str+ "Y_32.pvd")
    P_Y_32.rename("$Y_{32}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_3p2")
    vtkfile_32 << (P_Y_32, 0)                    
    Y_3p3 = Expression("sqrt(35/(16*pi))*x[0]*(pow(x[0], 2) - 3*pow(x[1], 2))/pow(r, 3)", r=r, degree=degree)
    P_Y_33 = project(Y_3p3, H)
    vtkfile_33 = File(output_folder_str+ "Y_33.pvd")
    P_Y_33.rename("$Y_{33}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_3p3")
    vtkfile_33 << (P_Y_33, 0)                    
    # n=4
    Y_40 = Expression("sqrt(9/(256*pi))*(35*pow(x[2], 4) - 30*pow(x[2], 2)*pow(r, 2) + 3*pow(r, 4))/pow(r, 4)", r=r, degree=degree)
    P_Y_40 = project(Y_40, H)
    vtkfile_40 = File(output_folder_str+ "Y_40.pvd")
    P_Y_40.rename("$Y_{40}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_40")
    vtkfile_40 << (P_Y_40, 0)                    
    Y_4p1 = Expression("sqrt(45/(32*pi))*x[0]*(7*pow(x[2], 3) - 3*x[2]*pow(r, 2))/pow(r, 4)", r=r, degree=degree)
    P_Y_41 = project(Y_4p1, H)
    vtkfile_41 = File(output_folder_str+ "Y_41.pvd")
    P_Y_41.rename("$Y_{41}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_4p1")
    vtkfile_41 << (P_Y_41, 0)                    
    Y_4p2 = Expression("sqrt(45/(64*pi))*(pow(x[0], 2) - pow(x[1], 2))*(7*pow(x[2], 2) - pow(r, 2))/pow(r, 4)", r=r, degree=degree)
    P_Y_42 = project(Y_4p2, H)
    vtkfile_42 = File(output_folder_str+ "Y_42.pvd")
    P_Y_42.rename("$Y_{42}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_4p2")
    vtkfile_42 << (P_Y_42, 0)                    
    Y_4p3 = Expression("sqrt(315/(32*pi))*x[0]*x[2]*(pow(x[0], 2) - 3*pow(x[1], 2))/pow(r, 4)", r=r, degree=degree)
    P_Y_43 = project(Y_4p3, H)
    vtkfile_43 = File(output_folder_str+ "Y_43.pvd")
    P_Y_43.rename("$Y_{43}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_4p3")
    vtkfile_43 << (P_Y_43, 0)                    
    Y_4p4 = Expression("sqrt(315/(256*pi))*(pow(x[0], 2)*(pow(x[0], 2) - 3*pow(x[1], 2)) - pow(x[1], 2)*(3*pow(x[0], 2) - pow(x[1], 2)))/pow(r, 4)", r=r, degree=degree)    
    P_Y_44 = project(Y_4p4, H)
    vtkfile_44 = File(output_folder_str+ "Y_44.pvd")
    P_Y_44.rename("$Y_{44}(\mathbf{x}),\;\mathbf{x}\in S^2$","Y_4p4")
    vtkfile_44 << (P_Y_44, 0)                    
