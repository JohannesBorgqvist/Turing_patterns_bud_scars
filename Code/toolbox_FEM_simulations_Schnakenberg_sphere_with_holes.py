# =================================================================================
# =================================================================================
# Script:"toolbox_FEM_simulations_Schnakenberg_sphere_with_holes"
# Date: 2021-12-17
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
# Import numpy as well
import numpy as np
# Write less trick
np.fac = np.math.factorial
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
            mesh_str += "r_" + str(round(radius,3)).replace(".","p") + "_"
    # Now, we add the file format to the string as well
    mesh_str += ".xdmf"
    # And tidy the file name up in necessary
    mesh_str.replace("_.xdmf",".xdmf")
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
# 7. The FEM function on the function space H stored in U in which we set the initial conditions. 
#------------------------------------------------------------------

def initial_conditions_Schnakenberg_sphere_with_holes(H,mesh,mf_subdomains,num_holes,steady_states,sigma,u,v):
    #----------------------------------------------------------------
    # DETAILED DESCRIPTION OF WHAT THE FUNCTION DOES
    # We assign the following initial conditions to the sphere:
    #1. u(t=0,*)=u0+epsilon
    #2. v(t=0,*)=v0+epsilon
    # where u0 and v0 are the provided steady states stored in the list steady_states and epsilon is a normally distributed variable that is drawn from the standard normal distribution with variance determined by the input sigma. If the mesh has a hole in it we set the initial condition u(t=0,x)=v(t=0,x)=0 for all spatial coordinates x in the hole. We use this solution as the hole is implemented as a sub-region in the mesh at hand.
    # In order to access the three different parts of the mesh (i.e the hole, adjacent region to the hole and the rest of the sphere) the so called dof_map of dolfin is used. This map enumerates all nodes in the mesh. Also, the provided mesh funcion in mf_subdomains keeps track of the various subdomains of the mesh which enable us to categorise all nodes as either "on the rest of the sphere" so to speak or "in the hole". Using the list of these nodes we are able to set the initial conditions properly by assigning values to the provided function U which is a vector valued function containing both states u and v in the Schnakenberg model.
    #----------------------------------------------------------------

    # Create a dofmap
    dofmap = H.dofmap()

    # Classify the dofs as either in the hole or on the rest of the sphere
    hole_dofs = []
    sphere_dofs = []

    # Loop over the cells in the mesh and add the dofs in the respective list
    for cell in cells(mesh):
        # Two cases: meshes without holes and with holes
        if num_holes == 0: # No holes
            # Add all dofs to the sphere
            sphere_dofs.extend(dofmap.cell_dofs(cell.index()))            
        else: # Mesh contains holes
            if mf_subdomains[cell] == 3: # The dofs in the hole
                hole_dofs.extend(dofmap.cell_dofs(cell.index()))
            else: # The dofs on the rest of the sphere
                sphere_dofs.extend(dofmap.cell_dofs(cell.index()))

    # Find the unique dofs
    hole_dofs = list(set(hole_dofs))
    sphere_dofs = list(set(sphere_dofs))

    # Access the component steady states as well
    u0 = steady_states[0]
    v0 = steady_states[1]    

    # Define an expression for the initial conditions based on the steady states
    ic_u = Expression("ic_u",ic_u=u0,degree=1)
    ic_v = Expression("ic_v",ic_v=v0,degree=1)

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

    # Assign the value zero to the holes if a hole exists
    if len(hole_dofs)>0:
        u.vector()[hole_dofs]=0*(ss_u.vector()[hole_dofs]+epsilon_u.vector()[hole_dofs])
        v.vector()[hole_dofs]=0*(ss_v.vector()[hole_dofs]+epsilon_v.vector()[hole_dofs])
  

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
# 6. mesh containing the FEM-mesh of the sphere with a hole (potentially),
# 7. activation_parameters containing the factors with which the reaction terms are increased or decreased in the region adjacent to the hole.
# The function returns the following output:
# 1. The mass_matrix in the LHS of the matrix equation resulting from FEM,
# 2.The stiffness_matrix in the LHS of the matrix equation resulting from FEM,
# 3.The reaction_matrix in the LHS of the matrix equation resulting from FEM,
# 4. The mass_load_vector_rhs in the RHS of the matrix equation resulting from FEM,
# 5.The stiffness_load_vector_rhs in the RHS of the matrix equation resulting from FEM,
# 6.The reaction_load_vector_rhs in the RHS of the matrix equation resulting from FEM.
#------------------------------------------------------------------

def VF_and_FEM_Schnakenberg_sphere_with_holes(parameters, u, v, phi_1, phi_2, u_prev, v_prev, dx_list, mesh, activation_parameters):

    # RE-DEFINE THE PARAMETERS AS CONSTANTS
    # Extract the parameters of the Schnakenberg model
    a = parameters[0]
    b = parameters[1]
    d = parameters[2]    
    gamma = parameters[3]

    # SPLIT THE SOLUTION AT THE PREVIOUS TIME STEP INTO ITS COMPONENT PARTS
    #u_prev,v_prev = U_prev.split()

    # SPLIT THE SOLUTION AT THE CURRENT TIME STEP INTO ITS COMPONENT PARTS
    #u_curr,v_curr = U_curr.split()

    # THE BELOW DOESNT WORK: Integral of type cell cannot contain a ReferenceNormal.
    #n = FacetNormal(mesh) 
    #def grad_T(w):
    #    return grad(w) - dot(grad(w), n)*n

    # DEFINE THE REACTION TERMS IN THE SCHNAKENBERG MODEL IN A MIXED IMPLICIT EXPLICIT FASHION
    # (here we exclude the gamma factor as we include it later in the time stepping)
    #-----------------------------------------------------
    # DISCRETISATION 1 (Heuristics from cdc42-project)
    #f = a - u + u_prev*u_prev*v_prev
    #g = b - u_prev*u_prev*v_prev
    #-----------------------------------------------------
    # DISCRETISATION 2 (1-SBEM)
    f = a - u + u*u_prev*v_prev
    g = b - u_prev*u_prev*v
    #-----------------------------------------------------
    # DISCRETISATION 3 (Fully explicit reaction (should be crap compared to others))
    #f = a - u_prev + u_prev*u_prev*v_prev
    #g = b - u_prev*u_prev*v_prev    
    #-----------------------------------------------------
    # DISCRETISATION 1 (Heuristics 2 from cdc42-project, maxed mixed)
    #f = a - u + u_prev*u_prev*v
    #g = b - u*u_prev*v_prev

    # Define the modified reaction terms if we have a hole
    if len(dx_list)>1:
        f_adjacent = a*activation_parameters[0] - u + u*u_prev*v_prev
        g_adjacent = b*activation_parameters[1] - u_prev*u_prev*v

    # DEFINE OUR THREE TYPE OF TERMS IN THE VARIATIONAL FORMULATION (VF):
    # NOTE THAT ALL TERMS SHOULD BE CONSIDERED AS BEING ON THE LHS, i.e., LHS=0
    # 1. Mass_form: originating from the time derivatives in the PDEs
    # 2. Stiffness form: originating from the Laplacian in the PDEs
    # 3. Reaction form: originating from the reaction terms in the PDEs
    # So we will define these terms by looping over our beloved measures
    # in the dx_list

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
    #stiffness_form_u = dot(grad_T(u), grad_T(phi_1))*dx
    #stiffness_form_v = d*dot(grad_T(v), grad_T(phi_2))*dx
    # Reactions
    reaction_form_u = -f*phi_1*dx
    reaction_form_v = -g*phi_2*dx

    # THE ADJACENT REGION
    if len(dx_list)>1:
        # Extract our second measure for the adjacent region
        dx = dx_list[1]
        # Time evolution
        mass_form_u += (u-u_prev)*phi_1*dx
        mass_form_v += (v-v_prev)*phi_2*dx
        # Diffusion
        stiffness_form_u += dot(grad(u), grad(phi_1))*dx
        stiffness_form_v += d*dot(grad(v), grad(phi_2))*dx
        # Reaction
        reaction_form_u += -f_adjacent*phi_1*dx
        reaction_form_v += -g_adjacent*phi_2*dx
    
    return mass_form_u, mass_form_v, stiffness_form_u, stiffness_form_v, reaction_form_u, reaction_form_v 

#------------------------------------------------------------------
# Function 6: "residual_Schnakenberg_sphere_with_holes"
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

def residual_Schnakenberg_sphere_with_holes(parameters,phi_1,phi_2,u_prev,v_prev,u_curr,v_curr,dx_list,mesh,activation_parameters,d,k):

    # RE-DEFINE THE PARAMETERS AS CONSTANTS

    # Extract the parameters of the Schackenberg model
    a = parameters[0]
    b = parameters[1]
    d = parameters[2]    
    gamma = parameters[3]

    # SPLIT THE SOLUTION AT THE PREVIOUS TIME STEP INTO ITS COMPONENT PARTS
    #u_prev,v_prev = U_prev.split()
    # SPLIT THE SOLUTION AT THE CURRENT TIME STEP INTO ITS COMPONENT PARTS
    #u_curr,v_curr = U_curr.split()

    # THE BELOW DOESNT WORK: Integral of type cell cannot contain a ReferenceNormal.
    #n = FacetNormal(mesh) 
    #def grad_T(w):
    #    return grad(w) - dot(grad(w), n)*n

    # DEFINE THE REACTION TERMS IN THE SCHNAKENBERG MODEL IN A MIXED IMPLICIT EXPLICIT FASHION
    # (here we exclude the gamma factor as we include it later in the time stepping)
    f = a - u_curr + u_curr*u_prev*v_prev
    g = b - u_prev*u_prev*v_curr

    # Define the modified reaction terms if we have a hole
    if len(dx_list)>1:
       f_adjacent = a*activation_parameters[0] - u_curr + u_curr*u_prev*v_prev
       g_adjacent = b*activation_parameters[1] - u_prev*u_prev*v_curr       

    # DEFINE OUR THREE TYPE OF TERMS IN THE VARIATIONAL FORMULATION (VF):
    # NOTE THAT ALL TERMS SHOULD BE CONSIDERED AS BEING ON THE LHS, i.e., LHS=0
    # 1. Mass_form: originating from the time derivatives in the PDEs
    # 2. Stiffness form: originating from the Laplacian in the PDEs
    # 3. Reaction form: originating from the reaction terms in the PDEs
    # So we will define these terms by looping over our beloved measures
    # in the dx_list

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
    #stiffness_form_u = dot(grad_T(u_curr), grad_T(phi_1))*dx
    #stiffness_form_v = d*dot(grad_T(v_curr), grad_T(phi_2))*dx
    # Reactions
    reaction_form_u = -f*phi_1*dx
    reaction_form_v = -g*phi_2*dx
    
    # THE ADJACENT REGION
    if len(dx_list)>1:
        # Extract our second measure for the adjacent region
        dx = dx_list[1]
        # Time evolution
        mass_form_u = (u_curr-u_prev)*phi_1*dx
        mass_form_v = (v_curr-v_prev)*phi_2*dx
        # Diffusion
        stiffness_form_u = dot(grad(u_curr), grad(phi_1))*dx
        stiffness_form_v = d*dot(grad(v_curr), grad(phi_2))*dx
        # Reactions
        reaction_form_u = -f_adjacent*phi_1*dx
        reaction_form_v = -g_adjacent*phi_2*dx

    # Lastly define the residual as well
    residual_form_u = mass_form_u + k*(stiffness_form_u + gamma*reaction_form_u)
    residual_form_v = mass_form_v + k*(stiffness_form_v + gamma*reaction_form_v)
    
    # Return the residual form
    return residual_form_u + residual_form_v

#------------------------------------------------------------------
# Function 7: "FEMFD_simulation_Schnakenberg_sphere_with_holes"
# 
#------------------------------------------------------------------
# The functions solves the Schnakenberg RD model on the sphere with potential holes and potentially the parameters are altered in the regions adjacent to the holes. The function does not return any output but it writes the concentration profiles of u and v respectively to vtk files which are stored in an appropriately named sub folder of the folder named "../Output". The function takes the following inputs:
# 1. The parameter num_holes determining which mesh that is read as they are classified according to the number of holes that are added on the sphere,
# 2. The list parameters=[a,b,d,gamma] containing all the parameters of the Schnakenberg model,
# 3. The list steady_states=[u0,v0] containing the two states of the Schnakenberg model,
# 4. The list numerical_parameters=[sigma,T] where sigma determines the perturbation in the initial condition and T determines the end time for the FD time stepping scheme,
# 5. The list activation parameters = [a_f,b_f]  which tells us how much the two rate parameters a and b are enhanced or weakened in the adjacent region around the hole,
# 6. A logical variable called cell_growth indicating whether the sphere increases in size over time or not.

def FEMFD_simulation_Schnakenberg_sphere_with_holes(num_holes,parameters,steady_states,numerical_parameters,activation_parameters,cell_growth):

    #--------------------------------------------------------------
    # STEP 1 OUT OF : EXTRACT PARAMETERS
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

    # We define a growth factor for gamma if needed
    gamma_growth = 0.05

    # Extract the numerical parameters
    sigma = numerical_parameters[0] # Determining the perturbation in the initial conditions
    T = numerical_parameters[1] # Determining the end time for the FD time stepping scheme

    #--------------------------------------------------------------
    # STEP 2 OUT OF : DEFINE MESHES, INTEGRATION MEASURES AND
    # DEFINE THE HILBERT SPACE FOR THE FEM FORMULATION
    #--------------------------------------------------------------    

    # Read in the mesh depending on the number of holes on the sphere
    mesh, mvc_subdomains, mf_subdomains, dx_list = read_mesh_Schnakenberg_sphere_with_holes(num_holes)

    # Define the finite element space for the Schackenberg model using the mesh
    #H = define_Hilbert_space_Schnakenberg_sphere_with_holes(mesh)    
    H = FunctionSpace(mesh, "P", 1)
    
    #--------------------------------------------------------------
    # STEP 3 OUT OF : DEFINE TEST FUNCTIONS AND TRIAL FUNCTIONS (I.E.
    # THE SOLUTION TO THE SCHNAKENBERG MODEL)
    #--------------------------------------------------------------    

    # Define test functions for the variational formulation (VF) 
    #phi_1, phi_2 = TestFunctions(H)
    phi_1 = TestFunction(H)
    phi_2 = TestFunction(H)

    # Define the trial functions (i.e. solutions of the Schnakenberg model) for the VF.
    # Think of them as the "analytical solution" which we want to approximate.
    #u, v = TrialFunctions(H)
    u = TrialFunction(H)
    v = TrialFunction(H)

    # Define the previous time step as a function of the function space H
    #U_prev = Function(H) # Previous time step in the FD time stepping scheme
    #U_curr = Function(H) # Current time step in the FD time stepping scheme
    u_prev = Function(H) # Previous time step in the FD time stepping scheme
    u_curr = Function(H) # Current time step in the FD time stepping scheme
    v_prev = Function(H) # Previous time step in the FD time stepping scheme
    v_curr = Function(H) # Current time step in the FD time stepping scheme

    #--------------------------------------------------------------
    # STEP 4 OUT OF : SET THE INITIAL CONDITIONS OF THE SYSTEM,
    # CREATE AN OUTPUT FOLDER WHERE THE RESULTS ARE STORED
    # AND SAVE THE INITIAL CONDITIONS.
    #--------------------------------------------------------------    

    # We wish to have the most explanatory name of our output folder explaining the chosen parameters. So we create a series of strings which gives us all the information we need in the very name of the folder.
    folder_str = "../Output/"
    hole_str = "h_" + str(num_holes) + "_"
    a_str = "a_" + str(round(a,2)).replace(".","p") + "_"
    b_str = "b_" + str(round(b,2)).replace(".","p") + "_"
    d_str = "d_" + str(round(d,2)).replace(".","p") + "_"
    gamma_str = "gamma_" + str(round(gamma,2)).replace(".","p") + "_"
    sigma_str = "sigma_" + str(round(sigma,2)).replace(".","p") + "_"
    T_str = "T_" + str(round(T,2)).replace(".","p") + "_"
    activation_str_a = "laa_" + str(round(activation_parameters[0],2)).replace(".","p") + "_"
    activation_str_b = "lab_" + str(round(activation_parameters[1],2)).replace(".","p") + "_"
    if cell_growth:
        cell_growth_str = "cg_yes/"
    else:
        cell_growth_str = "cg_no/"        

    # Gather all these substrings into one giant string where we will save the output files
    output_folder_str = folder_str + hole_str + a_str + b_str + d_str + gamma_str + sigma_str + T_str + activation_str_a + activation_str_b + cell_growth_str

    # Define two output files based on this giant result folder where we have one output file for each of the two states
    vtkfile_u = File(output_folder_str+"u.pvd")
    vtkfile_v = File(output_folder_str+"v.pvd")        

    # Set the time to zero as we are looking at the initial conditions
    t = 0.0

    # Calculate the initial conditions
    initial_conditions_Schnakenberg_sphere_with_holes(H,mesh,mf_subdomains,num_holes,steady_states,sigma,u_prev,v_prev)
    ## Split the initial conditions into its components part
    #u_prev,v_prev = U_prev.split()
    # Save the two initial conditions in the output folder
    u_prev.rename("Concentration profile, $u(\mathbf{x},t)$","u")
    vtkfile_u << (u_prev, t)
    v_prev.rename("Concentration profile, $v(\mathbf{x},t)$","v")    
    vtkfile_v << (v_prev, t)

    #--------------------------------------------------------------
    # STEP 5 OUT OF : DEFINE THE MATRICES RESULTING FROM THE
    # VF AND THE FEM
    #--------------------------------------------------------------

    # Compute the forms from the VF
    mass_form_u, mass_form_v, stiffness_form_u, stiffness_form_v, reaction_form_u, reaction_form_v = VF_and_FEM_Schnakenberg_sphere_with_holes(parameters, u, v, phi_1, phi_2, u_prev, v_prev, dx_list, mesh, activation_parameters)

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
    # STEP 7 OUT OF : THE ADAPTIVE TIME STEPPING USING FD IN TIME
    # AND FEM IN SPACE
    #--------------------------------------------------------------

    # Define the adaptive time step
    #dt = 1e-12 # we start with a really small step
    #dt = 1e0 # Doesn't work
    #dt = 1e-1 # Doesn't work
    dt = 1e-2
    
    k = Constant(dt) # For the fem solver as well
    
    # Define an iterator for the time stepping keeping track of
    # how many iterations that has passed
    t_it = 0
    
    # Define two tolerances for the adaptive time stepping
    #TOL_tadapt = 1e-5 # Used for choosing the adaptive step
    #dt_max = 1e-5 # An upper limit on the maximum time step
    
    #--------------------------------------------------------------
    # STEP 6 OUT OF : CALCULATE THE RESIDUAL FORMS NEEDED FOR THE
    #ADAPTIVE TIME STEPPING
    #--------------------------------------------------------------

    residual_form = residual_Schnakenberg_sphere_with_holes(parameters_as_constants, phi_1, phi_2, u_prev, v_prev, u_curr, v_curr, dx_list, mesh, activation_parameters, d, k)

    #----------------------------------------------------------------------------------

    # We solve the time stepping adaptively until the end time is reached 
    while t < T:

        # Update current time and iteration number 
        t_it += 1
        t += dt
        k = Constant(dt)

        #print("\t\tIteration %d, t\t=\t%0.15f out of %0.3f"%(t_it,t,T))

        # Check if we have cell growth or not. In this case we update gamma as well
        if cell_growth:
            gamma_time_dependent = gamma + gamma_growth*t
            gamma_time_dependent_const = gamma_const + Constant(gamma_growth*t)
        else:
            gamma_time_dependent = gamma
            gamma_time_dependent_const = gamma_const                    

        # Assemble time-dependent matrices for LHS
        reaction_matrix_u = assemble(reaction_form_u_lhs, keep_diagonal=True)
        reaction_matrix_v = assemble(reaction_form_v_lhs, keep_diagonal=True)

        # Assemble time-dependent vectors for RHS
        mass_vector_u = assemble(mass_form_u_rhs)
        mass_vector_v = assemble(mass_form_v_rhs)
            
        # Assemble system matrices and rhs vector
        # SYSTEM WITH REACTIONS
        A_u = mass_matrix_u + k*(stiffness_matrix_u + gamma_time_dependent_const*reaction_matrix_u)
        A_v = mass_matrix_v + k*(stiffness_matrix_v + gamma_time_dependent_const*reaction_matrix_v)
        b_u = mass_vector_u + k*gamma_time_dependent_const*reaction_vector_u
        b_v = mass_vector_v + k*gamma_time_dependent_const*reaction_vector_v

        # PURE DIFFUSION
        #A_u = mass_matrix_u + k*stiffness_matrix_u
        #A_v = mass_matrix_v + k*stiffness_matrix_v
        #b_u = mass_vector_u
        #b_v = mass_vector_v 

        # Solve linear variational problems for time step
        solve(A_u, u_curr.vector(), b_u)
        solve(A_v, v_curr.vector(), b_v)
        
        # Compute next timestep adaptively with residual
        #dt_old = dt
        #R = assemble(residual_form)
        #l2_norm_R = norm(R, 'l2')
        #dt = TOL_tadapt/l2_norm_R
        #dt = min(2.0*dt_old*dt/(dt_old + dt), dt_max)

        # Save and check the solution (every whatever iteration)
        if t_it % 10 == 0:

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
            if u_curr.vector().min() < 0.0 or v_curr.vector().min() < 0.0:
                print("\n\t\tINSTABILITY DETECTED! ### TERMINATE SIMULATION ###\n")
                break

        #if t_it > 2: #1000000: # Just in case we never finish, we break
        #    break
            
        # Update old solution
        u_prev.assign(u_curr)
        v_prev.assign(v_curr)
        
    print("\n\n\t\tALL IS FINE AND DANDY HERE!\n\n")
    print("Iterations are finished!")
