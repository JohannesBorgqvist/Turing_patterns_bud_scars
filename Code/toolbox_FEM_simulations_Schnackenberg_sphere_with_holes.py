# =================================================================================
# =================================================================================
# Script:"toolbox_FEM_simulations_Schnackenberg_sphere_with_holes"
# Date: 2021-09-02
# Implemented by: Johannes Borgqvist
# Description:
# This is the main script containing all important functions needed in order to
# run the FEM simulations of the Schnackenberg model on the sphere with holes.
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
# Most good things are contained in FEniCS
from fenics import *
# Import numpy as well
import numpy as np
# =================================================================================
# =================================================================================
# Functions for conducting the FEM simulations
# =================================================================================
# =================================================================================
#------------------------------------------------------------------
# Function 1: "calculate_critical_parameters_Schnackenberg"
# The function takes the two parameters a and b of the
# Schnackenberg RD modelas well as the squared wavenumber
# k_squared as input. It returns the criticald diffusion
# parameter dc and the critical wavenumber gamma_c.
#------------------------------------------------------------------
def calculate_steady_states_and_critical_parameters_Schnackenberg(a,b,k_squared):
    # Calculate the steady states
    u_0 = a + b
    v_0 = ( (b) / ( (a + b)** 2))
    # Calculate the four partial derivatives of
    # the Jacobian matrix used in the linear
    # stability analysis
    f_u = ( (b - a) / (b + a) )
    f_v = ( (a + b)**2 )
    g_u = ( (-2*b) / (a+b) )
    g_v = -f_v
    # Calculate the critical diffusion value reported
    # in Chaplain 2001. Since the expression is quite
    # messy, we divide it into three parts.
    d_c = ( (f_u*g_v) - (2*f_v*g_u) )
    d_c +=  np.sqrt( d_c**2 - ((f_u**2)*(g_v**2)) )
    d_c = ( (d_c) / (f_u**2) )
    # Using the critical value of the diffusion d_c,
    # we also calculate the critical wave length number
    # gamma_c reported in Chaplain 2001.
    gamma_c = ( ( 2 * d_c * k_squared) / ( (d_c*f_u) + g_v ) )
    # Return the parameters 
    return u_0, v_0, d_c, gamma_c
#------------------------------------------------------------------
# Function 2: "read_mesh_Schnackenberg_sphere_with_holes"
# The function takes the number of wholes as inputs and
# returns the following four outputs:
# 1. The mesh "mesh",
# 2. A mesh value collection "mvc_subdomains" containing
# the subdomains,
# 3. A mesh_function "mf_subdomain",
# 4. A list of integration measures "dx_list" used in the
# variational formulation.
#------------------------------------------------------------------
def read_mesh_Schnackenberg_sphere_with_holes(num_holes):
    # Define the string in which we read the mesh
    # depending on the number of holes
    if num_holes == 0: # No holes on the sphere
        mesh_str = "../Meshes/sphere_with_no_holes_surfaces.xdmf"
    elif num_holes == 1:# One hole on the sphere
        mesh_str = "../Meshes/sphere_with_1_hole_surfaces.xdmf"
    # Allocate memory for the mesh and the mesh value collection
    mesh = Mesh()
    mvc_subdomains = MeshValueCollection("size_t", mesh, 2)
    # Read in the mesh and the subdomains into these two variables
    with XDMFFile(mesh_str) as infile:
        infile.read(mesh)
        infile.read(mvc_subdomains, "name_to_read")
    # Define a mesh function taking the subdomains into account
    mf_subdomains = cpp.mesh.MeshFunctionSizet(mesh, mvc_subdomains)
    # Allocate memory for a list containing all the integration
    # measures involved in the variational formulation
    dx_list = []
    # Add measures to this list depending on if we have holes
    # or not
    if num_holes == 0:
        # In this case, we only have one subdomain and thus only one
        # integration measure
        dx_sphere = Measure("dx", domain=mesh, subdomain_data=mf_subdomains, subdomain_id=1)            
        # Append this measure to the list
        dx_list.append(dx_sphere)
    else:
        # In the case of holes, three subdomains are marked in the mesh.
        # These are the holes with subdomain_id=2, the adjacent region
        # with subdomain_id=3 and the sphere with subdomain_id=4. The
        # subdomain_id=1 is reserved for the boundary curve around the
        # holes which is not of interest in this application as we
        # apply Neumann conditions around the holes. 
        dx_hole = Measure("dx", domain=mesh, subdomain_data=mf_subdomains, subdomain_id=2)
        dx_adjacent = Measure("dx", domain=mesh, subdomain_data=mf_subdomains, subdomain_id=3)
        dx_sphere = Measure("dx", domain=mesh, subdomain_data=mf_subdomains, subdomain_id=4)            
        # Append these measures to the list in the order sphere, adjacent and hole:
        dx_list.append(dx_sphere)
        dx_list.append(dx_adjacent)
        dx_list.append(dx_hole)
    # Lastly, return our mesh, the mesh value collection, the mesh
    # function and the integration measures.
    return mesh, mvc_subdomains, mf_subdomains, dx_list
#------------------------------------------------------------------
# Function 3: "set_up_FEM_Schnackenberg_sphere_with_holes"
# The function sets up the FEM framework for the Schnackenberg
# model on the sphere with holes. It takes a mesh as input
# and then it returns the following output:
# 1. The function space H,
# 2. The test function tuple called test_function_tuple storing the testfunctions phi_1 and phi_2,
# 3. The trial function tuple called trial_function_tuple storing the trial function u and v (i.e. the solution),
# 4. The function value at the previous time step in the FD time stepping stored in the variable u_prev,
#------------------------------------------------------------------
def define_Hilbert_space_Schnackenberg_sphere_with_holes(mesh):
    # Define the basis finite elements and their order. We pick linear bases functions. 
    # IMPORTANT WITH FEM: should be ''tetrahedon'' in three dimensions,     ''triangle'' in two and ''line'' in one
    #P1 = FiniteElement('P', triangle, 1)
    # For some reason there is an error message with using triangle here so we replaced it by "mesh.ufl_cell()" which seems to work
    P1 = FiniteElement('P', mesh.ufl_cell(), 1)    
    # Define a mixed element since we have a coupled system of two PDEs
    element = MixedElement([P1, P1])
    # Collect everything in the function space (name it H for a Hilbert space)
    H = FunctionSpace(mesh, element)
    # Return all of these elements
    return H
#------------------------------------------------------------------
# Function 4: "perturbations_ic__Schnackenberg_sphere_with_holes"

#------------------------------------------------------------------
def perturbations_ic__Schnackenberg_sphere_with_holes(epsilon,sigma):
    epsilon.vector().set_local(np.random.normal(0,sigma,epsilon.vector().local_size()))
    return epsilon
#------------------------------------------------------------------
# Function 5: "initial_conditions_Schnackenberg_sphere_with_holes"
# The function sets the initial conditions for the system and writes these two the function u_prev which is given as an input. The initial conditions for the two states in the Schnackenberg are set two a small perturbation around the steady state. The perturbation is given by a normal distribution with zero mean and where the variance is determined by the input sigma. Hence, no output is returned and the following input must be provided:
# 1. The Hilbert space H,
# 2. The mesh stored in the variable mesh,
# 3. The mesh function mf_subdomains indicating the subdomains in the mesh,
# 4. The number of holes on the sphere stored in num_holes,
# 5. The steady states of the Schnackenberg model stored in steady_states,
# 6. The variance of the perturbation determined by sigma,
# 7. The FEM function on the function space H stored in U in which we set the initial conditions. 
#------------------------------------------------------------------
def initial_conditions_Schnackenberg_sphere_with_holes(H,mesh,mf_subdomains,num_holes,steady_states,sigma,U):
    #----------------------------------------------------------------
    # DETAILED DESCRIPTION OF WHAT THE FUNCTION DOES
    # We assign the following initial conditions to the sphere:
    #1. u(t=0,*)=u0+epsilon
    #2. v(t=0,*)=v0+epsilon
    # where u0 and v0 are the provided steady states stored in the list steady_states and epsilon is a normally distributed variable that is drawn from the standard normal distribution with variance determined by the input sigma. If the mesh has a hole in it we set the initial condition u(t=0,x)=v(t=0,x)=0 for all spatial coordinates x in the hole. We use this solution as the hole is implemented as a sub-region in the mesh at hand.
    # In order to access the three different parts of the mesh (i.e the hole, adjacent region to the hole and the rest of the sphere) the so called dof_map of dolfin is used. This map enumerates all nodes in the mesh. Also, the provided mesh funcion in mf_subdomains keeps track of the various subdomains of the mesh which enable us to categorise all nodes as either "on the rest of the sphere" so to speak or "in the hole". Using the list of these nodes we are able to set the initial conditions properly by assigning values to the provided function U which is a vector valued function containing both states u and v in the Schnackenberg model.
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
            if mf_subdomains[cell] < 3: # The dofs in the hole
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
    ic = Expression(('ic_u','ic_v'),ic_u=u0, ic_v=v0,degree=1)
    # Interpolate onto the function space
    ss = interpolate(ic,H)
    # Define a function for the random perturbations
    epsilon = Function(H)
    # Fill this vector with a bunch of errors drawn
    # from a standard normal distribution with variance sigma
    epsilon = perturbations_ic__Schnackenberg_sphere_with_holes(epsilon,sigma)
    # Assign the value to our input vector U which will correspond to the initial conditions
    U.vector()[sphere_dofs]=ss.vector()[sphere_dofs]+epsilon.vector()[sphere_dofs]
    # Assign the value zero to the holes if a hole exists
    if len(hole_dofs)>0:
        U.vector()[hole_dofs]=0*(ss.vector()[hole_dofs]+epsilon.vector()[hole_dofs])
#------------------------------------------------------------------
# Function 7: "FEMFD_simulation_Schnackenberg_sphere_with_holes"
# 
#------------------------------------------------------------------
# The functions solves the Schnackenberg RD model on the sphere with potential holes and potentially the parameters are altered in the regions adjacent to the holes. The function does not return any output but it writes the concentration profiles of u and v respectively to vtk files which are stored in an appropriately named sub folder of the folder named "../Output". The function takes the following inputs:
# 1. The parameter num_holes determining which mesh that is read as they are classified according to the number of holes that are added on the sphere,
# 2. The list parameters=[a,b,d,gamma] containing all the parameters of the Schnackenberg model,
# 3. The list steady_states=[u0,v0] containing the two states of the Schnackenberg model,
# 4. The list numerical_parameters=[sigma,T] where sigma determines the perturbation in the initial condition and T determines the end time for the FD time stepping scheme,
# 5. The list activation parameters = [a_f,b_f]  which tells us how much the two rate parameters a and b are enhanced or weakened in the adjacent region around the hole,
# 6. A logical variable called cell_growth indicating whether the sphere increases in size over time or not.
def FEMFD_simulation_Schnackenberg_sphere_with_holes(num_holes,parameters,steady_states,numerical_parameters,activation_parameters,cell_growth):
    #--------------------------------------------------------------
    # STEP 1 OUT OF : EXTRACT PARAMETERS
    #--------------------------------------------------------------    
    # Extract the parameters of the Schackenberg model
    a = parameters[0]
    b = parameters[1]
    d = parameters[2]    
    gamma = parameters[3]
    # Extract the numerical parameters
    sigma = numerical_parameters[0] # Determining the perturbation in the initial conditions
    T = numerical_parameters[1] # Determining the end time for the FD time stepping scheme
    #--------------------------------------------------------------
    # STEP 2 OUT OF : DEFINE MESHES, INTEGRATION MEASURES AND
    # DEFINE THE HILBERT SPACE FOR THE FEM FORMULATION
    #--------------------------------------------------------------    
    # Read in the mesh depending on the number of holes on the sphere
    mesh, mvc_subdomains, mf_subdomains, dx_list = read_mesh_Schnackenberg_sphere_with_holes(num_holes)
    # Define the Hilbert space for the Schackenberg model on the mesh
    H = define_Hilbert_space_Schnackenberg_sphere_with_holes(mesh)    
    #--------------------------------------------------------------
    # STEP 3 OUT OF : DEFINE TEST FUNCTIONS AND TRIAL FUNCTIONS (I.E.
    # THE SOLUTION TO THE SCHNACKENBERG MODEL)
    #--------------------------------------------------------------    
    # Define test functions for the variational formulation (VF) 
    phi_1, phi_2 = TestFunctions(H)
    # Define the trial functions (i.e. solutions of the Schnackenberg model) for the VF
    u, v = TrialFunctions(H) # Current time step in the FD time stepping scheme
    # Define the previous time step as a function of the function space H
    u_prev = Function(H) # Previous time step in the FD time stepping scheme
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
    initial_conditions_Schnackenberg_sphere_with_holes(H,mesh,mf_subdomains,num_holes,steady_states,sigma,u_prev)
    # Split the initial conditions into its components part
    _u,_v = u_prev.split()
    # Save the two initial conditions in the output folder2
    vtkfile_u << (_u, t)
    vtkfile_v << (_v, t)
    print("\n\n\t\tALL IS FINE AND DANDY HERE!\n\n")
    print("The steady states are:")
    print(steady_states)
