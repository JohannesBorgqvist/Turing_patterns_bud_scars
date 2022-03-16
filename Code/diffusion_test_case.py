# Most good things are contained in FEniCS
from fenics import *
# Import numpy as well
import numpy as np

def gradient_sphere(w,mesh):
    n = FacetNormal(mesh)
    #return n
    return grad(w) - dot(grad(w), n)*n
    #return ((dot(w,n) + abs(dot(w,n)))/(2.0))

def FEM_FD_simulation_pure_diffusion_sphere_with_holes(case_indicator):
    # Define the number of holes
    num_holes = 5
    #-------------------------------------------------------------------------
    #  READ MESH
    #-------------------------------------------------------------------------    
    # Allocate memory for the mesh and the mesh value collection
    mesh = Mesh()
    mvc_subdomains = MeshValueCollection("size_t", mesh, 2)
    # Define a mesh function taking the subdomains into account
    mf_subdomains = 0
    # Allocate memory for a list containing all the integration
    # measures involved in the variational formulation
    dx_list = []    
    # Define the string in which we read the mesh
    # depending on the number of holes
    if num_holes == 0 or num_holes==1 or num_holes==2 or num_holes==5: # No holes on the sphere
        mesh_str = "../Meshes/s_h_" + str(num_holes) +".xdmf"
    else:
        print("ERROR: No such mesh exist!")
        return mesh, mvc_subdomains, mf_subdomains, dx_list        
    # Read in the mesh and the subdomains into these two variables
    with XDMFFile(mesh_str) as infile:
        infile.read(mesh)
        infile.read(mvc_subdomains, "name_to_read")
    # The mesh value function keeping track of the subdomains    
    mf_subdomains = cpp.mesh.MeshFunctionSizet(mesh, mvc_subdomains)        
    #-------------------------------------------------------------------------
    #  DEFINE INTEGRATION MEASURES
    #-------------------------------------------------------------------------            
    # In the case of holes, three subdomains are marked in the mesh.
    # These are the holes with subdomain_id=3, the adjacent region
    # with subdomain_id=2 and the sphere with subdomain_id=1.
    dx_hole = Measure("dx", domain=mesh, subdomain_data=mf_subdomains, subdomain_id=3)
    dx_adjacent = Measure("dx", domain=mesh, subdomain_data=mf_subdomains, subdomain_id=2)
    dx_sphere = Measure("dx", domain=mesh, subdomain_data=mf_subdomains, subdomain_id=1)         
    #-------------------------------------------------------------------------
    #  DEFINE FUNCTION SPACE, TEST AND TRIAL FUNCTIONS
    #-------------------------------------------------------------------------            
    H = FunctionSpace(mesh, "CG", 1)
    phi   = TestFunction(H)
    u = TrialFunction(H)
    # Also, define a function over the mesh which will take the previous
    # value in the FD time stepping scheme
    U_prev = Function(H) # Previous time step in the FD time stepping scheme
    U_curr = Function(H) # Previous time step in the FD time stepping scheme      
    #-------------------------------------------------------------------------
    #  SET INITIAL CONDITIONS
    #-------------------------------------------------------------------------            
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
            if case_indicator == 1: # Simulations without a hole
                if mf_subdomains[cell] == 3: # The dofs in the hole
                    hole_dofs.extend(dofmap.cell_dofs(cell.index()))
                else: # The dofs on the rest of the sphere
                    sphere_dofs.extend(dofmap.cell_dofs(cell.index()))
            elif case_indicator == 2: # Simulations with a hole and a Neumann BC
                if mf_subdomains[cell] == 2: # The dofs in the adjacent region
                    hole_dofs.extend(dofmap.cell_dofs(cell.index()))
                elif mf_subdomains[cell] == 1: # The dofs on the rest of the sphere
                    sphere_dofs.extend(dofmap.cell_dofs(cell.index()))                
    # Find the unique dofs
    hole_dofs = list(set(hole_dofs))
    sphere_dofs = list(set(sphere_dofs))
    # Define an expression for the initial conditions which is just 1
    ic = Expression(('ic_u'),ic_u=1,degree=1)
    # Interpolate onto the function space
    ic_diffusion = interpolate(ic,H)
    # Assign the value 1 in the hole and zero everywhere else
    U_prev.vector()[sphere_dofs]=0*(ic_diffusion.vector()[sphere_dofs])
    U_prev.vector()[hole_dofs]=1*(ic_diffusion.vector()[hole_dofs])
    #-------------------------------------------------------------------------
    #  SAVE INITIAL CONDITION AND THE ORIGINAL MASS
    #-------------------------------------------------------------------------                
    # Save the previously defined initial condition in the solution
    #u = u_prev
    # Initiate the time to zero
    t = 0.0
    if case_indicator == 1:# Hole sphere no BC 
        # Define the residual and the mass of the system
        Residual_form = (U_curr-U_prev)*phi *dx_sphere + (U_curr-U_prev)*phi *dx_adjacent + (U_curr-U_prev)*phi *dx_hole
        mass = U_prev *dx_sphere + U_prev *dx_adjacent + U_prev *dx_hole
        # Create a vtk file where we can save the solution
        vtkfile_u = File("../Output/pure_diffusion_test_1/u.pvd")
    elif case_indicator == 2: # With a hole and Neumann BC
        # Define the residual and the mass of the system
        Residual_form = (U_curr-U_prev)*phi *dx_sphere + (U_curr-U_prev)*phi *dx_adjacent
        mass = U_prev *dx_sphere + U_prev *dx_adjacent
        # Create a vtk file where we can save the solution
        vtkfile_u = File("../Output/pure_diffusion_test_2/u.pvd")        



    # Save the initial condition in the Output folder
    U_prev.rename("Concentration profile, $u(\mathbf{x},t)$","u")
    vtkfile_u << (U_prev, t)
    # Allocate memory for a list of the time and the mass
    mass_list = []
    time_list = []
    # Save the initial mass and the time
    time_list.append(t)
    M = assemble(mass)
    mass_list.append(M)
    # Test to add a scaling
    d = 1
    #d = 5
    #d = 10
    #d = 15    
    #-------------------------------------------------------------------------
    #  DEFINE MATRIX EQUATION SYSTEM FROM FINITE ELEMENT METHOD (FEM)
    #-------------------------------------------------------------------------                    
    # DEFINE THE MATRICES WITH THE TIME DERIVATIVE (MASS FORM) AND DIFFUSION (STIFFNESS FORM)
    if case_indicator == 1:
        # Mass form
        mass_form = (u-U_prev)*phi *dx_sphere + (u-U_prev)*phi *dx_adjacent + (u-U_prev)*phi *dx_hole
        #mass_form = (U-u_prev)*phi *dx_sphere + (U-u_prev)*phi *dx_adjacent + (U-u_prev)*phi *dx_hole        
        # Stiffness form
        #stiffness_form =  dot(gradient_sphere(u,mesh), gradient_sphere(phi,mesh))*dx_hole + dot(gradient_sphere(u,mesh), gradient_sphere(phi,mesh))*dx_adjacent + dot(gradient_sphere(u,mesh), gradient_sphere(phi,mesh))*dx_sphere
        stiffness_form =  d*dot(grad(u), grad(phi))*dx_hole + d*dot(grad(u), grad(phi))*dx_adjacent + d*dot(grad(u), grad(phi))*dx_sphere
        #stiffness_form =  dot(grad(U), grad(phi))*dx_hole + dot(grad(U), grad(phi))*dx_adjacent + dot(grad(U), grad(phi))*dx_sphere
    elif case_indicator == 2:
        # Mass form
        mass_form = (u-U_prev)*phi *dx_sphere + (u-U_prev)*phi *dx_adjacent
        #mass_form = (U-u_prev)*phi *dx_sphere + (U-u_prev)*phi *dx_adjacent + (U-u_prev)*phi *dx_hole        
        # Stiffness form
        #stiffness_form =  dot(gradient_sphere(u,mesh), gradient_sphere(phi,mesh))*dx_hole + dot(gradient_sphere(u,mesh), gradient_sphere(phi,mesh))*dx_adjacent + dot(gradient_sphere(u,mesh), gradient_sphere(phi,mesh))*dx_sphere
        stiffness_form =  d*dot(grad(u), grad(phi))*dx_adjacent + d*dot(grad(u), grad(phi))*dx_sphere
        #stiffness_form =  dot(grad(U), grad(phi))*dx_hole + dot(grad(U), grad(phi))*dx_adjacent + dot(grad(U), grad(phi))*dx_sphere        
    # LHS: MATRICES
    # Mass matrix
    mass_matrix = assemble(lhs(mass_form), keep_diagonal = True)
    # Stiffness matrix
    stiffness_matrix = assemble(lhs(stiffness_form), keep_diagonal = True)
    # RHS: LOAD VECTORS
    mass_load_vector_rhs = rhs(mass_form) # Mass (Time derivative)
    stiffness_load_vector_rhs = rhs(stiffness_form) # Stiffness (Diffusion)
    #-------------------------------------------------------------------------
    #  SET UP TIME STEPPING USING FINITE DIFFERENCES (FD)
    #-------------------------------------------------------------------------                    
    # Define the adaptive time step
    dt = 1e-12 # we start with a really small step
    k = Constant(dt) # For the fem solver as well
    # Define an end time for the time stepping scheme
    T = 0.05
    # Define an iterator for the time stepping keeping track of
    # how many iterations that has passed
    t_it = 1
    # Define two tolerances for the adaptive time stepping
    TOL_tadapt = 1e-2 # Used for choosing the adaptive step
    dt_max = T/200 # An upper limit on the maximum time step
    # Update current time and iteration number
    t += dt
    k = Constant(dt)
    # CHECKING PROPERTIES OF THE RHS
    #print("Stiffness load, hey?")
    #print(norm(assemble(stiffness_load_vector_rhs)))
    print("Mass load, hey?")
    print(norm(assemble(mass_load_vector_rhs)))
    #----------------------------------------------------------------------------------
    # CHECKING PROPERTIES OF THE LHS
    print("Stiffness matrix, hey?")
    print("Det\t=\t%0.19f\n"%(np.linalg.det(stiffness_matrix.array())))
    print(stiffness_matrix.array()[0:10,0:10])    
    print("Mass matrix, hey?")
    print("Det\t=\t%0.19f\n"%(np.linalg.det(mass_matrix.array())))
    print(mass_matrix.array()[0:10,0:10])        
    #-------------------------------------------------------------------------
    #  SOLVE THE HEAT EQUATION ADAPTIVELY USING THE FD-FEM ALGORITHM
    #-------------------------------------------------------------------------                                
    while t < T:
        # Update current time and iteration number
        t_it += 1
        t += dt
        k = Constant(dt)
        # Assemble system matrix and  the load vectorvector
        A = mass_matrix + dt*stiffness_matrix
        b = assemble(mass_load_vector_rhs + k*(stiffness_load_vector_rhs))        #
        
        # Solve linear variational problem for time step
        solve(A,  U_curr.vector(), b)
        # Save the solution (every second iteration)
        if t_it %4 == 0:
            print("\t\tIteration %d, t\t=\t%0.15f out of %0.3f\tM(t)\t=\t%0.5f"%(t_it,t,T,M))
            U_curr.rename("Concentration profile, $u(\mathbf{x},t)$","u")
            vtkfile_u << (U_curr, t)
        elif t_it > 200:
            break
        # Compute next timestep adaptively with residual
        dt_old = dt
        R = assemble(Residual_form)
        dt = TOL_tadapt/norm(R, 'l2')
        dt = min(2.0*dt_old*dt/(dt_old + dt), dt_max)
        # Update old solution
        U_prev.assign(U_curr)
        # Save the initial mass and the time
        time_list.append(t)
        M = assemble(mass)
        mass_list.append(M)
    return time_list, mass_list
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
# =========================================================================================================================
# Actually doing the simulations 
# =========================================================================================================================
# CASE 1: DIFFUSION ON THE WHOLE SPHERE WITHOUT A HOLE
# Define the case we are looking at
case_indicator = 1
# Run the simulation on the sphere without a hole
t, M = FEM_FD_simulation_pure_diffusion_sphere_with_holes(case_indicator)
# Save the mass over time
# Define the string were the plot will be stored
file_str = "../Figures/FigS1/Input/mass_" + str(case_indicator) + ".tex"
# Define the string defining the settings for the plot
plot_str = "color=mass_1,line width=2pt,"
# Define the string with the legend
legend_str = "Sphere without a hole"
# Define the plot
plot_LaTeX_2D(t,M,file_str,plot_str,legend_str)
# CASE 1: DIFFUSION ON THE WHOLE SPHERE WITHOUT A HOLE
# Define the case we are looking at
case_indicator = 2
# Run the simulation on the sphere without a hole
t, M = FEM_FD_simulation_pure_diffusion_sphere_with_holes(case_indicator)
# Save the mass over time
# Define the string were the plot will be stored
file_str = "../Figures/FigS1/Input/mass_" + str(case_indicator) + ".tex"
# Define the string defining the settings for the plot
plot_str = "color=mass_2,line width=2pt,"
# Define the string with the legend
legend_str = "Sphere with five holes"
# Define the plot
plot_LaTeX_2D(t,M,file_str,plot_str,legend_str)
