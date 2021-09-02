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
