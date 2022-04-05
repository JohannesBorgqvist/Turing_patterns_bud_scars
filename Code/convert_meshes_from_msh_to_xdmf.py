# =================================================================================
# =================================================================================
# Script:"convert_mesh_from_msh_to_xdmf"
# Date: 2022-03-30
# Implemented by: Johannes Borgqvist
# Description:
# The script reads in the spherical mesh that was generated using the pyhton
# version of gmsh into fenics using the meshio function.
# =================================================================================
# =================================================================================
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import meshio # To extract the important parts of the mesh
import numpy as np # For numerical calculations
# The main thing is dolfin allowing us to do the FEM calculations
from dolfin import *
# =================================================================================
# =================================================================================
# Functions
# =================================================================================
# =================================================================================
# Function 1: "create_mesh"
# The function was written by Jorgen Dokken, and what it does essentially is that
# it extracts the information from the msh mesh and then it returns a xdmf mesh
# which later can be used by dolfin and fenics. What is nice here that depending
# on the cell_type (i.e. "line", "triangle" or "tetra") a unique mesh with the
# lines or curves, surfaces and volumes is generated which can later be accessed
# by python. 
def create_mesh(mesh, cell_type, prune_z=False):
    # Retrieve 
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    out_mesh = meshio.Mesh(points=mesh.points, cells={
                           cell_type: cells}, cell_data={"name_to_read": [cell_data]})
    if prune_z:
        out_mesh.prune_z_0()
        
    return out_mesh

# =================================================================================
# =================================================================================
# Reading in the mesh and studying its properties
# =================================================================================
# =================================================================================
print("---------------------------------------------------------------------------------------------------")
print("\tCONVERT MSH MESHES TO XDMF")
print("---------------------------------------------------------------------------------------------------")
# Loop over all radii and convert the meshes
for hole_radius in np.arange(0,0.75,0.05):
    # Define the mesh name
    if hole_radius == 0:
        print("\tConvert mesh with zero holes")
        mesh_name = "s_h_0"
    else:
        print("\tConvert mesh with one hole and radius:\t\tr\t=\t%0.2f"%(hole_radius))        
        mesh_name = "s_h_1_r_"+str(round(hole_radius,3)).replace(".","p")
    # Read the msh file
    msh = meshio.read("../Meshes/" + mesh_name + ".msh")
    # Create a triangle mesh using the function "create_mesh"
    triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
    # Save the triangle mesh
    meshio.write("../Meshes/" + mesh_name + ".xdmf", triangle_mesh)
    print("\t\tDone!")    
print("---------------------------------------------------------------------------------------------------")
