# =================================================================================
# =================================================================================
# Script:"convert_mesh_from_msh_to_xdmf"
# Date: 2021-08-30
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
#import fenics
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
print("-----------------------------------------------------------")
print("\tCONVERT MSH MESHES TO XDMF")
print("-----------------------------------------------------------")
#------------------------------------------------------------
# NO HOLES
#------------------------------------------------------------
print("\t\tConvert mesh with zero holes")
# Read the mesh without holes using meshio
msh = meshio.read("../Meshes/s_h_0.msh")
# Create a triangle mesh using the function "create_mesh"
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
# Save the triangle mesh
meshio.write("../Meshes/s_h_0.xdmf", triangle_mesh)
print("\t\t\tDone!")
#------------------------------------------------------------
# 1 HOLE RADIUS 0.05
#------------------------------------------------------------
print("\t\tConvert mesh with 1 hole r=0.05")
# Read the mesh with one hole using meshio
msh = meshio.read("../Meshes/s_h_1_r_0p05.msh")
# Create a triangle mesh using the function "create_mesh"
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
# Save the triangle mesh
meshio.write("../Meshes/s_h_1_r_0p05.xdmf", triangle_mesh)
print("\t\t\tDone!")
#------------------------------------------------------------
# 1 HOLE RADIUS 0.10
#------------------------------------------------------------
print("\t\tConvert mesh with 1 hole r=0.10")
# Read the mesh with one hole using meshio
msh = meshio.read("../Meshes/s_h_1_r_0p10.msh")
# Create a triangle mesh using the function "create_mesh"
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
# Save the triangle mesh
meshio.write("../Meshes/s_h_1_r_0p10.xdmf", triangle_mesh)
print("\t\t\tDone!")
#------------------------------------------------------------
# 1 HOLE RADIUS 0.15
#------------------------------------------------------------
print("\t\tConvert mesh with 1 hole r=0.15")
# Read the mesh with one hole using meshio
msh = meshio.read("../Meshes/s_h_1_r_0p15.msh")
# Create a triangle mesh using the function "create_mesh"
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
# Save the triangle mesh
meshio.write("../Meshes/s_h_1_r_0p15.xdmf", triangle_mesh)
print("\t\t\tDone!")
#------------------------------------------------------------
# 1 HOLE RADIUS 0.20
#------------------------------------------------------------
print("\t\tConvert mesh with 1 hole r=0.20")
# Read the mesh with one hole using meshio
msh = meshio.read("../Meshes/s_h_1_r_0p20.msh")
# Create a triangle mesh using the function "create_mesh"
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
# Save the triangle mesh
meshio.write("../Meshes/s_h_1_r_0p20.xdmf", triangle_mesh)
print("\t\t\tDone!")
#------------------------------------------------------------
# 1 HOLE RADIUS 0.25
#------------------------------------------------------------
print("\t\tConvert mesh with 1 hole r=0.25")
# Read the mesh with one hole using meshio
msh = meshio.read("../Meshes/s_h_1_r_0p25.msh")
# Create a triangle mesh using the function "create_mesh"
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
# Save the triangle mesh
meshio.write("../Meshes/s_h_1_r_0p25.xdmf", triangle_mesh)
print("\t\t\tDone!")

