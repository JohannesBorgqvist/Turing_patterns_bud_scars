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
# Surpress all the output from FEniCS
set_log_level(LogLevel.ERROR)
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
# 1 HOLE RADIUS 0.3
#------------------------------------------------------------
print("\t\tConvert mesh with 1 hole r=0.3")
# Read the mesh with one hole using meshio
msh = meshio.read("../Meshes/s_h_1_r_0p3.msh")
# Create a triangle mesh using the function "create_mesh"
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
# Save the triangle mesh
meshio.write("../Meshes/s_h_1_r_0p3.xdmf", triangle_mesh)
print("\t\t\tDone!")
#------------------------------------------------------------
# 1 HOLE RADIUS 0.35
#------------------------------------------------------------
print("\t\tConvert mesh with 1 hole r=0.35")
# Read the mesh with one hole using meshio
msh = meshio.read("../Meshes/s_h_1_r_0p35.msh")
# Create a triangle mesh using the function "create_mesh"
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
# Save the triangle mesh
meshio.write("../Meshes/s_h_1_r_0p35.xdmf", triangle_mesh)
print("\t\t\tDone!")
#------------------------------------------------------------
# 1 HOLE RADIUS 0.4
#------------------------------------------------------------
print("\t\tConvert mesh with 1 hole r=0.4")
# Read the mesh with one hole using meshio
msh = meshio.read("../Meshes/s_h_1_r_0p4.msh")
# Create a triangle mesh using the function "create_mesh"
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
# Save the triangle mesh
meshio.write("../Meshes/s_h_1_r_0p4.xdmf", triangle_mesh)
print("\t\t\tDone!")
#------------------------------------------------------------
# 1 HOLE RADIUS 0.45
#------------------------------------------------------------
print("\t\tConvert mesh with 1 hole r=0.45")
# Read the mesh with one hole using meshio
msh = meshio.read("../Meshes/s_h_1_r_0p45.msh")
# Create a triangle mesh using the function "create_mesh"
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
# Save the triangle mesh
meshio.write("../Meshes/s_h_1_r_0p45.xdmf", triangle_mesh)
print("\t\t\tDone!")
#------------------------------------------------------------
# 1 HOLE RADIUS 0.5
#------------------------------------------------------------
print("\t\tConvert mesh with 1 hole r=0.5")
# Read the mesh with one hole using meshio
msh = meshio.read("../Meshes/s_h_1_r_0p5.msh")
# Create a triangle mesh using the function "create_mesh"
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
# Save the triangle mesh
meshio.write("../Meshes/s_h_1_r_0p5.xdmf", triangle_mesh)
print("\t\t\tDone!")
#------------------------------------------------------------
# 1 HOLE RADIUS 0.55
#------------------------------------------------------------
print("\t\tConvert mesh with 1 hole r=0.55")
# Read the mesh with one hole using meshio
msh = meshio.read("../Meshes/s_h_1_r_0p55.msh")
# Create a triangle mesh using the function "create_mesh"
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
# Save the triangle mesh
meshio.write("../Meshes/s_h_1_r_0p55.xdmf", triangle_mesh)
print("\t\t\tDone!")
