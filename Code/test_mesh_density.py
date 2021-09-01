# =================================================================================
# =================================================================================
# Script:"FEM_simulations"
# Date: 2021-08-31
# Implemented by: Johannes Borgqvist
# Description:
# This is the script where we will run the FEM simulations
# =================================================================================
# =================================================================================
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
# DOLFIN RELATED LIBRARIES 
from fenics import *
import numpy as np
# =================================================================================
# =================================================================================
# Conduct calculations
# =================================================================================
# =================================================================================
# SANITY CHECK: MESH WITH A HOLE
# Allocate memory for the main mesh with surfaces
mesh = Mesh()
mvc_subdomain = MeshValueCollection("size_t", mesh, 2)
mvc_curves = MeshValueCollection("size_t", mesh,1)
# Read in the mesh with subdomains
with XDMFFile("../Meshes/sphere_with_1_hole_surfaces.xdmf") as infile:
    infile.read(mesh)
    infile.read(mvc_subdomain, "name_to_read")
    #tag_info = infile.read("information_int")
# Read in the curves as well
with XDMFFile("../Meshes/sphere_with_1_hole_curves.xdmf") as infile:
    infile.read(mvc_curves, "name_to_read")
# Define meshfunctions for the surfaces and curves respectively
mf_triangle_sphere = cpp.mesh.MeshFunctionSizet(mesh, mvc_subdomain)
mf_curves = cpp.mesh.MeshFunctionSizet(mesh, mvc_curves)
# Define measures based on this as well
dx_sphere = Measure("dx", domain=mesh, subdomain_data=mf_triangle_sphere, subdomain_id=2)
dx_adjacent = Measure("dx", domain=mesh, subdomain_data=mf_triangle_sphere, subdomain_id=3)
dx_hole = Measure("dx", domain=mesh, subdomain_data=mf_triangle_sphere, subdomain_id=4)
ds_hole = Measure("ds", domain=mesh, subdomain_data=mf_curves, subdomain_id=1)
# Calculate the area of stuff
area_sphere = Constant(1.0)*dx_sphere
area_adjacent = Constant(1.0)*dx_adjacent
area_hole = Constant(1.0)*dx_hole
circumference_hole = Constant(1.0)*ds_hole
print("Area of hole")
print(assemble(area_sphere))
print("Area of adjacent")
print(assemble(area_adjacent))
print("Area of sphere")
print(assemble(area_hole))
print("Circumference of hole (curve integral)")
print(circumference_hole)
# SANITY CHECK: MESH WITHOUT A HOLE
# Allocate memory for the main mesh with surfaces
mesh_unit = Mesh()
mvc_subdomain_unit = MeshValueCollection("size_t", mesh_unit, 2)
# Read in the mesh with subdomains
with XDMFFile("../Meshes/sphere_with_no_holes_surfaces.xdmf") as infile:
    infile.read(mesh_unit)
    infile.read(mvc_subdomain_unit, "name_to_read")
# Define meshfunctions for the surfaces and curves respectively
mf_triangle_unit_sphere = cpp.mesh.MeshFunctionSizet(mesh_unit, mvc_subdomain_unit)
# Define measures based on this as well
dx_unit_sphere = Measure("dx", domain=mesh_unit, subdomain_data=mf_triangle_unit_sphere, subdomain_id=1)
# Calculate the area of stuff
area_unit_sphere = Constant(1.0)*dx_unit_sphere
# Print as a comaprison
print("The surface area of the numerical unit sphere:")
print(assemble(area_unit_sphere))
print("The analytical surface area of the unit sphere A=4*pi*1^2:")
print(4*np.pi)
