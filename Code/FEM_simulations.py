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
