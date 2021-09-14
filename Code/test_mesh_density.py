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
# CONDUCT SANITY CHECK
# =================================================================================
# =================================================================================
# We check the meshes by importing them to FEniCS and then we calculate the total area of all surfaces in the mesh. This is then compared to the analytical value of the area of the unit sphere given by 4*pi. 
#-------------------------------------------------------------
# ANALYTICAL AREA
#-------------------------------------------------------------
# Calculate the analytical area
area_analytical = 4*np.pi
#-------------------------------------------------------------
# NO HOLES
#-------------------------------------------------------------
# Allocate memory for the main mesh with surfaces
mesh = Mesh()
mvc_subdomain = MeshValueCollection("size_t", mesh, 2)
# Read in the mesh with subdomains
with XDMFFile("../Meshes/s_h_0.xdmf") as infile:
    infile.read(mesh)
    infile.read(mvc_subdomain, "name_to_read")
# Define meshfunctions for the surfaces and curves respectively
mf_h_0 = cpp.mesh.MeshFunctionSizet(mesh, mvc_subdomain)
# Define measures based on this as well
dx_sphere_h_0 = Measure("dx", domain=mesh, subdomain_data=mf_h_0, subdomain_id=1)
# Calculate the area of stuff
area_sphere_0 = assemble(Constant(1.0)*dx_sphere_h_0)
#-------------------------------------------------------------
# 1 HOLE
#-------------------------------------------------------------
# Allocate memory for the main mesh with surfaces
mesh = Mesh()
mvc_subdomain = MeshValueCollection("size_t", mesh, 2)
# Read in the mesh with subdomains
with XDMFFile("../Meshes/s_h_1.xdmf") as infile:
    infile.read(mesh)
    infile.read(mvc_subdomain, "name_to_read")
# Define meshfunctions for the surfaces and curves respectively
mf_h_1 = cpp.mesh.MeshFunctionSizet(mesh, mvc_subdomain)
# Define measures based on this as well
dx_sphere_h_1 = Measure("dx", domain=mesh, subdomain_data=mf_h_1, subdomain_id=1)
dx_adjacent_h_1 = Measure("dx", domain=mesh, subdomain_data=mf_h_1, subdomain_id=2)
dx_hole_h_1 = Measure("dx", domain=mesh, subdomain_data=mf_h_1, subdomain_id=3)
# Calculate the area of stuff
area_sphere_1 = assemble(Constant(1.0)*dx_sphere_h_1)
area_adjacent_1 = assemble(Constant(1.0)*dx_adjacent_h_1)
area_hole_1 = assemble(Constant(1.0)*dx_hole_h_1)
area_total_1 = area_sphere_1 + area_adjacent_1 + area_hole_1
#-------------------------------------------------------------
# 2 HOLES
#-------------------------------------------------------------
# Allocate memory for the main mesh with surfaces
mesh = Mesh()
mvc_subdomain = MeshValueCollection("size_t", mesh, 2)
# Read in the mesh with subdomains
with XDMFFile("../Meshes/s_h_2.xdmf") as infile:
    infile.read(mesh)
    infile.read(mvc_subdomain, "name_to_read")
# Define meshfunctions for the surfaces and curves respectively
mf_h_2 = cpp.mesh.MeshFunctionSizet(mesh, mvc_subdomain)
# Define measures based on this as well
dx_sphere_h_2 = Measure("dx", domain=mesh, subdomain_data=mf_h_2, subdomain_id=1)
dx_adjacent_h_2 = Measure("dx", domain=mesh, subdomain_data=mf_h_2, subdomain_id=2)
dx_hole_h_2 = Measure("dx", domain=mesh, subdomain_data=mf_h_2, subdomain_id=3)
# Calculate the area of stuff
area_sphere_2 = assemble(Constant(1.0)*dx_sphere_h_2)
area_adjacent_2 = assemble(Constant(1.0)*dx_adjacent_h_2)
area_hole_2 = assemble(Constant(1.0)*dx_hole_h_2)
area_total_2 = area_sphere_2 + area_adjacent_2 + area_hole_2
#-------------------------------------------------------------
# 5 HOLES
#-------------------------------------------------------------
# Allocate memory for the main mesh with surfaces
mesh = Mesh()
mvc_subdomain = MeshValueCollection("size_t", mesh, 2)
# Read in the mesh with subdomains
with XDMFFile("../Meshes/s_h_5.xdmf") as infile:
    infile.read(mesh)
    infile.read(mvc_subdomain, "name_to_read")
# Define meshfunctions for the surfaces and curves respectively
mf_h_5 = cpp.mesh.MeshFunctionSizet(mesh, mvc_subdomain)
# Define measures based on this as well
dx_sphere_h_5 = Measure("dx", domain=mesh, subdomain_data=mf_h_5, subdomain_id=1)
dx_adjacent_h_5 = Measure("dx", domain=mesh, subdomain_data=mf_h_5, subdomain_id=2)
dx_hole_h_5 = Measure("dx", domain=mesh, subdomain_data=mf_h_5, subdomain_id=3)
# Calculate the area of stuff
area_sphere_5 = assemble(Constant(1.0)*dx_sphere_h_5)
area_adjacent_5 = assemble(Constant(1.0)*dx_adjacent_h_5)
area_hole_5 = assemble(Constant(1.0)*dx_hole_h_5)
area_total_5 = area_sphere_5 + area_adjacent_5 + area_hole_5
#-------------------------------------------------------------
# PRINT EVERYTHING
#-------------------------------------------------------------
print("\n\t\tANALYTICAL AREA:\n\t\t\tA_ana\t=\t%0.5f\n"%(area_analytical))
print("\n\t\tNO HOLES\n")
print("\t\t\tArea sphere:\tA_sphere\t=\t%0.5f\n"%(area_sphere_0))
print("\n\t\tONE HOLE\n")
print("\t\t\tArea sphere:\tA_sphere\t=\t%0.5f\n"%(area_sphere_1))
print("\t\t\tArea adjacent:\tA_adjacent\t=\t%0.5f\n"%(area_adjacent_1))
print("\t\t\tArea sphere:\tA_hole\t=\t%0.5f\n"%(area_hole_1))
print("\t\t\tArea sphere:\tA_total\t=\t%0.5f\n"%(area_total_1))
print("\n\t\tTWO HOLES\n")
print("\t\t\tArea sphere:\tA_sphere\t=\t%0.5f\n"%(area_sphere_2))
print("\t\t\tArea adjacent:\tA_adjacent\t=\t%0.5f\n"%(area_adjacent_2))
print("\t\t\tArea sphere:\tA_hole\t=\t%0.5f\n"%(area_hole_2))
print("\t\t\tArea sphere:\tA_total\t=\t%0.5f\n"%(area_total_2))
print("\n\t\tFIVE HOLES\n")
print("\t\t\tArea sphere:\tA_sphere\t=\t%0.5f\n"%(area_sphere_5))
print("\t\t\tArea adjacent:\tA_adjacent\t=\t%0.5f\n"%(area_adjacent_5))
print("\t\t\tArea sphere:\tA_hole\t=\t%0.5f\n"%(area_hole_5))
print("\t\t\tArea sphere:\tA_total\t=\t%0.5f\n"%(area_total_5))
