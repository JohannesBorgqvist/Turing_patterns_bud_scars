"""
Test for DG advection on sphere. Adapted from
dg-advection-diffusion demo by Johan Hake.
"""

# Copyright (C) 2008 Johan Hake
# Copyright (C) 2013 Marie E. Rognes, Imperial College London, Andrew McRae
# Modified by Colin J. Cotter 2013.

from dolfin import *
from ufl.log import info_blue
import numpy

# Define global normal
global_normal = Expression(("x[0]", "x[1]", "x[2]"))

# Define meshes
meshes = ["sphere_ico3", "sphere_ico4", "sphere_ico5", "sphere_ico6"]

# Set-up some storage for the errors. NB: Storing one error at a time
# in case things fail
errorfilename = "results/errors_DGadvection.py"
info_blue("Storing errors to %s" % errorfilename)
errorfile = open(errorfilename, "w")
errorfile.write("# Results from DG advection testcase.\n")
errorfile.write("L2_error_D = [0.0 for i in range(%d)];\n" % len(meshes))
errorfile.write("h_max = [0.0 for i in range(%d)];\n" % len(meshes))
errorfile.write("conservation_error = [0.0 for i in range(%d)];\n" % len(meshes))

errors = numpy.zeros((4,))

for (i,meshid) in enumerate(meshes):
    meshname = "../mixed-poisson/hdiv-l2/meshes/%s.xml.gz" % meshid
    mesh = Mesh(meshname)
    mesh.init_cell_orientations(global_normal)

    # Define function spaces and basis functions
    V_dg = FunctionSpace(mesh, "DG", 1)
    M = VectorFunctionSpace(mesh, "CG", 1)

    # advecting velocity
    u0 = Expression(('-x[1]','x[0]','0'))
    u = interpolate(u0,M)
    
    # Test and trial functions
    phi   = TestFunction(V_dg)
    D = TrialFunction(V_dg)
    
    # Mesh-related functions
    n = FacetNormal(mesh)

    # ( dot(v, n) + |dot(v, n)| )/2.0
    un = (dot(u, n) + abs(dot(u, n)))/2.0

    # un is positive if u is pointing outwards and negative if pointing inwards
    # un('+') = -un('-') if u is continuous.

    # parameters
    dt = numpy.pi/3*0.01
    
    # Bilinear form
    a_mass = phi*D*dx
    a_int = dot(grad(phi), -u*D)*dx
    a_flux = (dot(jump(phi), un('+')*D('+') - un('-')*D('-') )*dS  
              + dot(phi, un*D)*ds)

    M = assemble(a_mass)
    arhs = -dt*(a_int + a_flux)
    
    D0 = Expression("exp(-pow(x[2],2) - pow(x[1],2))")
    D = interpolate(D0, V_dg)
    
    initial_tracer_mass = assemble(D*dx)

    dD1 = Function(V_dg)
    D1 = Function(V_dg)
    dD2 = Function(V_dg)
    D2 = Function(V_dg)
    dD3 = Function(V_dg)
    
    t = 0.0
    T = 2*pi
    k = 0
    dumpfreq = 20

    Dfile = File("results/pvds/D.pvd")
    while(t < (T-dt/2)):
        L = assemble(action(arhs,D))
        solve(M,dD1.vector(),L)

        D1.vector()[:] = D.vector().copy()
        D1.vector().axpy(1.0, dD1.vector())
        L = assemble(action(arhs,D1))
        solve(M,dD2.vector(),L)

        D2.vector()[:] = D.vector().copy()
        D2.vector().axpy(0.25, dD1.vector())
        D2.vector().axpy(0.25, dD2.vector())
        L = assemble(action(arhs,D2))
        solve(M,dD3.vector(),L)

        D.vector().axpy((1.0/6.0), dD1.vector())
        D.vector().axpy((1.0/6.0), dD2.vector())
        D.vector().axpy((2.0/3.0), dD3.vector())
        # Store solutions to xml and pvd
        t +=dt
        k += 1
        if(k % dumpfreq == 0):
            solutionfile = File("results/xmls/D_%d.xml.gz" % k)
            solutionfile << D
            Dfile << D

    # Compute errors
    info_blue("Computing errors")
    L2_error_D = errornorm(D0,D)
    conservation_error = assemble(D*dx) - initial_tracer_mass
    
    # Store errors
    info_blue("Storing errors")
    errorfile.write("h_max[%d] = %g;\n" % (i, mesh.hmax()))
    errorfile.write("L2_error_D[%d] = %g;\n" % (i, L2_error_D))
    errorfile.write("conservation_error[%d] = %g;\n" % (i, conservation_error))
    errorfile.flush()

errorfile.close()

#h = numpy.array([1.0,1.0/2,1.0/4,1.0/8])
#rate = numpy.log(errors[:3]/errors[1:])/numpy.log(h[:3]/h[1:])
#print rate

