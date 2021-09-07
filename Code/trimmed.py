import gmsh
import sys
import math

gmsh.initialize(sys.argv)

# create sphere and get its boundary
v0 = gmsh.model.occ.addSphere(0,0,0, 2)
gmsh.model.occ.synchronize()
s0 = gmsh.model.getBoundary([(3, v0)])[0][1]
gmsh.model.occ.remove([(3, v0)])

# cut with a cylinder to generate first spherical patch
v1 = gmsh.model.occ.addCylinder(0,0,0, 0,2,2, 0.5)
gmsh.model.occ.intersect([(2, s0)], [(3, v1)], removeObject=False)

# create a wire in the parametric plane of the spherical surface
# [-pi,pi]x[-pi/2,pi/2] to create the other spherical patch
s3 = gmsh.model.occ.addRectangle(math.pi/2-0.5,0.5,0, 1,0.75)
gmsh.model.occ.synchronize()
b3 = gmsh.model.getBoundary([(2, s3)])
w3 = gmsh.model.occ.addWire([p[1] for p in b3])
gmsh.model.occ.addTrimmedSurface(s0, [w3])
gmsh.model.occ.remove([(2, s3)], recursive=True)

# fragment all surfaces to make everything conformal
gmsh.model.occ.removeAllDuplicates()

gmsh.model.occ.synchronize()

# set mesh size
gmsh.option.setNumber('Mesh.MeshSizeMax', 0.1)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
