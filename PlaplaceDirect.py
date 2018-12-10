from dolfin import *

mesh = UnitSquareMesh(32,32)
element = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, element)
u = Function(V)
v = TestFunction(V)

f = Constant(1.0)
g = Constant(0.0)
p = Constant(5.0)
epsilon = Constant(1.0e-5)
gamma = (epsilon**2+0.5*inner(grad(u),grad(u)))**((p-2)/2)
bc = DirichletBC(V, g, DomainBoundary())
F = inner(grad(v),gamma*grad(u))*dx-inner(f,v)*dx
u.interpolate(Constant(0))
solve(F == 0, u, bc)
File("plapace.pvd") << u
