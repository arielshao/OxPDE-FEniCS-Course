from dolfin import *
mesh = UnitSquareMesh(32,32)
element = FiniteElement("Lagrange", mesh.ufl_cell(), 2)
V = FunctionSpace(mesh, element)
u = Function(V)
v = TestFunction(V)
f = Constant(-6.0)
g = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", element=V.ufl_element())
bc = DirichletBC(V, g, DomainBoundary())
F = inner(grad(u), grad(v))*dx - f*v*dx
solve(F == 0, u, bc)
File("poisson.pvd") << u
print(errornorm(g, u, 'L2'))
print(errornorm(g, u, 'H1'))