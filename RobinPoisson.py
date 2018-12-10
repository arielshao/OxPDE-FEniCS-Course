from dolfin import *

mesh = UnitSquareMesh(32,32)
element = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, element)
u = Function(V)
v = TestFunction(V)

f = Constant(1.0)
g = Constant(0.0)
bc = DirichletBC(V, g, "x[0] >0 && on_boundary")

dim = 2
colors = MeshFunction("size_t", mesh, dim -1)
colors.set_all(0) 
CompiledSubDomain("x[0] == 0").mark(colors, 1)
File("colors.pvd") << colors

ds = Measure("ds", subdomain_data=colors)

F = inner(grad(u), grad(v))*dx + inner(u,v)*ds(1)- f*v*dx
    
solve(F == 0, u, bc)
File("RobinPoisson.pvd") << u
