from dolfin import *

mesh= UnitIntervalMesh(500)
element = FiniteElement("Lagrange", interval, 1)
V=FunctionSpace(mesh, element)


u=Function(V)
v=TestFunction(V)
lmbda = Constant(2.0)
F=-inner(grad(u),grad(v))*dx+lmbda*exp(u)*v*dx
bc=DirichletBC(V,0, DomainBoundary())

pvd=File("Bratu_equation.pvd")
for u_init in [0, 3]:
	u.interpolate(Constant(u_init))
	solve(F == 0, u , bc) # apply Newton-Kantorovich
	#plot(u, interactive = True)
	pvd<<u
