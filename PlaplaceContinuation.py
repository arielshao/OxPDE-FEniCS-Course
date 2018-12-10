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
for p_val in [2, 3, 4, 5]:
    p.assign(p_val)
    solve(F == 0, u, bc, solver_parameters={"newton_solver": {"maximum_iterations": 100}})
    File("plapace.pvd") << u
