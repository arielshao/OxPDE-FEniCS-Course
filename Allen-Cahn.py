from dolfin import *

mesh = UnitSquareMesh(100, 100)

V = FunctionSpace(mesh, 'CG',1)
u = Function(V)
v = TestFunction(V)
f = Constant(0.0)
delta=Constant(0.04)

F = delta* inner(grad(u), grad(v))*dx +1/delta*inner(u**3-u,v)*dx

bclr=DirichletBC(V,+1, "x[0]==0.0 || x[0]==1.0")
bctb=DirichletBC(V,-1, "x[1]==0.0 || x[1]==1.0")


pvd=File("Allen-Cahn.pvd")

for ig in [0,1,-1]:
	u.interpolate(Constant(ig))
	solve(F==0, u, [bclr,bctb])
	pvd<<u

