from dolfin import *
from mshr   import *
from math   import atan2


class BoundaryData(UserExpression):
	def eval(self, values, x):
		r = sqrt(x[0]**2+x[1]**2)
		theta =atan2(x[1], x[0])

		if theta < 0: theta += 2* pi
		values[0]=r**(2.0/3.0)*sin((2.0/3.0)*theta)
	
square=Rectangle (Point(-1, -1), Point (1,1))
cutout=Rectangle(Point(+0, -1), Point(1,0))
domain = square -cutout
for mesh_size in [50, 100, 200]:
	mesh = generate_mesh(domain,mesh_size)
	element = FiniteElement("Lagrange", mesh.ufl_cell(), 2)
	V = FunctionSpace(mesh, element)
	u = Function(V)
	v = TestFunction(V)
	f = Constant(0.0)
	g=BoundaryData(degree=5)
	bc = DirichletBC(V, g, DomainBoundary())
	F = inner(grad(u), grad(v))*dx - f*v*dx
	solve(F == 0, u, bc)
	print( errornorm (g, u, 'L2'))
	print( errornorm (g, u, 'H1')) 
	
	#File("LshapePoisson.pvd") << u