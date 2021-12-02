from dolfin.mesh import subdomain
import matplotlib.pyplot as plt
import math
from dolfin import *
from mshr import *
import numpy as np
# Create mesh and define function space
mesh = Mesh("circle.xml")
Q = FunctionSpace(mesh, "CG", 1)
# Define parameters:
T = 10
h = mesh.hmin()
print(h)
r = .2
R = .5
dt = h
alpha = 0.01


# Create subdomain for Dirichlet boundary


def boundary(x):
    return np.isclose(np.sqrt(x[0]**2+x[1]**2), r) or np.isclose(np.sqrt(x[0]**2+x[1]**2), R)


# Set up boundary condition
g = Constant(0.0)
bc = DirichletBC(Q, g, boundary)
# Define initial condition
rho = 100
indata = Expression(
    "abs(R - pow(pow(x[0], 2) + pow(x[1], 2), 0.5)) <= r ? rho:0", rho=rho, R=R, r=r, degree=2)
u0 = Function(Q)
u0 = interpolate(indata, Q)


u = TrialFunction(Q)
v = TestFunction(Q)
# Create bilinear and linear forms
f = Constant(0)
a = u * v * dx + (dt * alpha/2) * inner(grad(u), grad(v)) * dx
L = a

problem = LinearVariationalProblem(lhs(a), rhs(a), u, bc)
solver = LinearVariationalSolver(problem)

t = 0
# Set an output file
file = File("test.pvd")

# Set initial condition
u = Function(Q)
u.assign(u0)

# Time-stepping
A = assemble(a)
while t < T:
    # # assign u0
    u.assign(u)
    b = assemble(L)

    # solving
    solver.solve()
    t += dt
    print(t)

file << (u, t)
plt.figure()

plot(u)
plt.show()
