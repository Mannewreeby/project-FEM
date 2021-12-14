import numpy as np
from dolfin import *
import matplotlib.pyplot as plt

mesh = Mesh("2d_mesh/circle1.xml")

V = VectorFunctionSpace(mesh, "P", 2)

T = 20
dt = 1/20
alpha_1 = 0.01
alpha_2 = 0.005
cf = 0.024
ck = 0.055
rho = 2


class InitialCondition(UserExpression):
    def eval(self, values, x):
        values[0] = rho/2*np.random.rand()
        values[1] = rho*(1 + 0.5*np.random.rand())

    def value_shape(self):
        return (2, )


indata = InitialCondition(degree=2)
u0 = Function(V)
u0 = interpolate(indata, V)
S = u0[0]*u0[1]*u0[1]

u = TrialFunction(V)
v = TestFunction(V)

M1 = u[0]*v[0]*dx
M2 = u[1]*v[1]*dx
A1 = inner(grad(u[0]), grad(v[0]))*dx
A2 = inner(grad(u[1]), grad(v[1]))*dx
b = cf*v[0]*dx


a1 = 2 * M1 + dt*(cf*M1 + alpha_1*A1)
a2 = 2 * M2 + dt*((cf + ck)*M2 + alpha_2*A2)

L1 = 2*dt*(b - S*u0[0]*v[0]*dx) + (2*u0[0]*v[0]*dx - dt *
                                 (cf*u0[0]*v[0]*dx + alpha_1*inner(grad(u0[0]), grad(v[0]))*dx))
L2 = 2*dt*S*u0[1]*v[1]*dx + (2*u0[1]*v[1]*dx - dt*((cf + ck) * u0[1]*v[1]*dx + alpha_2*inner(grad(u0[1]), grad(v[1]))*dx))

a = a1 + a2
L = L1 + L2

T = 20
t = 0

dt = .5

u = Function(V)
u.assign(u0)


while t < T:
    u0.assign(u)
    solve(a == L, u)
    t += dt

plot(u)

plt.show()
