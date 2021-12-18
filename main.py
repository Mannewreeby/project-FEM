import numpy as np
from dolfin import *
import matplotlib.pyplot as plt

mesh = Mesh("2d_mesh/circle2.xml")

V = VectorFunctionSpace(mesh, "P", 1)

T = 20
dt = 1/20
alpha_1 = Constant(0.01)
alpha_2 = Constant(0.005)
cf = Constant(0.024)
ck = Constant(0.055)
rho = Constant(2)
R = 0.5
r = 0.3


class InitialCondition(UserExpression):
    def eval(self, values, x):
        if (abs(R - sqrt(x[0]**2 + x[1]**2)) <= r):
            values[0] = rho/2*np.random.rand()
            values[1] = rho*(1 + 0.5*np.random.rand())
        else:
            values[0] = 0
            values[1] = 0

    def value_shape(self):
        return (2, )


indata = InitialCondition(degree=2)
# Homogenous Neumann is standard boundary condition in FEniCS
u0 = Function(V)
u0 = interpolate(indata, V)

plt.ion()
plt.figure("u initial")
initU = plot(u0[0])
plt.colorbar(initU)
plt.figure("v initial")
initV = plot(u0[1])
plt.colorbar(initV)
plt.show()

u = TrialFunction(V)
v = TestFunction(V)


# a1 = 2 * u[0]*v[0]*dx + dt * \
#     (cf*u[0]*v[0] + alpha_1*inner(grad(u[0]), grad(v[0])))*dx
# a2 = 2 * u[1]*v[1]*dx + dt*((cf + ck)*u[1]*v[1] +
#                             alpha_2*inner(grad(u[1]), grad(v[1])))*dx
# L1 = 2*dt*(cf*v[0]*dx - S*v[0]*dx) + 2*u0[0]*v[0]*dx - dt * \
#     (cf*u0[0]*v[0] + alpha_1*inner(grad(u0[0]), grad(v[0])))*dx
# L2 = 2*dt*(S*v[1]*dx) + 2*u0[1]*v[1]*dx - dt*((cf + ck)*u0[1]
#                                               * v[1] + alpha_2*inner(grad(u0[1]), grad(v[1])))*dx

# Treat non linear term explicitly
S = u0[0]*u0[1]*u0[1]

# Bilinear and linear forms
a1 = u[0]*v[0]*dx + dt/2*(cf*u[0]*v[0]*dx + alpha_1 *
                          inner(grad(u[0]), grad(v[0]))*dx)
a2 = u[1]*v[1]*dx + dt/2*((cf + ck)*u[1]*v[1]*dx +
                          alpha_2*inner(grad(u[1]), grad(v[1]))*dx)

L1 = u0[0]*v[0]*dx - dt/2 * \
    (cf*u0[0]*v[0]*dx + alpha_1*inner(grad(u0[0]), grad(v[0]))
     * dx) + dt*(cf*v[0]*dx - S*v[0]*dx)
L2 = u0[1]*v[1]*dx - dt/2*((cf + ck)*u0[1]*v[1]*dx +
                           alpha_2*inner(grad(u0[1]), grad(v[1]))*dx) + dt*S*v[1]*dx

a = a1 + a2
L = L1 + L2


# Setup start and stop times
T = 20
t = 0

# Time step
dt = .05
time = np.linspace(0, T, int(T/dt))

# Setup u and assign initial condition
u = Function(V)
u.assign(u0)

# Prepare mass loss integrals
Mu = (u0[0] - u[0])*dx
Mv = (u0[1] - u[1])*dx
massLoss = np.zeros((2, len(time)))
# Assemble for t=0
massLoss[:, 0] = [assemble(Mu), assemble(Mv)]


# Time stepping
for i, step in enumerate(time[1:]):
    u0.assign(u)
    solve(a == L, u)
    massLoss[:, i + 1] = [assemble(Mu), assemble(Mv)]
    # plot(u[0])
    # plt.show()
    t += dt

# Save u to file
file = File("result.pvd")
file << u


# Plot result
plt.figure("u final")
p = plot(u[0])
plt.colorbar(p)
plt.figure("v final")
z = plot(u[1])
plt.colorbar(z)
plt.show()


# Plot mass loss
plt.figure("mass loss")
plt.xlabel('time')
plt.ylabel('u')
plt.title('u and v mass loss')
uPlot = plt.plot(time, massLoss[0, :], 'b', label='u')
vPlot = plt.plot(time, massLoss[1, :], 'r', label='v')

plt.show()


close = input("press enter to close plots and terminate program")
plt.close("all")
