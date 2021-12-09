from dolfin import *
from ufl.operators import contraction

mesh = Mesh("circle.xml")
V = VectorFunctionSpace(mesh, "P", 1)

T = 20
h = 1/20
dt = h / 2
alpha_1 = 0.01
alpha_2 = 0.005
...


class InitialConditions(UserExpression):
    def eval(self, values, x):
        values[0] = 1
        values[1] = 2

    def value_shape(self):
        return (2,)

indata = InitialConditions(2)
u0 = Function(V)
u0 = interpolate(indata, V)

a = ...
L = ...

file = File("output.pvd")


t = 0
while t<T:
    u0.assign(u)

