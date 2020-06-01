import numpy
from sympy import symbols, exp, sqrt, diff, N


def function(x):
    return numpy.exp(-numpy.sqrt(x))


def define_step(a, b, eps):
    x = symbols('x')
    return float(N(2 * ((eps / (abs(diff(exp(-sqrt(x)), x, 2)).subs(x, a))) ** (1 / 4))))