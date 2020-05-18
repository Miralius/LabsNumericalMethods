import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols, pi, diff, N

x = symbols('x')


def function(point):
    return numpy.math.pi * point / (10 + point / 2)


def define_step(func):
    return N(2 * ((eps / abs(diff(func, x, 4)).subs(x, a)) ** (1 / 4)))


def discretize():
    n = numpy.ceil((b - a) / h)
    return numpy.array([a + i * h for i in range(n)]), numpy.array([function(a + i * h) for i in range(n)])


def interpolate_with_the_Lagrange_polynomial(x, otr):
    i = 0
    L = 0
    while i < 4:
        j = 0
        mult = 1
    while j < 4:
        if i != j:
            mult *= (x - otr[j]) / (otr[i] - otr[j])
        j += 1
        # L += f(otr[i]) * mult
        i += 1
    return L


if __name__ == "__main__":
    print("Лабораторная работа №3, вариант #14")
    a = 0
    b = 10
    eps = 10 ** (-4)
    print("a = " + str(a))
    print("b = " + str(b))
    print("Точность = " + str(eps))
    h = define_step(pi * x / (10 + x / 2))
    print("Шаг = " + str(h))
    x, y = discretize()
    print("x = " + str(x))
