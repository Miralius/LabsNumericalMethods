import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols, pi, diff, N

x = symbols('x')


def function(point):
    return numpy.math.pi * point / (10 + point / 2)


def define_step(func):
    return N(2 * ((eps / abs(diff(func, x, 4)).subs(x, a)) ** (1 / 4)))


# noinspection SpellCheckingInspection
def discretize():
    n = numpy.ceil((b - a) / h)
    return numpy.array([a + i * h for i in range(n)]), numpy.array([function(a + i * h) for i in range(n)])


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
