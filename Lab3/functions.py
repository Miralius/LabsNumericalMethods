import numpy
from sympy import symbols, exp, sqrt, diff, N, solve, integrate


# noinspection SpellCheckingInspection
def m_n_plus_one(n, a, b, func):
    x = symbols('x')
    extremums = solve(diff(func, x, n + 2))
    maximum = 0
    maximum = max(maximum, abs(N(diff(func, x, n + 1).subs(x, a))))
    maximum = max(maximum, abs(N(diff(func, x, n + 1).subs(x, b))))
    for number in extremums:
        maximum = max(maximum, abs(N(diff(func, x, n + 1).subs(x, number))))
    return float(maximum)


def function(x):
    return numpy.exp(-numpy.sqrt(x))


def define_step(a, b, eps):
    x = symbols('x')
    h = float(numpy.sqrt(12 * eps / (m_n_plus_one(1, a, b, exp(-sqrt(x))) * (b - a))))
    n = numpy.ceil((b - a) / h)
    n_remainder = n % 4
    if n_remainder != 0:
        n += 4 - n_remainder
    h = (b - a) / n
    return h


def trapezes_integrate(a, b, h):
    n = numpy.round((b - a) / h)
    integral = 0
    i = 1
    while i <= n:
        integral += (h / 2) * (function(a + (i - 1) * h) + function(a + i * h))
        i += 1
    return integral


def simpson_integrate(a, b, h):
    n = numpy.round((b - a) / h)
    integral = 0
    i = 1
    while i <= n:
        integral += (h / 6) * (function(a + (i - 1) * h) + 4 * function(a + (i - 1 / 2) * h) + function(a + i * h))
        i += 1
    return integral


def newton_leibniz_integrate(a, b):
    x = symbols('x')
    return integrate(exp(-sqrt(x)), (x, a, b))
