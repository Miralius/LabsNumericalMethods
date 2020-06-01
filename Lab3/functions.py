import numpy
from sympy import symbols, exp, sqrt, diff, N, solve, integrate


def function(x):
    return numpy.exp(-numpy.sqrt(x))


def function_symbolic():
    return exp(-sqrt(symbols('x')))


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


def define_step(a, b, eps):
    h = numpy.sqrt(12 * eps / (m_n_plus_one(1, a, b, function_symbolic()) * (b - a)))
    n = numpy.ceil((b - a) / h)
    n_remainder = n % 4
    if n_remainder != 0:
        n += 4 - n_remainder
    return (b - a) / n


def trapezes_integrate(a, b, h):
    n = numpy.ceil((b - a) / h)
    integral = 0
    i = 1
    while i <= n:
        integral += (h / 2) * (function(a + (i - 1) * h) + function(a + i * h))
        i += 1
    return integral


def simpson_integrate(a, b, h):
    n = numpy.ceil((b - a) / h)
    integral = 0
    i = 1
    while i <= n:
        integral += (h / 6) * (function(a + (i - 1) * h) + 4 * function(a + (i - 1 / 2) * h) + function(a + i * h))
        i += 1
    return integral


def newton_leibniz_integrate(a, b):
    return integrate(function_symbolic(), (symbols('x'), a, b))
