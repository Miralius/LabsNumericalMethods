import numpy
from sympy import symbols, exp, sqrt, diff, N, solve


# noinspection SpellCheckingInspection
def m_n_plus_one(n, a, b, func):
    x = symbols('x')
    extremums = solve(diff(func, x, n + 2))
    maximum = 0
    maximum = max(maximum, abs(N(diff(func, x, n + 1).subs(x, a))))
    maximum = max(maximum, abs(N(diff(func, x, n + 1).subs(x, b))))
    i = 0
    while i < len(extremums):
        maximum = max(maximum, abs(N(diff(func, x, n + 1).subs(x, extremums.args[i]))))
        i += 1
    return float(maximum)


def function(x):
    return numpy.exp(-numpy.sqrt(x))


def define_step(a, b, eps):
    x = symbols('x')
    return float(numpy.sqrt(12 * eps / (m_n_plus_one(1, a, b, exp(-sqrt(x))) * (b - a))))
