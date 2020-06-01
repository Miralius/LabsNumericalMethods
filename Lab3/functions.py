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
    h = float(numpy.sqrt(12 * eps / (m_n_plus_one(1, a, b, exp(-sqrt(x))) * (b - a))))
    n = numpy.ceil((b - a) / h)
    n_remainder = (n - 1) % 4
    if n_remainder != 0:
        n += 4 - n_remainder
    h = (b - a) / n
    return h


#def trapezoidal(a, b, h):
#    n = (b-a) / h
#    n += 4 - n % 4
#    value = 0
#    for   i = 0:n-1
#        value = value + (subs(f, a + i*h) + subs(f, a + (i+1)*h))*h/2
#    value = eval(value)
#    return value
