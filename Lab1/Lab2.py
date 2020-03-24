import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols, diff, solveset, S, Eq, cos, sin, solve

x, y = symbols('x y')


def plot_functions(func1, func2, a, b):
    x_num = numpy.linspace(a, b, 1000)
    func1 = solve(func1, y)
    func2 = solve(func2, x)
    func1_numpy = numpy.zeros(1000)
    func2_numpy = numpy.zeros(1000)
    i = 0
    while i < 1000:
        func1_numpy[i] = func1[0].subs(x, x_num[i])
        i += 1
    i = 0
    while i < 1000:
        func2_numpy[i] = func2[0].subs(y, x_num[i])
        i += 1
    xy = matplotlib.pyplot.subplot()
    xy.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator((b - a) / 5))
    xy.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
    xy.grid(which='major', color='k')
    xy.plot(x_num, func1_numpy, color="blue", label="f1(x)")
    xy.plot(func2_numpy, x_num, color="red", label="f2(x)")
    xy.set_xlabel("x")
    xy.set_ylabel("y")
    xy.legend()
    matplotlib.pyplot.show()


plot_functions(Eq(cos(x + 0.5) - y, 2), Eq(sin(y) - 2 * x, 1), -10, 10)
plot_functions(Eq(cos(x + 0.5) - y, 2), Eq(sin(y) - 2 * x, 1), -4, 0)