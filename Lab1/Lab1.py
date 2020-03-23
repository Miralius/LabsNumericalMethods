import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols
from sympy import S

x = symbols('x')


def plot_function(func, a, b):
    x_num = numpy.linspace(a, b, 1000)
    func_numpy = numpy.zeros(1000)
    i = 0
    while i < 1000:
        func_numpy[i] = func.subs(x, x_num[i])
        i += 1
    y_max = func_numpy[999]
    y_min = func_numpy[0]
    xy = matplotlib.pyplot.subplot()
    xy.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator((b - a) / 5))
    xy.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator((b - a) / 25))
    xy.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator((y_max - y_min) / 5))
    xy.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator((y_max - y_min) / 25))
    xy.grid(which='major', color='k')
    xy.minorticks_on()
    xy.grid(which="minor", color='gray', linestyle=':')
    xy.plot(x_num, func_numpy, color="black", label="y(x)")
    xy.set_xlabel("x")
    xy.set_ylabel("y")
    xy.legend()
    matplotlib.pyplot.show()


def find_solve(func, a, b, eps):
    root = 0
    count = 0
    while True:
        count += 1
        f_a = func.subs(x, a)
        c = 0.5 * (a + b)
        if (b - a) < 2 * eps:
            root = c
            break
        f_c = func.subs(x, c)
        if f_c == 0:
            root = c
            break
        if f_a * f_c < 0:
            b = c
        else:
            a = c
            f_a = f_c
    print(root)


plot_function(x ** 3 + 8.5 * x ** 2 + 21.8 * x + S(15.6), -10, 10)
