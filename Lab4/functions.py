import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols


def function(x, y):
    return y * (0.5 * x * y - 1)


def func():
    x, y = symbols('x y')
    return y * (0.5 * x * y - 1)


def delta(x, y, h):
    phi0 = h * function(x, y)
    phi1 = h * function(x + h / 2, y + phi0 / 2)
    phi2 = h * function(x + h / 2, y + phi1 / 2)
    phi3 = h * function(x + h, y + phi2)
    return (phi0 + 2 * phi1 + 2 * phi2 + phi3) / 6


# noinspection DuplicatedCode
def define_step(h, a, b, y0, x0):
    print("Определение шага…")
    h0 = h
    while True:
        y1 = y0 + delta(x0, y0, h0)
        y2 = y1 + delta(x0 + h0, y1, h0)
        y2_ = y0 + delta(x0, y0, 2 * h0)
        dyy2 = abs(y2 - y2_) / 15
        print("Шаг h = " + str(h0))
        print("y1 = " + str(y1))
        print("y2 = " + str(y2))
        print("Двойной шаг, y2_ = " + str(y2_))
        print("Погрешность y2 - y2_ = " + str(dyy2))
        if dyy2 < 10 ** (-4):
            print("Увеличиваем шаг…")
            h0 *= 2
        else:
            print("Шаг для решения выбран!")
            break
    n = numpy.ceil((b - a) / h0)
    print("Определим n = " + str(n))
    if n % 2:
        print("n нечетное, выбираем n = " + str(n + 1))
        n += 1
    h0 = (b - a) / n
    y1 = y0 + delta(x0, y0, h0)
    y2 = y1 + delta(x0 + h0, y1, h0)
    y2_ = y0 + delta(x0, y0, 2 * h0)
    dyy2 = abs(y2 - y2_) / 15
    print("Шаг h = " + str(h0))
    print("y1 = " + str(y1))
    print("y2 = " + str(y2))
    print("Двойной шаг, y2_ = " + str(y2_))
    print("Погрешность y2 - y2_ = " + str(dyy2))
    return h0


def discretize(a, b, h):
    return numpy.array([a + i * h for i in range(int(round((b - a) / h)) + 1)])


# noinspection SpellCheckingInspection
def runge_kutta_dsolve(x, y0, a, b, h):
    n = numpy.ceil((b - a) / h) + 1
    y = numpy.zeros(int(n))
    y[0] = y0
    i = 1
    while i < n:
        y[i] = y[i - 1] + delta(x[0] + (i - 1) * h, y[i - 1], h)
        i += 1
    return y


# noinspection SpellCheckingInspection
def eiler_dsolve(x, y0, a, b, h):
    n = numpy.ceil((b - a) / h) + 1
    y = numpy.zeros(int(n))
    y[0] = y0
    i = 1
    while i < n:
        y[i] = y[i - 1] + h * function(x[0] + (i - 1) * h, y[i - 1])
        i += 1
    return y


def plot(x, functions, names, axis_values=None):
    xy = matplotlib.pyplot.subplot()
    xy.grid(which='major', color='k')
    xy.minorticks_on()
    xy.grid(which="minor", color='gray', linestyle=':')
    i = 0
    while i < len(functions):
        xy.plot(x, functions[i], label=names[i])
        i += 1
    xy.set_xlabel("x")
    xy.set_ylabel("y")
    if axis_values is not None:
        matplotlib.pyplot.axis(axis_values)
    xy.legend()
    matplotlib.pyplot.show()
