import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols, pi, diff, N


def function(x):
    return numpy.math.pi * x / (10 + x / 2)


# noinspection PyUnusedLocal
def define_step(a, b, eps):
    x = symbols('x')
    func = pi * x / (10 + x / 2)
    return float(N(2 * ((eps / (abs(diff(func, x, 4)).subs(x, a))) ** (1 / 4))))


def discretize(a, b, h):
    n = int(round((b - a) / h)) + 1
    return numpy.array([a + i * h for i in range(n)]), numpy.array([function(a + i * h) for i in range(n)])


def interpolate_with_lagrange_polynomial(x, interval):
    i = 0
    interpolated_function = 0
    while i < 4:
        j = 0
        multiplier = 1
        while j < 4:
            if i != j:
                multiplier *= (x - interval[j]) / (interval[i] - interval[j])
            j += 1
        interpolated_function += function(interval[i]) * multiplier
        i += 1
    return interpolated_function


def interpolate_with_newton_method(x, x_int, y_int_l, h):
    t = (x - x_int[0]) / h
    dy0 = y_int_l[1] - y_int_l[0]
    dy1 = y_int_l[2] - y_int_l[1]
    dy2 = y_int_l[3] - y_int_l[2]
    dy0_2 = dy1 - dy0
    dy1_2 = dy2 - dy1
    dy0_3 = dy1_2 - dy0_2
    return y_int_l[0] + dy0 * t + dy0_2 * t * (t - 1) / 2 + dy0_3 * t * (t - 1) * (t - 2) / 6


def interpolate_linear_spline(x, x_int_l, y_int_n, h):
    matrix = numpy.array([[1, 0, 0, 0, 0, 0],
                          [1, h, 0, 0, 0, 0],
                          [0, 0, 1, 0, 0, 0],
                          [0, 0, 1, h, 0, 0],
                          [0, 0, 0, 0, 1, 0],
                          [0, 0, 0, 0, 1, h]])
    vector = numpy.array([[y_int_n[0]],
                          [y_int_n[1]],
                          [y_int_n[1]],
                          [y_int_n[2]],
                          [y_int_n[2]],
                          [y_int_n[3]]])
    solved = numpy.linalg.solve(matrix, vector)
    return evaluate_spline(x, x_int_l, solved.ravel(), 1)


def interpolate_parabolic_spline(x, x_int_l, y_int_n, h):
    matrix = numpy.array([[1, 0, 0, 0, 0, 0, 0, 0, 0],
                          [1, h, h * h, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 1, h, h * h, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 1, h, h * h],
                          [1, h, h * h, -1, 0, 0, 0, 0, 0],
                          [0, 0, 0, 1, h, h * h, -1, 0, 0],
                          [0, 1, 2 * h, 0, -1, 0, 0, 0, 0],
                          [0, 0, 0, 0, 1, 2 * h, 0, -1, 0],
                          [0, 0, 1, 0, 0, 0, 0, 0, 0]])
    vector = numpy.array([[y_int_n[0]],
                          [y_int_n[1]],
                          [y_int_n[2]],
                          [y_int_n[3]],
                          [0],
                          [0],
                          [0],
                          [0],
                          [0]])
    solved = numpy.linalg.solve(matrix, vector)
    return evaluate_spline(x, x_int_l, solved.ravel(), 2)


def interpolate_cubic_spline(x, x_int_l, y_int_n, h):
    matrix = numpy.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [1, h, h * h, h * h * h, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 1, h, h * h, h * h * h, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 1, h, h * h, h * h * h],
                          [0, 1, 2 * h, 3 * h * h, 0, -1, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 1, 2 * h, 3 * h * h, 0, -1, 0, 0],
                          [0, 0, 2, 6 * h, 0, 0, -2, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 2, 6 * h, 0, 0, -2, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 6 * h],
                          [0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    vector = numpy.array([[y_int_n[0]],
                          [y_int_n[1]],
                          [y_int_n[1]],
                          [y_int_n[2]],
                          [y_int_n[2]],
                          [y_int_n[3]],
                          [0],
                          [0],
                          [0],
                          [0],
                          [0],
                          [0]])
    solved = numpy.linalg.solve(matrix, vector)
    return evaluate_spline(x, x_int_l, solved.ravel(), 3)


def evaluate_spline(x, x_int_l, coefficients, order):
    value = 0
    if x_int_l[0] <= x <= x_int_l[3]:
        i = 0
        while i <= order:
            n = 0 if x < x_int_l[1] else order + 1 if x_int_l[1] <= x < x_int_l[2] else 2 * (order + 1)
            n += i
            value += coefficients[n] * numpy.power(x - x_int_l[int(n / (order + 1))], n % (order + 1))
            i += 1
    return value


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
