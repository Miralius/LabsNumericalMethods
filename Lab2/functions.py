import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols, pi, diff, N, solve, Eq


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


def interpolate_linear_spline(point, x_int_l, y_int_n):
    a, b, x0, x, y = symbols('a b x0, x y')
    spline = a + b * (x - x0)
    ak = numpy.zeros(len(x_int_l) - 1)
    bk = numpy.zeros(len(x_int_l) - 1)
    k = numpy.array([])
    i = 0
    while i < len(ak):
        ak[i] = solve(
            Eq(spline, y).subs(x0, x_int_l[i]).subs(x, x_int_l[i]).subs(y, y_int_n[i]), a)[0]
        k = numpy.append(k, ak[i])
        bk[i] = solve(
            Eq(spline, y).subs(a, ak[i]).subs(x0, x_int_l[i]).subs(x, x_int_l[i + 1]).subs(y, y_int_n[i + 1]), b)[0]
        k = numpy.append(k, bk[i])
        i += 1
    return evaluate_spline(point, x_int_l, k, 1)


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


def interpolate_parabolic_spline(point, x_int_l, y_int_n):
    a, b, c, x0, x, y = symbols('a b c x0, x y')
    spline = a + b * (x - x0) + c * (x - x0) ** 2
    ak = numpy.zeros(len(x_int_l) - 1)
    bk = numpy.zeros(len(x_int_l) - 1)
    ck = numpy.zeros(len(x_int_l) - 1)
    k = numpy.zeros(3)
    k[0] = ak[0] = solve(Eq(spline, y).subs(x0, x_int_l[0]).subs(x, x_int_l[0]).subs(y, y_int_n[0]), a)[0]
    k[2] = ck[0] = solve(Eq(diff(spline, x, 2), 0).subs(x0, x_int_l[0]).subs(x, x_int_l[0]).subs(y, y_int_n[0]), c)[0]
    k[1] = bk[0] = solve(
        Eq(spline, y).subs(a, ak[0]).subs(c, ck[0]).subs(x0, x_int_l[0]).subs(x, x_int_l[1]).subs(y, y_int_n[1]), b)[0]
    i = 1
    while i < len(ak):
        ak[i] = solve(Eq(spline.subs(a, ak[i - 1]).subs(b, bk[i - 1]).subs(c, ck[i - 1]).subs(x0, x_int_l[i - 1]).
                         subs(x, x_int_l[i]), spline.subs(x0, x_int_l[i]).subs(x, x_int_l[i])), a)[0]
        k = numpy.append(k, ak[i])
        bk[i] = solve(Eq(diff(spline, x).subs(b, bk[i - 1]).subs(c, ck[i - 1]).subs(x0, x_int_l[i - 1]).
                         subs(x, x_int_l[i]), diff(spline, x).subs(x0, x_int_l[i]).subs(x, x_int_l[i])), b)[0]
        k = numpy.append(k, bk[i])
        ck[i] = solve(Eq(spline, y).subs(a, ak[i]).subs(b, bk[i]).subs(x0, x_int_l[i]).subs(x, x_int_l[i + 1]).subs(
            y, y_int_n[i + 1]), c)[0]
        k = numpy.append(k, ck[i])
        i += 1
    return evaluate_spline(point, x_int_l, k, 2)


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


def define_index_of_spline(x, x_int_l, order):
    j = 1
    while x >= x_int_l[j]:
        if j + 1 == len(x_int_l):
            break
        else:
            j += 1
    return (j - 1) * (order + 1)


def evaluate_spline(x, x_int_l, coefficients, order):
    value = 0
    if x_int_l[0] <= x <= x_int_l[len(x_int_l) - 1]:
        i = 0
        while i <= order:
            n = define_index_of_spline(x, x_int_l, order)
            n += i
            value += coefficients[n] * numpy.power(x - x_int_l[n // (order + 1)], n % (order + 1))
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
