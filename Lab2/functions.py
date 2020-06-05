import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols, pi, diff, N, solve, Eq


def function(x):
    return numpy.math.pi * x / (10 + x / 2)


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
    x = symbols('x')
    func = pi * x / (10 + x / 2)
    return float(N(2 * ((eps / m_n_plus_one(3, a, b, func)) ** (1 / 4))))


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
    while i < len(x_int_l) - 1:
        ak[i] = solve(Eq(spline, y).subs(x0, x_int_l[i]).subs(x, x_int_l[i]).subs(y, y_int_n[i]), a)[0]
        k = numpy.append(k, ak[i])
        bk[i] = solve(
            Eq(spline, y).subs(a, ak[i]).subs(x0, x_int_l[i]).subs(x, x_int_l[i + 1]).subs(y, y_int_n[i + 1]), b)[0]
        k = numpy.append(k, bk[i])
        i += 1
    return evaluate_spline(point, x_int_l, k, 1)


def interpolate_parabolic_spline(point, x_int_l, y_int_n):
    a, b, c, x0, x, y = symbols('a b c x0, x y')
    spline = a + b * (x - x0) + c * (x - x0) ** 2
    ak = numpy.zeros(len(x_int_l) - 1)
    bk = numpy.zeros(len(x_int_l) - 1)
    ck = numpy.zeros(len(x_int_l) - 1)
    k = numpy.zeros(3)
    k[0] = ak[0] = solve(Eq(spline, y).subs(x0, x_int_l[0]).subs(x, x_int_l[0]).subs(y, y_int_n[0]), a)[0]
    k[2] = ck[0] = solve(Eq(diff(spline, x, 2), 0).subs(x0, x_int_l[0]).subs(x, x_int_l[0]), c)[0]
    k[1] = bk[0] = solve(
        Eq(spline, y).subs(a, ak[0]).subs(c, ck[0]).subs(x0, x_int_l[0]).subs(x, x_int_l[1]).subs(y, y_int_n[1]), b)[0]
    i = 1
    while i < len(x_int_l) - 1:
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


def interpolate_cubic_spline(point, x_int_l, y_int_n):
    a, b, c, d, x0, x, y = symbols('a b c d x0, x y')
    spline = a + b * (x - x0) + c * (x - x0) ** 2 + d * (x - x0) ** 3
    splines = numpy.array([])
    a_s = numpy.array([])
    b_s = numpy.array([])
    c_s = numpy.array([])
    d_s = numpy.array([])
    i = 1
    while i < len(x_int_l):
        a_s = numpy.append(a_s, symbols('a' + str(i)))
        b_s = numpy.append(b_s, symbols('b' + str(i)))
        c_s = numpy.append(c_s, symbols('c' + str(i)))
        d_s = numpy.append(d_s, symbols('d' + str(i)))
        splines = numpy.append(
            splines, spline.subs(a, a_s[i - 1]).subs(b, b_s[i - 1]).subs(c, c_s[i - 1]).subs(d, d_s[i - 1]))
        i += 1
    ak = numpy.array([])
    bk = numpy.array([])
    ck = numpy.array([])
    dk = numpy.array([])
    ck = numpy.append(ck, solve(Eq(diff(spline, x, 2), 0).subs(x0, x_int_l[0]).subs(x, x_int_l[0]), c)[0])
    i = 0
    while i < len(x_int_l) - 1:
        ak = numpy.append(ak,
                          solve(Eq(splines[i], y).subs(x0, x_int_l[i]).subs(x, x_int_l[i]).subs(y, y_int_n[i]), a_s[i])[
                              0])
        bk = numpy.append(bk, solve(
            Eq(splines[i], y).subs(c_s[i], ck[i]).subs(a_s[i], ak[i]).subs(x0, x_int_l[i]).subs(x, x_int_l[i + 1]).subs(
                y, y_int_n[i + 1]), b_s[i])[0])
        if i != len(x_int_l) - 2:
            ck = numpy.append(ck, solve(
                Eq(diff(splines[i], x, 2).subs(c_s[i], ck[i]).subs(x0, x_int_l[i]).subs(x, x_int_l[i + 1]),
                   diff(splines[i + 1], x, 2).subs(x0, x_int_l[i + 1]).subs(x, x_int_l[i + 1])), c_s[i + 1])[0])
            dk = numpy.append(dk, solve(Eq(
                diff(splines[i], x).subs(b_s[i], bk[i]).subs(c_s[i], ck[i]).subs(x0, x_int_l[i]).subs(x,
                                                                                                      x_int_l[i + 1]),
                diff(splines[i + 1], x).subs(x0, x_int_l[i + 1]).subs(x, x_int_l[i + 1])), d_s[i])[0])
        i += 1
    dk = numpy.append(dk, solve(Eq(
        diff(splines[len(x_int_l) - 2], x, 2).subs(c_s[len(x_int_l) - 2], ck[len(x_int_l) - 2]).subs(x0, x_int_l[
            len(x_int_l) - 2]).subs(x, x_int_l[len(x_int_l) - 1]), 0), d_s[len(x_int_l) - 2])[0])
    i = 0
    while i < len(x_int_l) - 2:
        dk[i] = solve(Eq(dk[i].subs(b_s[i + 1], bk[i + 1]), d_s[i]), d_s[i])[0]
        i += 1
    i = 0
    while i < len(x_int_l) - 2:
        j = i + 1
        while j < len(x_int_l) - 1:
            dk[j] = solve(Eq(dk[j].subs(d_s[i], dk[i]), d_s[j]), d_s[j])[0]
            j += 1
        i += 1
    while i >= 0:
        j = len(x_int_l) - 2
        while j >= 0:
            bk[j] = solve(Eq(bk[j].subs(d_s[i], dk[i]), b_s[j]), b_s[j])[0]
            ck[j] = solve(Eq(ck[j].subs(d_s[i], dk[i]), c_s[j]), c_s[j])[0]
            dk[j] = solve(Eq(dk[j].subs(d_s[i], dk[i]), d_s[j]), d_s[j])[0]
            j -= 1
        i -= 1
    i = 0
    k = numpy.array([])
    while i < len(x_int_l) - 1:
        k = numpy.append(k, ak[i])
        k = numpy.append(k, bk[i])
        k = numpy.append(k, ck[i])
        k = numpy.append(k, dk[i])
        i += 1
    return evaluate_spline(point, x_int_l, k, 3)


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
