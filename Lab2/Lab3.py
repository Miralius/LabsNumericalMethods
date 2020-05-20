import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols, pi, diff, N

xx = symbols('x')


def function(point):
    return numpy.math.pi * point / (10 + point / 2)


def define_step(func):
    return float(N(2 * ((eps / (abs(diff(func, xx, 4)).subs(xx, a))) ** (1 / 4))))


def discretize():
    n = int(numpy.ceil((b - a) / h))
    return n, numpy.array([a + i * h for i in range(n)]), numpy.array([function(a + i * h) for i in range(n)])


def interpolate_with_the_lagrange_polynomial(point, interval):
    i = 0
    interpolated_function = 0
    while i < 4:
        j = 0
        multiplier = 1
        while j < 4:
            if i != j:
                multiplier *= (point - interval[j]) / (interval[i] - interval[j])
            j += 1
        interpolated_function += function(interval[i]) * multiplier
        i += 1
    return interpolated_function


def interpolate_with_the_newton_method(point):
    t = (point - x_L[0]) / h
    dy0 = y_interpolated_by_lagrange[1] - y_interpolated_by_lagrange[0]
    dy1 = y_interpolated_by_lagrange[2] - y_interpolated_by_lagrange[1]
    dy2 = y_interpolated_by_lagrange[3] - y_interpolated_by_lagrange[2]
    dy0_2 = dy1 - dy0
    dy1_2 = dy2 - dy1
    dy0_3 = dy1_2 - dy0_2
    return y_interpolated_by_lagrange[0] + dy0 * t + dy0_2 * t * (t - 1) / 2 + dy0_3 * t * (t - 1) * (t - 2) / 6


def interpolate_linear_spline(point, y_n):
    matrix = numpy.array([[1, 0, 0, 0, 0, 0],
                          [1, h, 0, 0, 0, 0],
                          [0, 0, 1, 0, 0, 0],
                          [0, 0, 1, h, 0, 0],
                          [0, 0, 0, 0, 1, 0],
                          [0, 0, 0, 0, 1, h]])
    vector = numpy.array([[y_n[0]],
                          [y_n[1]],
                          [y_n[1]],
                          [y_n[2]],
                          [y_n[2]],
                          [y_n[3]]])
    solved = numpy.linalg.solve(matrix, vector)
    return evaluate_spline(point, solved.ravel(), 1)


def interpolate_parabolic_spline(point, y_n):
    matrix = numpy.array([[1, 0, 0, 0, 0, 0, 0, 0, 0],
                          [1, h, h * h, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 1, h, h * h, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 1, h, h * h],
                          [1, h, h * h, -1, 0, 0, 0, 0, 0],
                          [0, 0, 0, 1, h, h * h, -1, 0, 0],
                          [0, 1, 2 * h, 0, -1, 0, 0, 0, 0],
                          [0, 0, 0, 0, 1, 2 * h, 0, -1, 0],
                          [0, 0, 1, 0, 0, 0, 0, 0, 0]])
    vector = numpy.array([[y_n[0]],
                          [y_n[1]],
                          [y_n[2]],
                          [y_n[3]],
                          [0],
                          [0],
                          [0],
                          [0],
                          [0]])
    solved = numpy.linalg.solve(matrix, vector)
    return evaluate_spline(point, solved.ravel(), 2)


def interpolate_cubic_spline(point, y_n):
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
    vector = numpy.array([[y_n[0]],
                          [y_n[1]],
                          [y_n[1]],
                          [y_n[2]],
                          [y_n[2]],
                          [y_n[3]],
                          [0],
                          [0],
                          [0],
                          [0],
                          [0],
                          [0]])
    solved = numpy.linalg.solve(matrix, vector)
    return evaluate_spline(point, solved.ravel(), 3)


def evaluate_spline(point, coefficients, order):
    value = 0
    if x_L[0] < point < x_L[3]:
        i = 0
        while i <= order:
            if point < x_L[1]:
                n = 0
            elif x_L[1] <= point < x_L[2]:
                n = order + 1
            else:
                n = 2 * (order + 1)
            n += i
            degree = n % (order + 1)
            index = int(n / (order + 1))
            base = point - x_L[index]
            value += coefficients[n] * numpy.power(base, degree)
            i += 1
    return value


def plot_function(points, func, name):
    xy = matplotlib.pyplot.subplot()
    xy.grid(which='major', color='k')
    xy.minorticks_on()
    xy.grid(which="minor", color='gray', linestyle=':')
    xy.plot(points, func, color="deeppink", label=name)
    xy.set_xlabel("x")
    xy.set_ylabel("y")
    xy.legend()
    matplotlib.pyplot.show()


def plot_functions(points, func1, func2, name1, name2):
    xy = matplotlib.pyplot.subplot()
    xy.grid(which='major', color='k')
    xy.minorticks_on()
    xy.grid(which="minor", color='gray', linestyle=':')
    xy.plot(points, func1, color="blue", label=name1)
    xy.plot(points, func2, color="deeppink", label=name2)
    xy.set_xlabel("x")
    xy.set_ylabel("y")
    xy.legend()
    matplotlib.pyplot.show()


def plot_three_functions(points, func1, func2, func3, name1, name2, name3):
    xy = matplotlib.pyplot.subplot()
    xy.grid(which='major', color='k')
    xy.minorticks_on()
    xy.grid(which="minor", color='gray', linestyle=':')
    xy.plot(points, func1, color="blue", label=name1)
    xy.plot(points, func2, color="deeppink", label=name2)
    xy.plot(points, func3, color="green", label=name3)
    xy.set_xlabel("x")
    xy.set_ylabel("y")
    xy.legend()
    matplotlib.pyplot.show()


if __name__ == "__main__":
    print("Лабораторная работа №3, вариант #14")
    a = 0
    b = 10
    eps = 10 ** (-4)
    print("a = " + str(a))
    print("b = " + str(b))
    print("Точность = " + str(eps))

    print("Первое задание ---------")
    h = define_step(pi * xx / (10 + xx / 2))
    print("Шаг = " + str(h))
    number, x, y = discretize()
    plot_function(x, y, "f(x) = pi * x / (10 + x / 2)")
    print("x = " + str(x))
    number_x = int(input("Выберите начало отрезка x для интерполяции полиномами Лагранжа: "))
    x_L = numpy.array([x[number_x + i] for i in range(4)])
    y_interpolated_by_lagrange = \
        numpy.array([interpolate_with_the_lagrange_polynomial(x_L[i], x_L) for i in range(len(x_L))])
    a = x_L[0]
    b = x_L[3]
    h = (x_L[3] - x_L[0]) / number
    number, x, y = discretize()
    y_L = numpy.array([interpolate_with_the_lagrange_polynomial(x[i], x_L) for i in range(len(x))])
    plot_functions(x, y, y_L, "f(x)", "interpolated function by Lagrange polynomial")
    plot_function(x, numpy.array([abs(y_L[i] - y[i]) for i in range(len(x))]), "absolute error")

    print("Второе задание ---------")
    a = 0
    b = 10
    h = define_step(pi * xx / (10 + xx / 2))
    y_N = numpy.array([interpolate_with_the_newton_method(x_L[i]) for i in range(len(x_L))])
    y_interpolated_by_newton = numpy.array([interpolate_with_the_newton_method(x[i]) for i in range(len(x))])
    plot_three_functions(x, y, y_L, y_interpolated_by_newton, "f(x)", "Lagrange polynomial interpolation",
                         "Newton polynomial interpolation")
    plot_function(x, numpy.array([abs(y_interpolated_by_newton[i] - y[i]) for i in range(len(x))]),
                  "absolute error f(x) & f_N(x)")
    plot_function(x, numpy.array([abs(y_interpolated_by_newton[i] - y_L[i]) for i in range(len(x))]),
                  "absolute error f_L(x) & f_N(x)")

    print("Третье задание ---------")
    a = x_L[0]
    b = x_L[3]
    h = (x_L[3] - x_L[0]) / 100
    number, x, y = discretize()
    y_interpolated_linear_spline = numpy.array([interpolate_linear_spline(x[i], y_N) for i in range(len(x))])
    y_interpolated_parabolic_spline = numpy.array([interpolate_parabolic_spline(x[i], y_N) for i in range(len(x))])
    y_interpolated_cubic_spline = numpy.array([interpolate_cubic_spline(x[i], y_N) for i in range(len(x))])
    plot_three_functions(x, y_interpolated_linear_spline, y_interpolated_parabolic_spline, y_interpolated_cubic_spline,
                         "linear spline", "parabolic spline", "cubic spline")
