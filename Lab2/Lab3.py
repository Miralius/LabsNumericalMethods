import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols, pi, diff, N

x = symbols('x')


def function(point):
    return numpy.math.pi * point / (10 + point / 2)


def define_step(func):
    return N(2 * ((eps / (abs(diff(func, x, 4)).subs(x, a))) ** (1 / 4)))


def discretize():
    n = round((b - a) / h)
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


def interpolate_with_the_newton_method(x, x_L, y_L):
    return x_L


def plot_function(x, func, name):
    xy = matplotlib.pyplot.subplot()
    xy.grid(which='major', color='k')
    xy.minorticks_on()
    xy.grid(which="minor", color='gray', linestyle=':')
    xy.plot(x, func, color="blue", label=name)
    xy.set_xlabel("x")
    xy.set_ylabel("y")
    xy.legend()
    matplotlib.pyplot.show()


def plot_functions(x, func1, func2, name1, name2):
    xy = matplotlib.pyplot.subplot()
    xy.grid(which='major', color='k')
    xy.minorticks_on()
    xy.grid(which="minor", color='gray', linestyle=':')
    xy.plot(x, func1, color="blue", label=name1)
    xy.plot(x, func2, color="red", label=name2)
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
    h = define_step(pi * x / (10 + x / 2))
    print("Шаг = " + str(h))
    number, x, y = discretize()
    plot_function(x, y, "f(x) = pi * x / (10 + x / 2)")
    print("x = " + str(x))
    number_x = int(input("Выберите начало отрезка x для интерполяции полиномами Лагранжа: "))
    x_L = numpy.array([x[number_x + i] for i in range(4)])
    y_interpolated_by_lagrange = numpy.array([interpolate_with_the_lagrange_polynomial(x_L[number_x + i], x_L) for i in range(len(x_L))])
    a = x_L[0]
    b = x_L[3]
    h = (x_L[3] - x_L[0]) / number
    number, x, y = discretize()
    y_L = numpy.array([interpolate_with_the_lagrange_polynomial(x[i], x_L) for i in range(len(x))])
    plot_functions(x, y, y_L, "f(x)", "interpolated function by Lagrange polynomial")
    plot_function(x, numpy.array([abs(y_L[i] - y[i]) for i in range(len(x))]), "absolute error")

    print("Второе задание ---------")
    y_interpolated_by_newton = numpy.array([interpolate_with_the_newton_method(x_L[i], y_interpolated_by_lagrange[i]) for i in range(len(x_L))])
