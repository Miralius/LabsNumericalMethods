import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols, pi, diff, N

x = symbols('x')


def function(point):
    return numpy.math.pi * point / (10 + point / 2)


def define_step(func):
    return N(2 * ((eps / (abs(diff(func, x, 4)).subs(x, a))) ** (1 / 4)))


def discretize(left, right, h):
    n = round((right - left) / h)
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


def plot_function(x, func):
    xy = matplotlib.pyplot.subplot()
    xy.grid(which='major', color='k')
    xy.minorticks_on()
    xy.grid(which="minor", color='gray', linestyle=':')
    xy.plot(x, func, color="blue", label="f1(x)")
    xy.set_xlabel("x")
    xy.set_ylabel("y")
    xy.legend()
    matplotlib.pyplot.show()


def plot_functions(ox, func1, func2):
    xy = matplotlib.pyplot.subplot()
    xy.grid(which='major', color='k')
    xy.minorticks_on()
    xy.grid(which="minor", color='gray', linestyle=':')
    xy.plot(ox, func1, color="blue", label="f1(x)")
    xy.plot(ox, func2, color="red", label="f2(x)")
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
    step = define_step(pi * x / (10 + x / 2))
    print("Шаг = " + str(step))
    number, x, y = discretize(a, b, step)
    plot_function(x, y)
    print("x = " + str(x))
    number_x = int(input("Выберите начало отрезка x для интерполяции полиномами Лагранжа: "))
    number_L, x_L, y_L = discretize(x[number_x], x[number_x + 3], (x[number_x + 3] - x[number_x]) / number)
    y_interpolated_by_lagrange = numpy.array([interpolate_with_the_lagrange_polynomial(x_L[j], numpy.array([x[number_x + i] for i in range(4)])) for j in range(len(x_L))])
    plot_functions(x_L, y_L, y_interpolated_by_lagrange)
    plot_function(x_L, numpy.array([abs(y_L[i] - y_interpolated_by_lagrange[i]) for i in range(len(x_L))]))
