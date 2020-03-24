from math import sqrt
import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols, diff, S, Eq, cos, sin, solve

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


def find_solve_system_with_newton_method(f1, f2, x0, y0, eps):
    count = 0
    x_k = x0
    y_k = y0
    x_k1 = (x_k - diff(f2, y).subs(x, x_k).subs(y, y_k) * f1.subs(x, x_k).subs(y, y_k) + diff(f1, y).subs(x, x_k).subs(y, y_k) * f2.subs(x, x_k).subs(y, y_k)).n(16)
    y_k1 = (y_k - diff(f2, x).subs(x, x_k).subs(y, y_k) * f1.subs(x, x_k).subs(y, y_k) + diff(f1, x).subs(x, x_k).subs(y, y_k) * f2.subs(x, x_k).subs(y, y_k)).n(16)
    while sqrt((x_k1 - x_k) ** 2 + (y_k1 - y_k) ** 2) > eps:
        count += 1
        x_k_1 = x_k
        y_k_1 = y_k
        x_k = x_k1
        y_k = y_k1
        x_k1 = (x_k - diff(f2, y).subs(x, x_k).subs(y, y_k) * f1.subs(x, x_k).subs(y, y_k) + diff(f1, y).subs(x,x_k).subs(y, y_k) * f2.subs(x, x_k).subs(y, y_k)).n(16)
        y_k1 = (y_k - diff(f2, x).subs(x, x_k).subs(y, y_k) * f1.subs(x, x_k).subs(y, y_k) + diff(f1, x).subs(x, x_k).subs(y, y_k) * f2.subs(x, x_k).subs(y, y_k)).n(16)
        if sqrt((x_k1 - x_k) ** 2 + (y_k1 - y_k) ** 2) > sqrt((x_k - x_k_1) ** 2 + (y_k - y_k_1) ** 2):
            print("Метод расходится! :(")
            break
    print("МЕТОД НЬЮТОНА")
    print("x* = " + str(x_k1))
    print("y* = " + str(y_k1))
    print("Количество итераций = " + str(count))
    print("f1(x*) = " + str(f1.subs(x, x_k1).subs(y, y_k1)))
    print("f2(x*) = " + str(f2.subs(x, x_k1).subs(y, y_k1)))
    print("Погрешность вычислений для f1 = " + str(f1.subs(x, x_k1 - x_k).subs(y, y_k1 - y_k)))
    print("Погрешность вычислений для f2 = " + str(f2.subs(x, x_k1 - x_k).subs(y, y_k1 - y_k)))

def find_solve_system_with_modified_newton_method(f1, f2, x0, y0, eps):
    count = 0
    x_k = x0
    y_k = y0
    df1_dx = (diff(f1, x).subs(x, x_k).subs(y, y_k)).n(16)
    df1_dy = (diff(f1, y).subs(x, x_k).subs(y, y_k)).n(16)
    df2_dx = (diff(f2, x).subs(x, x_k).subs(y, y_k)).n(16)
    df2_dy = (diff(f2, y).subs(x, x_k).subs(y, y_k)).n(16)
    x_k1 = (x_k - df2_dy * f1.subs(x, x_k).subs(y, y_k) + df1_dy * f2.subs(x, x_k).subs(y, y_k)).n(16)
    y_k1 = (y_k - df2_dx * f1.subs(x, x_k).subs(y, y_k) + df1_dx * f2.subs(x, x_k).subs(y, y_k)).n(16)
    while sqrt((x_k1 - x_k) ** 2 + (y_k1 - y_k) ** 2) > eps:
        count += 1
        x_k_1 = x_k
        y_k_1 = y_k
        x_k = x_k1
        y_k = y_k1
        x_k1 = (x_k - df2_dy * f1.subs(x, x_k).subs(y, y_k) + df1_dy * f2.subs(x, x_k).subs(y, y_k)).n(16)
        y_k1 = (y_k - df2_dx * f1.subs(x, x_k).subs(y, y_k) + df1_dx * f2.subs(x, x_k).subs(y, y_k)).n(16)
        if sqrt((x_k1 - x_k) ** 2 + (y_k1 - y_k) ** 2) > sqrt((x_k - x_k_1) ** 2 + (y_k - y_k_1) ** 2):
            print("Метод расходится! :(")
            break
    print("МОДИФИЦИРОВАННЫЙ МЕТОД НЬЮТОНА")
    print("x* = " + str(x_k1))
    print("y* = " + str(y_k1))
    print("Количество итераций = " + str(count))
    print("f1(x*) = " + str(f1.subs(x, x_k1).subs(y, y_k1)))
    print("f2(x*) = " + str(f2.subs(x, x_k1).subs(y, y_k1)))
    print("Погрешность вычислений для f1 = " + str(f1.subs(x, x_k1 - x_k).subs(y, y_k1 - y_k)))
    print("Погрешность вычислений для f2 = " + str(f2.subs(x, x_k1 - x_k).subs(y, y_k1 - y_k)))

plot_functions(Eq(cos(x + 0.5) - y, 2), Eq(sin(y) - 2 * x, 1), -10, 10)
plot_functions(Eq(cos(x + 0.5) - y, 2), Eq(sin(y) - 2 * x, 1), -4, 0)
find_solve_system_with_newton_method(cos(x + 0.5) - y - S(2), sin(y) - 2 * x - S(1), -1, -1, 10 ** (-3))
find_solve_system_with_modified_newton_method(cos(x + 0.5) - y - S(2), sin(y) - 2 * x - S(1), -1, -1, 10 ** (-3))