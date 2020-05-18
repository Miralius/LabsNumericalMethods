import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols
from sympy import S, sqrt
from sympy import solveset
from sympy import diff
from sympy import bessely
from sympy import besselj

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
    xy.plot(x_num, func_numpy, color="black", label="f(x)")
    xy.set_xlabel("x")
    xy.set_ylabel("y")
    xy.legend()
    matplotlib.pyplot.show()


def find_solve_with_secant_method(func, a, b, eps):
    root = 0
    count = 0
    exact_solution = 0
    i = 0
    # solutions = solveset(func)
    # while i < len(solutions):
    #    if a < solutions.args[i] < b:
    #        exact_solution = solutions.args[i]
    #    i += 1
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
    print("МЕТОД ПОЛОВИННОГО ДЕЛЕНИЯ")
    print("Найденный корень x* = " + str(root))
    print("Невязка f(x*) = " + str(func.subs(x, root)))
    print("Кол-во итераций = " + str(count))
    #print("Абсолютная погрешность = " + str(abs(root - exact_solution)))
    #print("Относительная погрешность = " + str(abs(root - exact_solution) / abs(root) * 100) + "%")


def find_solve_with_hybrid_method(func, a, b, x0, eps):
    x_k = x0
    count = 0
    # exact_solution = 0
    i = 0
    solutions = solveset(func)
    # while i < len(solutions):
    #     if a < solutions.args[i] < b:
    #        exact_solution = solutions.args[i]
    #    i += 1
    while True:
        x_k1 = x_k - (func.subs(x, x_k)) / ((diff(func, x)).subs(x, x_k))
        while True:
            count += 1
            if abs(func.subs(x, x_k1)) < abs(func.subs(x, x_k)):
                break
            else:
                x_k_temp = x_k
                x_k = x_k1
                x_k1 = (x_k1 + x_k_temp) / 2
        if abs(x_k - x_k1) <= eps:
            break
        x_k = x_k1
    print("ГИБРИДНЫЙ МЕТОД")
    print("Найденный корень x* = " + str(x_k1))
    print("Невязка f(x*) = " + str(func.subs(x, x_k1)))
    print("Кол-во итераций = " + str(count))
    #print("Абсолютная погрешность = " + str(abs(x_k1 - exact_solution)))
    #print("Относительная погрешность = " + str(abs(x_k1 - exact_solution) / abs(x_k1) * 100) + "%")


if __name__ == "__main__":
    alpha = float(input("\nВведите отношение масс шарика и струны: "))
    # plot_function(besselj(0, 2 * sqrt(alpha + 1) * x) * (sqrt(alpha) * x * bessely(0, 2 * sqrt(alpha) * x) - bessely(1, 2 * sqrt(alpha) * x)) - bessely(0, 2 * sqrt(alpha + 1) * x) * (sqrt(alpha) * x * besselj(0, 2 * sqrt(alpha) * x) - besselj(1, 2 * sqrt(alpha) * x)), 0, 8)
    while True:
        left = float(input("\nВведите левую границу отрезка: "))
        right = float(input("Введите правую границу отрезка: "))
        init_x = float(input("Введите начальное приближение: "))
        print("Для метода половинного деления задан отрезок [" + str(left) + "; " + str(right) + "]")
        print("Для гибридного метода задано начальное приближение: " + str(init_x))
        if init_x == left or init_x == right:
            print("Начальное приближение совпадает с одним из концов отрезка")
        else:
            print("Начальное приближение не совпадает ни с одним из концов отрезка")
        find_solve_with_secant_method(besselj(0, 2 * sqrt(alpha + 1) * x) * (sqrt(alpha) * x * bessely(0, 2 * sqrt(alpha) * x) - bessely(1, 2 * sqrt(alpha) * x)) - bessely(0, 2 * sqrt(alpha + 1) * x) * (sqrt(alpha) * x * besselj(0, 2 * sqrt(alpha) * x) - besselj(1, 2 * sqrt(alpha) * x)), left, right, 10 ** (-3))
        # find_solve_with_hybrid_method(besselj(0, 2 * sqrt(alpha + 1) * x) * (sqrt(alpha) * x * bessely(0, 2 * sqrt(alpha) * x) - bessely(1, 2 * sqrt(alpha) * x)) - bessely(0, 2 * sqrt(alpha + 1) * x) * (sqrt(alpha) * x * besselj(0, 2 * sqrt(alpha) * x) - besselj(1, 2 * sqrt(alpha) * x)), left, right, init_x, 10 ** (-3))
        progress = int(input("Для продолжения нажмите введите 1, для завершения  0: "))
        if progress == 0:
            break
