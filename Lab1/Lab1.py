import matplotlib.pyplot
import matplotlib.ticker
import numpy
from sympy import symbols
from sympy import S
from sympy import solveset
from sympy import diff

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
    solutions = solveset(func)
    while i < len(solutions):
        if a < solutions.args[i] < b:
            exact_solution = solutions.args[i]
        i += 1
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
    print("Абсолютная погрешность = " + str(abs(root - exact_solution)))
    print("Относительная погрешность = " + str(abs(root - exact_solution) / abs(root) * 100) + "%")


def find_solve_with_hybrid_method(func, a, b, x0, eps):
    x_k = x0
    count = 0
    exact_solution = 0
    i = 0
    solutions = solveset(func)
    while i < len(solutions):
        if a < solutions.args[i] < b:
            exact_solution = solutions.args[i]
        i += 1
    while True:
        x_k1 = x0 - (func.subs(x, x_k)) / ((diff(func, x)).subs(x, x_k))
        while True:
            count += 1
            if abs(func.subs(x, x_k1)) < abs(func.subs(x, x_k)):
                root = x_k1
                x_k = (x_k1 + x_k) / 2
                break
            else:
                x_k_temp = x_k
                x_k = x_k1
                x_k1 = (x_k1 + x_k_temp) / 2
        if abs(root - x_k) <= eps:
            break
    print("ГИБРИДНЫЙ МЕТОД")
    print("Найденный корень x* = " + str(root))
    print("Невязка f(x*) = " + str(func.subs(x, root)))
    print("Кол-во итераций = " + str(count))
    print("Абсолютная погрешность = " + str(abs(root - exact_solution)))
    print("Относительная погрешность = " + str(abs(root - exact_solution) / abs(root) * 100) + "%")


plot_function(x ** 3 + 8.5 * x ** 2 + 21.8 * x + S(15.6), -10, 10)
plot_function(x ** 3 + 8.5 * x ** 2 + 21.8 * x + S(15.6), -8, 0)
plot_function(x ** 3 + 8.5 * x ** 2 + 21.8 * x + S(15.6), -4.8, 0)
find_solve_with_secant_method(x ** 3 + 8.5 * x ** 2 + 21.8 * x + S(15.6), -4.8, -3.84, 10 ** (-3))
find_solve_with_secant_method(x ** 3 + 8.5 * x ** 2 + 21.8 * x + S(15.6), -3.84, -2.88, 10 ** (-3))
find_solve_with_secant_method(x ** 3 + 8.5 * x ** 2 + 21.8 * x + S(15.6), -1.92, -0.96, 10 ** (-3))
find_solve_with_hybrid_method(x ** 3 + 8.5 * x ** 2 + 21.8 * x + S(15.6), -4.8, -3.84, -4, 10 ** (-3))
find_solve_with_hybrid_method(x ** 3 + 8.5 * x ** 2 + 21.8 * x + S(15.6), -3.84, -2.88, -3, 10 ** (-3))
find_solve_with_hybrid_method(x ** 3 + 8.5 * x ** 2 + 21.8 * x + S(15.6), -1.92, -0.96, -1, 10 ** (-3))