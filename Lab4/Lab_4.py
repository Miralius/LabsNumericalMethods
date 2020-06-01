from Lab4.functions import *

if __name__ == "__main__":
    print("Лабораторная работа №4, вариант №14")
    print("Решение задачи Коши: y' + y = 0.5xy², y(0) = 2, a = 0, b = 2")
    h = 0.1
    a = 0
    b = 2
    y0 = 2
    x0 = 0
    h = define_step(h, a, b, y0, x0)
    print("Оптимальный шаг h = " + str(h))
    x = discretize(a, b, h)
    y_runge_kutta_h = runge_kutta_dsolve(x, y0, a, b, h)
    y_runge_kutta_2h = runge_kutta_dsolve(x, y0, a, b, 2 * h)
    print("Точка остановыыыыыыыыы")

