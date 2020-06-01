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
    print("Найдём решение ОДУ с помощью метода Рунге-Кутта:")
    x = discretize(a, b, h)
    y_runge_kutta_h = runge_kutta_dsolve(x, y0, a, b, h)
    y_runge_kutta_2h = runge_kutta_dsolve(x, y0, a, b, 2 * h)
    print("x: " + str(x))
    print("y: " + str(y_runge_kutta_h))
    print("y с двойным шагом: " + str(y_runge_kutta_2h))
    plot(x, numpy.array([y_runge_kutta_h]), numpy.array(["Интегральная кривая, метод Рунге-Кутта (IV)"]))
    print("Найдём решение ОДУ с помощью метода Эйлера:")
    y_eiler_h = eiler_dsolve(x, y0, a, b, h)
    y_eiler_2h = eiler_dsolve(x, y0, a, b, 2 * h)
    print("x: " + str(x))
    print("y: " + str(y_eiler_h))
    print("y с двойным шагом: " + str(y_eiler_2h))
    plot(x, numpy.array([y_runge_kutta_h, y_eiler_h]),
         numpy.array(["Интегральная кривая, метод Рунге-Кутта (IV)", "Интегральная кривая, метод Эйлера"]))
    print("Найдём точное решение задачи Коши:")
    y = numpy.array([exact_solution(x[i]) for i in range(len(x))])
    print("y: " + str(y))
    plot(x, numpy.array([y, y_runge_kutta_h, y_eiler_h]), numpy.array(
        ["Точное решение", "Интегральная кривая, метод Рунге-Кутта (IV)", "Интегральная кривая, метод Эйлера"]))
    plot(x, numpy.array([y, y_runge_kutta_h, y_eiler_h]), numpy.array(
        ["Точное решение", "Интегральная кривая, метод Рунге-Кутта (IV)", "Интегральная кривая, метод Эйлера"]),
         numpy.array([1.995, 2, 0.666, 0.668]))

