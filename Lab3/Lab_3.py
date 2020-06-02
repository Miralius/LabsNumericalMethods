from Lab3.functions import *

if __name__ == "__main__":
    print("Лабораторная работа №3, вариант #14")
    a = 1
    b = 4
    eps = 10 ** (-3)
    h = define_step(a, b, eps)
    print("Исходная функция: f(x)=exp(-sqrt(x))")
    print("Границы интегрирования a и b: " + str(a) + " и " + str(b))
    print("Заданная точность вычислений: " + str(eps))
    print("1) Шаг интегрирования h: " + str(h))
    print("2) Значение интеграла по формуле трапеций c шагом 2h: " + str(trapezes_integrate(a, b, 2 * h)))
    print("Значение интеграла по формуле трапеций c шагом h: " + str(trapezes_integrate(a, b, h)))
    print("Уточненная погрешность при вычислении по формуле трапеций: " + str(abs(trapezes_integrate(a, b, 2 * h) -
                                                                                  trapezes_integrate(a, b, h)) / 3))
    print("3) Значение интеграла по формуле Симпсона c шагом 2h: " + str(simpson_integrate(a, b, 2 * h)))
    print("Значение интеграла по формуле Симпсона c шагом h: " + str(simpson_integrate(a, b, h)))
    print("Уточненная погрешность при вычислении по формуле Симпсона: " + str(abs(simpson_integrate(a, b, 2 * h) -
                                                                                  simpson_integrate(a, b, h)) / 15))
    print("4) Точное значение интеграла, посчитанного по формуле Ньютона-Лейбница: " +
          str(newton_leibniz_integrate(a, b)))
    print("Точное значение интеграла, посчитанного по формуле Ньютона-Лейбница (десятичной дробью): " +
          str(N(newton_leibniz_integrate(a, b))))
    print("Сравнение точного значения со значением, посчитанным методом трапеций, шаг 2h: " +
          str(abs(N(newton_leibniz_integrate(a, b) - trapezes_integrate(a, b, 2 * h)))))
    print("Сравнение точного значения со значением, посчитанным методом трапеций, шаг h: " +
          str(abs(N(newton_leibniz_integrate(a, b) - trapezes_integrate(a, b, h)))))
    print("Сравнение точного значения со значением, посчитанным методом Симпсона, шаг 2h: " +
          str(abs(N(newton_leibniz_integrate(a, b) - simpson_integrate(a, b, 2 * h)))))
    print("Сравнение точного значения со значением, посчитанным методом Симпсона, шаг h (самый точный результат!): " +
          str(abs(N(newton_leibniz_integrate(a, b) - simpson_integrate(a, b, h)))))
