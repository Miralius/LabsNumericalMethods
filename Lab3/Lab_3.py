from Lab3.functions import *

if __name__ == "__main__":
    print("Лабораторная работа №3, вариант #14")
    a = 1
    b = 4
    eps = 10 ** (-3)
    h = define_step(a, b, eps)
    print("Исходная функция: f(x)=exp(-sqrt(x))")
    print("Границы интегрирования a и b: " + str(a) + " и " + str(b))
    print("1) Шаг интегрирования h: " + str(h))
    print("2) Значение интеграла по формуле трапеций c шагом 2h: " + str(trapezes(a, b, 2 * h)))
    print("2) Значение интеграла по формуле трапеций c шагом h: " + str(trapezes(a, b, h)))
