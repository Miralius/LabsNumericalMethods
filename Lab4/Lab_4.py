from Lab4.functions import *

if __name__ == "__main__":
    print("Лабораторная работа №4, вариант №14")
    h = 0.1
    a = 0
    b = 2
    y0 = 2
    x0 = 0
    h = define_step(h, a, b, y0, x0)
    print("Оптимальный шаг h = " + str(h))
