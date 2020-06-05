from Lab2.functions import *

if __name__ == "__main__":
    print("Лабораторная работа №3, вариант #14")
    a = 0
    b = 10
    eps = 10 ** (-4)
    print("a = " + str(a))
    print("b = " + str(b))
    print("Точность = " + str(eps))

    print("Первое задание ---------")
    h = define_step(a, b, eps)
    print("Шаг = " + str(h))
    x, y = discretize(a, b, h)
    plot(x, numpy.array([y]), numpy.array(["f(x) = pi * x / (10 + x / 2)"]))
    print("x = " + str(x))
    print("Всего точек: " + str(len(x)))
    number_x = int(input("Выберите начало отрезка x для интерполяции полиномами Лагранжа: "))
    x_L = numpy.array([x[number_x + i] for i in range(4)])
    y_L = numpy.array([interpolate_with_lagrange_polynomial(x_L[i], x_L) for i in range(len(x_L))])
    x, y = discretize(x_L[0], x_L[3], (x_L[3] - x_L[0]) / 100)
    y_int_L = numpy.array([interpolate_with_lagrange_polynomial(x[i], x_L) for i in range(len(x))])
    plot(x, numpy.array([y, y_int_L]), numpy.array(["f(x)", "interpolated function by Lagrange polynomial"]))
    plot(x, numpy.array([numpy.array([abs(y_int_L[i] - y[i]) for i in range(len(x))])]),
         numpy.array(["absolute error y(x) & y_L(x)"]))

    y_N = numpy.array([interpolate_with_newton_method(x_L[i], x_L, y_L, h) for i in range(len(x_L))])
    y_int_N = numpy.array([interpolate_with_newton_method(x[i], x_L, y_L, h) for i in range(len(x))])
    plot(x, numpy.array([y, y_int_L, y_int_N]),
         numpy.array(["f(x)", "Lagrange polynomial interpolation", "Newton polynomial interpolation"]))
    plot(x, numpy.array([numpy.array([abs(y_int_N[i] - y[i]) for i in range(len(x))])]),
         numpy.array(["absolute error y(x) & y_N(x)"]))
    plot(x, numpy.array([numpy.array([abs(y_int_N[i] - y[i]) for i in range(len(x))]),
                         numpy.array([abs(y_int_L[i] - y[i]) for i in range(len(x))])]),
         numpy.array(["absolute error y(x) & y_N(x)", "absolute error y(x) & y_L(x)"]))
    plot(x, numpy.array([numpy.array([abs(y_int_N[i] - y_int_L[i]) for i in range(len(x))])]),
         numpy.array(["absolute error y_L(x) & y_N(x)"]))

    y_interpolated_linear_spline = numpy.array([interpolate_linear_spline(x[i], x_L, y_N) for i in range(len(x))])
    y_interpolated_parabolic_spline = numpy.array([interpolate_parabolic_spline(x[i], x_L, y_N) for i in range(len(x))])
    y_interpolated_cubic_spline = numpy.array([interpolate_cubic_spline(x[i], x_L, y_N) for i in range(len(x))])
    plot(x, numpy.array([y_interpolated_linear_spline, y_interpolated_parabolic_spline, y_interpolated_cubic_spline]),
         numpy.array(["linear spline", "parabolic spline", "cubic spline"]))
    plot(x, numpy.array([numpy.array([abs(y_interpolated_linear_spline[i] - y[i]) for i in range(len(x))]),
                         numpy.array([abs(y_interpolated_parabolic_spline[i] - y[i]) for i in range(len(x))]),
                         numpy.array([abs(y_interpolated_cubic_spline[i] - y[i]) for i in range(len(x))])]),
         numpy.array(["linear spline error", "parabolic spline error", "cubic spline error"]))
    plot(x, numpy.array([numpy.array([abs(y_interpolated_linear_spline[i] - y[i]) for i in range(len(x))]),
                         numpy.array([abs(y_interpolated_parabolic_spline[i] - y[i]) for i in range(len(x))]),
                         numpy.array([abs(y_interpolated_cubic_spline[i] - y[i]) for i in range(len(x))]),
                         numpy.array([abs(y_int_L[i] - y[i]) for i in range(len(x))]),
                         numpy.array([abs(y_int_N[i] - y[i]) for i in range(len(x))])]),
         numpy.array(["linear spline error", "parabolic spline error", "cubic spline error", "Lagrange error",
                      "Newton error"]))
    plot(x, numpy.array([y_interpolated_linear_spline, y_interpolated_parabolic_spline, y_interpolated_cubic_spline,
                         y_int_L, y_int_N, y]),
         numpy.array(["linear spline", "parabolic spline", "cubic spline", "Lagrange polynomial", "Newton polynomial",
                      "original function"]))
    plot(x, numpy.array([y_interpolated_linear_spline, y_interpolated_parabolic_spline, y_interpolated_cubic_spline,
                         y_int_L, y_int_N, y]),
         numpy.array(["linear spline", "parabolic spline", "cubic spline", "Lagrange polynomial", "Newton polynomial",
                      "original function"]), numpy.array([5.6995, 5.7, 1.39333, 1.39356]))
    plot(x, numpy.array([y_interpolated_linear_spline, y_interpolated_parabolic_spline, y_interpolated_cubic_spline,
                         y_int_L, y_int_N, y]),
         numpy.array(["linear spline", "parabolic spline", "cubic spline", "Lagrange polynomial", "Newton polynomial",
                      "original function"]), numpy.array([5.69998, 5.7, 1.393541, 1.393547]))
