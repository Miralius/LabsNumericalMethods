import numpy
from sympy import symbols, exp, sqrt, diff, N, solve, integrate


def function(x, y):
    return y * (0.5 * x * y - 1)


def func():
    x, y = symbols('x y')
    return y * (0.5 * x * y - 1)


def delta(x, y, h):
    phi0 = h * function(x, y)
    phi1 = h * function(x + h / 2, y + phi0 / 2)
    phi2 = h * function(x + h / 2, y + phi1 / 2)
    phi3 = h * function(x + h, y + phi2)
    return (phi0 + 2 * phi1 + 2 * phi2 + phi3) / 6


def define_step(h, a, b, y0, x0):
    print("Определение шага…")
    h0 = h
    y1 = y0 + delta(x0, y0, h0)
    x1 = x0 + h0
    y2 = y1 + delta(x1, y1, h0)
    y2_ = y0 + delta(x0, y0, 2 * h0)
