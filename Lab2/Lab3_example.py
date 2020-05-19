import math as m
import numpy as np
import matplotlib.pyplot as pl


def f(x):
    return (np.pi * x)/(10 + x/2)

def M4(x):
    return 48 * np.pi * (x / (x + 20) - 1)/(x + 20)**4


def f2(x0, x1):
    return (f(x1) - f(x0)) / (x1 - x0)

def fw(x,x0, x1,x2,x3):
    return (x-x0)*(x-x1)*(x-x2)*(x-x3)

def f3(x0, x1, x2):
    return (f2(x1, x2) - f2(x0, x1)) / (x2 - x0)


def f4(x0, x1, x2, x3):
    return (f3(x1, x2, x3) - f3(x0, x1, x2)) / (x3 - x0)


def lagrange(x, otr):
    i = 0
    L = 0
    while i < 4:
        j = 0
        mult = 1
    while j < 4:
        if i != j:
            mult *= (x - otr[j]) / (otr[i] - otr[j])
        j += 1
        L += f(otr[i]) * mult
        i += 1
    return L


def Newton(x, otr):
        return f(otr[0]) + f2(otr[0], otr[1]) * (x - otr[0]) + f3(otr[0], otr[1], otr[2]) * (x - otr[0]) * (
        x - otr[1]) + f4(otr[0], otr[1], otr[2], otr[3]) * (x - otr[0]) * (x - otr[1]) * (x - otr[2])


a = 0
b = 10
otr1 = np.array([ 1.2, 1.9, 2.9, 3.5])

w = np.linspace(0, 10, 100)
h = 0.1
pl.figure(1)
#pl.plot(w, M4(w)*np.abs(fw(w,otr1[0], otr1[1], otr1[2], otr1[3])), label=u'График производной 4-го порядка f(x)*w')
pl.plot(w, M4(w), label=u'график производной 4-ого порядка')

pl.grid()
pl.legend()





otr = np.array([b-3*h, b-2*h, b-h, b])

print('Отрезок: [' + str(b-3*h) + ',' + str(b) + ']')
print('Шаг: ' + str(h))

k = b-3*h
x = np.linspace(b-3*h, b, 100)
pl.figure(2)
lol = lagrange(x, otr)
pl.title("Графики исходной функции и многочлена Лагранжа")
pl.plot(x, f(x), linewidth=3, color='red')
pl.plot(x, lagrange(x, otr), linewidth=1, color='blue')
while k <= b:
    pl.scatter(k, f(k))
    k += h
    pl.grid()


a11 = f(otr[0])
a12 = f(otr[1])
a13 = f(otr[2])
b11 = (f(otr[1]) - a11) / h
b12 = (f(otr[2]) - a12) / h
b13 = (f(otr[3]) - a13) / h


def L_spline(v):
    res = []
    for x in v:
        if x <= otr[1]:
            res.append(a11 + (x - otr[0]) * b11)
        elif x <= otr[2]:
            res.append(a12 + (x - otr[1]) * b12)
        elif x <= otr[3]:
            res.append(a13 + (x - otr[2]) * b13)
    return np.array(res)


A = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[1, h, h ** 2, h ** 3, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 1, h, h ** 2, h ** 3, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 1, h, h ** 2, h ** 3],
[0, 1, 2 * h, 3 * h ** 2, 0, -1, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 1, 2 * h, 3 * h ** 2, 0, -1, 0, 0],
[0, 0, 2, 6 * h, 0, 0, -2, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 2, 6 * h, 0, 0, -2, 0],
[0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 6 * h]])

b = np.array([[f(otr[0])],
[f(otr[1])],
[f(otr[1])],
[f(otr[2])],
[f(otr[2])],
[f(otr[3])],
[0],
[0], [0],
[0],
[0],
[0]])

res = np.linalg.solve(A, b)
res = res.ravel()

a31, b31, c31, d31, a32, b32, c32, d32, a33, b33, c33, d33 = res


def C_spline(vect):
    res = []
    for v in vect:
        if v <= otr[1]:
            res.append(a31 + b31 * (v - otr[0]) + c31 * (v - otr[0]) ** 2 + d31 * (v - otr[0]) ** 3)
        elif v <= otr[2]:
            res.append(a32 + b32 * (v - otr[1]) + c32 * (v - otr[1]) ** 2 + d32 * (v - otr[1]) ** 3)
        elif v <= otr[3]:
            res.append(a33 + b33 * (v - otr[2]) + c33 * (v - otr[2]) ** 2 + d33 * (v - otr[2]) ** 3)
    return np.array(res)


A = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0],
[1, h, h ** 2, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 1, 0, 0, 0, 0, 0],
[0, 0, 0, 1, h, h ** 2, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 1, 0, 0],
[0, 0, 0, 0, 0, 0, 1, h, h ** 2],
[0, 1, 2 * h, 0, -1, 0, 0, 0, 0],
[0, 0, 0, 0, 1, 2 * h, 0, -1, 0],
[0, 0, 1, 0, 0, 0, 0, 0, 0]])

b = np.array([[f(otr[0])],
[f(otr[1])],
[f(otr[1])],
[f(otr[2])],
[f(otr[2])],
[f(otr[3])],
[0],
[0],
[0]])

res = np.linalg.solve(A, b)
res = res.ravel()

a21, b21, c21, a22, b22, c22, a23, b23, c23 = res


def P_spline(vect):
    res = []
    for v in vect:
        if v <= otr[1]:
            res.append(a21 + b21 * (v - otr[0]) + c21 * (v - otr[0]) ** 2)
        elif v <= otr[2]:
            res.append(a22 + b22 * (v - otr[1]) + c22 * (v - otr[1]) ** 2)
        elif v <= otr[3]:
            res.append(a23 + b23 * (v - otr[2]) + c23 * (v - otr[2]) ** 2)
    return np.array(res)


pl.figure(3)
pl.title("Абсолютные погрешности сплайнов")
pl.xlabel("x")
pl.ylabel("R_n(x)")
pl.plot(x, np.absolute(f(x) - L_spline(x)), label="S1(x)", linewidth=3, color='red')
pl.plot(x, np.absolute(f(x) - P_spline(x)), label="S2(x)", color='blue')
pl.plot(x, np.absolute(f(x) - C_spline(x)), label="S3(x)", color='green')
pl.legend()
pl.grid()

pl.figure(4);


x = np.linspace(1.7, 3.5, 100)
otr = np.array([ 1.7, 2.3, 2.9, 3.5])
pl.title("Абслоютные погрешности многочленов Лагранжа и Ньютона ")
pl.xlabel("x")
pl.ylabel("R_n(x)")
pl.plot(x, abs(f(x) - Newton(x, otr)), label="Newton",linewidth=3, color='pink')
pl.plot(x, abs(f(x) - lagrange(x, otr)), label="Lagrange", color='black')
pl.legend()
pl.grid()