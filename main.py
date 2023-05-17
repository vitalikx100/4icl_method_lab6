from numpy import *
from sympy import *
import sympy as sm
import numpy as np
import math

x = Symbol('x')

def task():
    f = x ** 2 * ln(x)
    a = 1
    b = 2
    M = max(abs(diff(f, x, x).subs(x, xi)) for xi in linspace(a, b, 10))
    print("M =", M)

    eps = 0.001
    h = ((12 * eps) / (M * abs(b - a))) ** 0.5
    print("Шаг интегрирования h =", h, "\n")

    n = int((b - a) / h)
    x_int_h = 0
    x_int_2h = 0

    for i in range(n + 1):
        if ((i == 0) or (i == (n + 1))):
            x_int_h += (f.evalf(subs={x: -0.5 + i * h}) / 2)
        else:
            x_int_h += f.evalf(subs={x: -0.5 + i * h})

    x_int_h *= h

    for i in range(n + 1):
        if ((i == 0) or (i == (n + 1))):
            x_int_2h += (f.evalf(subs={x: -0.5 + i * 2 * h}) / 2)
        else:
            x_int_2h += f.evalf(subs={x: -0.5 + i * 2 * h})

    x_int_2h *= (2 * h)

    print("Интеграл по формуле трапеций с шагом h: ", x_int_h, "\n")
    print("Интеграл по формуле трапеций с шагом 2h: ", x_int_2h, "\n")

    print("Погрешность по правилу Рунге: ", abs(x_int_2h - x_int_h) / 3, "\n")

    simpson_int_h = 0
    simpson_int_2h = 0

    for i in range(n + 1):
        if ((i == 0) or (i == (n + 1))):
            simpson_int_h += f.subs(x, -0.5 + i * h)
        elif (i % 2 == 0):
            simpson_int_h += (2 * f.subs(x, -0.5 + i * h))
        else:
            simpson_int_h += (4 * f.subs(x, -0.5 + i * h))

    simpson_int_h *= (h / 3)

    for i in range(n + 1):
        if ((i == 0) or (i == (n + 1))):
            simpson_int_2h += f.evalf(subs={x: -0.5 + i * 2 * h})
        elif (i % 2 == 0):
            simpson_int_2h += (2 * f.subs(x, -0.5 + i * 2 * h))
        else:
            simpson_int_2h += (4 * f.subs(x, -0.5 + i * 2 * h))

    simpson_int_2h *= (2 * h / 3)

    print("Интеграл методом Симпсона с шагом h: ", simpson_int_h, "\n")
    print("Интеграл методом Симпсона с шагом 2h: ", simpson_int_2h, "\n")

    print("Уточненная погрешность по правилу Рунге: ", abs(simpson_int_2h - simpson_int_h) / 15, "\n")

    print ("Библиотечный способ: ", sm.integrate(f, (x, 1, 2)))

if __name__ == '__main__':
    task()