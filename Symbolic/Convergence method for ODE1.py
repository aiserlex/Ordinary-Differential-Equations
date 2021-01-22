# -*- coding: utf-8 -*-

import sympy as sym
import matplotlib.pyplot as plt

x = sym.Symbol('x')
y = -1

colors = ('r', 'g', 'b', 'c', 'm', 'k')
lines = []

for i in range(1, 7):
    y = -1 + sym.integrate(y + x, (x, 0, x))
    print('=====', i, '=====')
    sym.pprint(y)
    lines.append((y, (x, -3, 3)))

curves = sym.plot(*lines)

for curve, color in zip(curves, colors):
    curve.line_color = color

curves.show()
plt.grid()
