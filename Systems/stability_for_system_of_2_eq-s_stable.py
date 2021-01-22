""" Пример устойчивой системы.
    Два решения, с практически одинаковыми начальными условиями
    с течением времени остаются похожими (при любом t>0).

"""

import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import ode


XLIM = (-4, 4)
YLIM = (-4, 4)
ZLIM = (-4, 4)


# правая часть первого уравнения
def g1(t, y1, y2):
    x = y1
    y = y2
    return y + np.sin(x)


# правая часть второго уравнения
def g2(t, y1, y2):
    x = y1
    y = y2
    a = -10
    b = -3
    return a*x+b*y


def f(t, y):
    return [g1(t, *y), g2(t, *y)]


def dsolve(func, y1_0, y2_0, t_0):
    """ Численное решение системы """

    de = ode(func)
    de.set_integrator('dop853')
    # de.set_integrator('zvode', method='bdf')

    dt = 0.01
    soln = [[t_0], [y1_0], [y2_0]]

    de.set_initial_value([y1_0, y2_0], t_0)
    while de.successful() and de.t <= XLIM[1]:
        de.integrate(de.t + dt)
        soln[0].append(de.t)
        soln[1].append(de.y[0])
        soln[2].append(de.y[1])

    return soln


def init_axes():
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax1 = fig.add_subplot(121, xlim=XLIM, ylim=YLIM, zlim=ZLIM, projection='3d')
    ax1.set_xlabel('t')
    ax1.set_ylabel('y1')
    ax1.set_zlabel('y2')
    ax1.set_title('Интегральные кривые')

    ax2 = fig.add_subplot(122, xlim=YLIM, ylim=ZLIM)
    ax2.set_aspect('equal')
    ax2.set_xlabel('y1')
    ax2.set_ylabel('y2')
    ax2.set_title('Траектории')

    return ax1, ax2


def dir_field(ax2, t_0):
    """ Поле направлений фазового пространства """
    y1, y2 = np.meshgrid(np.linspace(*YLIM, 25), np.linspace(*ZLIM, 25))
    hyp = np.hypot(g1(t_0, y1, y2), g2(t_0, y1, y2))
    ax2.quiver(y1, y2, g1(t_0, y1, y2)/hyp, g2(t_0, y1, y2)/hyp, angles='xy')


def plots(ax1, ax2):
    t_0 = 0
    dir_field(ax2, t_0)

    y1_0 = -0.1; y2_0 = 0
    soln1 = dsolve(f, y1_0, y2_0, t_0)

    ax1.plot(*soln1, 'b')
    ax2.plot(soln1[1], soln1[2], 'b')
    ax2.plot([y1_0], [y2_0], 'b.')

    y1_0 = 0.1; y2_0 = 0.1
    soln2 = dsolve(f, y1_0, y2_0, t_0)

    ax1.plot(*soln2, 'r')
    ax2.plot(soln2[1], soln2[2], 'r')
    ax2.plot([y1_0], [y2_0], 'r.')


def main():
    ax1, ax2 = init_axes()
    plots(ax1, ax2)
    plt.show()


if __name__ == "__main__":
    main()
