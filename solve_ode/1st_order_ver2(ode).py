import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import ode


XLIM = (-4, 4)
YLIM = (-4, 4)


def f(x, y):
    """ Правая часть ДУ y'=f(x, y) """
    return x / 4 - 1 / (1 + y**2)


def dsolve(func, x0, y0):
    """ Численное решение ДУ с помощью класса ode """

    de = ode(f)
    de.set_integrator('dop853')

    dt = 0.05  # шаг интегрирования
    soln = [[x0], [y0]]  # переменная для решения

    # интегрирование "вправо" от начальной точки
    de.set_initial_value(y0, x0)
    while de.successful() and de.t <= XLIM[1]:
        de.integrate(de.t + dt)
        soln[0].append(de.t)
        soln[1].append(de.y[0])

    # интегрирование "влево" от начальной точки
    de.set_initial_value(y0, x0)
    while de.successful() and de.t >= XLIM[0]:
        de.integrate(de.t - dt)
        soln[0].insert(0, de.t)
        soln[1].insert(0, de.y[0])

    return soln


def on_move(event, fig, ax, line):
    """ Обработчик событий мыши """

    # начальные данные
    x0 = event.xdata
    y0 = event.ydata

    # выход курсора за пределы системы координат
    if not x0 or not y0:
        line.set_data([], [])
        fig.canvas.draw()
        return

    soln = dsolve(f, x0, y0)

    if event.button:
        ax.plot(soln[0], soln[1], 'r')
    else:
        line.set_data(soln[0], soln[1])

    fig.canvas.draw()


def dir_field(ax):
    """ Построение поля направлений """

    x_points = np.linspace(*XLIM, 17)
    y_points = np.linspace(*YLIM, 17)
    x, y = np.meshgrid(x_points, y_points)

    # нормирующий множитель, чтобы все векторы поля
    # имели одинаковую длину
    norm = np.hypot(1, f(x, y))

    # поле направлений
    kwargs = {'angles':'xy', 'width':0.002, 'pivot':'mid'}
    ax.quiver(x, y, 1/norm, f(x, y)/norm, **kwargs)


def init_coord_system():
    """ Инициализация системы координат """

    fig = plt.figure()

    ax = plt.axes(xlim=XLIM, ylim=YLIM)
    ax.set_aspect('equal')
    ax.grid()

    # оси координат
    ax.hlines(0, *XLIM, lw=0.5)
    ax.vlines(0, *YLIM, lw=0.5)

    return fig, ax


def main():
    fig, ax = init_coord_system()

    dir_field(ax)

    # линия, которая будет отрисовывать график решения
    # при движении мыши
    line, = ax.plot([], [], 'm', lw=2)

    # подключение обработчика событий
    callback = lambda e: on_move(e, fig, ax, line)
    fig.canvas.mpl_connect('motion_notify_event', callback)
    fig.canvas.mpl_connect('button_press_event', callback)

    plt.show()


if __name__ == "__main__":
    main()
