import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import numpy as np

from scipy.integrate import ode


XLIM = (-4, 4)
YLIM = (-4, 4)


def f(x, ys):
    """ Правая часть ДУ y''=f(x, y, y'),
        представленная системой в нормальной форме.
        ys: (y, dy) - вектор неизвестных,
        x: независимая переменная.
    """

    y, dy = ys
    return [dy, np.sin(y)]
    # return [dy, 50 * (1 - y**2) * dy - y]


def dsolve(func, y0, dy0, x0):
    """ Численное решение ДУ с помощью класса ode """

    de = ode(func)
    de.set_integrator('dop853')
    # de.set_integrator('zvode', method='bdf')

    dt = 0.05
    soln = [[x0], [y0]]

    # интегрирование "вправо" от начальной точки
    de.set_initial_value([y0, dy0], x0)
    while de.successful() and de.t <= XLIM[1]:
        de.integrate(de.t + dt)
        soln[0].append(de.t)
        soln[1].append(de.y[0])

    # интегрирование "влево" от начальной точки
    de.set_initial_value([y0, dy0], x0)
    while de.successful() and de.t >= XLIM[0]:
        de.integrate(de.t - dt)
        soln[0].insert(0, de.t)
        soln[1].insert(0, de.y[0])

    return soln


def on_move(event, ax, line, tang, dot):
    """ Обработчик событий мыши """

    global data

    x1 = event.xdata
    y1 = event.ydata

    if x1 is None or y1 is None:  # мышь за пределами системы координат
        dot.set_data([], [])
        line.set_data([], [])
        tang.set_data([], [])
        ax.set_title("")
        ax.figure.canvas.draw_idle()
        return

    if len(data) == 0:  # первое начальное условие
        dot.set_data([x1], [y1])
        ax.set_title(f"y({x1:.2f})={y1:.2f}")

        if event.button == 1:
            data = (x1, y1)

    elif len(data) == 2:  # второе начальное условие
        x0, y0 = data

        delta_x = x0 - x1
        delta_y = y0 - y1
        dot.set_data([x0], [y0])
        tang.set_data([x1, x0 + delta_x], [y1, y0 + delta_y])

        if  delta_x == 0:  # деление на ноль невозможно
            line.set_data([], [])
            ax.figure.canvas.draw_idle()
            return

        dy0 = delta_y / delta_x
        soln = dsolve(f, y0, dy0, x0)

        line.set_data(soln[0], soln[1])

        ax.set_title(f"y({x0:.2f})={y0:.2f},  y'({x0:.2f})={dy0:.2f}")

        if event.button == 1:
            ax.plot(soln[0], soln[1], 'r')
            print(f"y({x0:.2f})={y0:.2f},  y'({x0:.2f})={dy0:.2f}")

        elif event.button == 3:
            dot.set_data([x1], [y1])
            tang.set_data([], [])
            line.set_data([], [])
            ax.set_title(f"y({x1:.2f})={y1:.2f}")
            data = ()

    ax.figure.canvas.draw_idle()


def init_axes():
    fig = plt.figure()
    ax = fig.add_subplot(xlim=XLIM, ylim=YLIM)

    ax.set_aspect('equal')
    ax.grid()

    ax.hlines(0, *XLIM, lw=0.5)
    ax.vlines(0, *YLIM, lw=0.5)

    dot, = ax.plot([], [], '.g')
    tang, = ax.plot([], [], 'g', lw=2)
    line, = ax.plot([], [], 'm', lw=2)

    callback = lambda event: on_move(event, ax, line, tang, dot)
    fig.canvas.mpl_connect('button_press_event', callback)
    fig.canvas.mpl_connect('motion_notify_event', callback)

    plt.show()


if __name__ == "__main__":
    data = ()
    init_axes()
