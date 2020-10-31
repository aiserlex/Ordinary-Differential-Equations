import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import numpy as np

from scipy.integrate import ode


XLIM = (-4, 4)
YLIM = (-4, 4)


def f(x, ys):
    """ Правая часть ДУ y'''=f(x, y, y', y''),
        представленная системой в нормальной форме.
        ys: (y, dy, d2y) - вектор неизвестных,
        x: независимая переменная.
    """

    y, dy, d2y = ys
    return [dy, d2y, 3 * dy * d2y**2 / (1 + dy**2)]


def dsolve(func, y0, dy0, d2y0, x0):
    """ Численное решение ДУ с помощью класса ode """

    de = ode(func)
    de.set_integrator('dop853')
    # de.set_integrator('zvode', method='bdf')

    dt = 0.05
    soln = [[x0], [y0]]

    # интегрирование "вправо" от начальной точки
    de.set_initial_value([y0, dy0, d2y0], x0)
    while de.successful() and de.t <= XLIM[1]:
        de.integrate(de.t + dt)
        soln[0].append(de.t)
        soln[1].append(de.y[0])

    # интегрирование "влево" от начальной точки
    de.set_initial_value([y0, dy0, d2y0], x0)
    while de.successful() and de.t >= XLIM[0]:
        de.integrate(de.t - dt)
        soln[0].insert(0, de.t)
        soln[1].insert(0, de.y[0])

    return soln


def on_move(event, ax, line, circ, tang, dot):
    """ Обработчик событий мыши """

    global data

    x2 = event.xdata
    y2 = event.ydata

    if x2 is None or y2 is None:  # мышь за пределами системы координат
        dot.set_data([], [])
        tang.set_data([], [])
        line.set_data([], [])
        circ.set_radius(0)
        ax.set_title("")
        ax.figure.canvas.draw_idle()
        return

    if len(data) == 0:  # первое начальное условие
        dot.set_data([x2], [y2])
        ax.set_title(f"y({x2:.2f})={y2:.2f}")

        if event.button == 1:
            data = (x2, y2)

    elif len(data) == 2:  # второе начальное условие
        x0, y0 = data

        dot.set_data([x0], [y0])
        tang.set_data([2*x0 - x2, x2], [2*y0 - y2, y2])

        if  x2 - x0 == 0:  # деление на ноль невозможно
            line.set_data([], [])
            ax.figure.canvas.draw_idle()
            return

        dy0 = (y2 - y0) / (x2 - x0)

        ax.set_title(f"y({x0:.2f})={y0:.2f},  y'({x0:.2f})={dy0:.2f}")

        if event.button == 1:
            data = data + (x2, y2, dy0)

    elif len(data) == 5:  # третье начальное условие
        x0, y0, x1, y1, dy0 = data

        dot.set_data([x0], [y0])
        tang.set_data([2*x0 - x1, x1], [2*y0 - y1, y1])

        if x2 - x0 == 0:  # деление на ноль запрещено
            line.set_data([], [])
            ax.figure.canvas.draw_idle()
            return

        if dy0 == 0:
            x_c = x0
            y_c = 0.5*(y2+y0+(x2-x0)**2/(y2-y0))
            R = abs(y0-y_c)
        else:
            # Угловой коэффициент прямой, на кот. должен лежать центр окр.
            k = -1/dy0

            # Координаты центра окружности
            x_c = .5 * (y2**2 - y0**2 + x2**2 - x0**2 + 2 * (y0 - k*x0) * (y0 - y2)) / (k*(y2 - y0) + x2 - x0)
            y_c = k * (x_c - x0) + y0

            # Радиус окружности
            R = np.hypot(y0 - y_c, x0 - x_c)

        # Отрисовка окружности
        circ.center = (x_c, y_c)
        circ.set_radius(R)
        print(x_c, y_c, R)

        # Начальное значение 2ой производной
        d2y0 = (1+dy0**2)**(3/2)/R*np.sign(y_c-y0)

        ax.set_title(f"y({x0:.2f})={y0:.2f},  y'({x0:.2f})={dy0:.2f},  y''({x0:.2f})={d2y0:.2f}")

        soln = dsolve(f, y0, dy0, d2y0, x0)
        line.set_data(soln[0], soln[1])

        if event.button == 1:
            ax.plot(soln[0], soln[1], 'r')
            print(f"y({x0:.2f})={y0:.2f},  y'({x0:.2f})={dy0:.2f}, y''({x0:.2f})={d2y0:.2f}")

        elif event.button == 3:
            dot.set_data([x2], [y2])
            tang.set_data([], [])
            line.set_data([], [])
            circ.set_radius(0)
            ax.set_title(f"y({x2:.2f})={y2:.2f}")
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
    tang, = ax.plot([], [], 'g', lw=1)
    circ = Circle((0, 0), 0, color='g', lw=2, fill=False)
    ax.add_patch(circ)
    line, = ax.plot([], [], 'm', lw=2)

    callback = lambda event: on_move(event, ax, line, circ, tang, dot)
    fig.canvas.mpl_connect('button_press_event', callback)
    fig.canvas.mpl_connect('motion_notify_event', callback)

    plt.show()


if __name__ == "__main__":
    data = ()
    init_axes()
