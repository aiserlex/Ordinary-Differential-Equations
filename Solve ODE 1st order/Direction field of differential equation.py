import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import numpy as np


XLIM = (-4, 4)
YLIM = (-4, 4)


def f(x, y):
    """ Правая часть ДУ y'=f(x, y) """
    return x + y**2


def main():
    fig = plt.figure()

    ax = plt.axes(xlim=XLIM, ylim=YLIM)
    ax.set_aspect('equal')

    # оси координат
    ax.hlines(0, *XLIM, lw=0.5)
    ax.vlines(0, *YLIM, lw=0.5)

    # здесь можно изменить плотность сетки, в узлах которой строятся векторы
    x = np.linspace(*XLIM, 21)
    y = np.linspace(*YLIM, 21)
    X, Y = np.meshgrid(x, y)

    # нормирующий множитель, приводит векторы поля к единой длине
    norm = np.hypot(1, f(X, Y))

    # поле направлений
    kwargs = {'angles':'xy', 'width':0.002, 'pivot':'mid'}
    ax.quiver(X, Y, 1/norm, f(X, Y)/norm, **kwargs)

    plt.show()


if __name__ == "__main__":
    main()
