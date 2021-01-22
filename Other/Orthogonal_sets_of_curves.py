import matplotlib.pyplot as plt
import numpy as np

from numpy import exp, sin, cos, log
from scipy.optimize import root  # для численного поиска корней уравнений


XLIM = (-4, 4)  # область системы
YLIM = (-4, 4)  # координат


def f1(k, x, y):
    """ Исходное семейство кривых, заданное в виде f1(k, x, y)=0 """
    return y - k * (x + 1) * exp(-x)


def f2(k, x, y):
    """ Ортогональные кривые, заданные в виде f2(k, x, y)=0 """
    return x**2 - k * exp(y**2 - 2*x)


def onmove(event, fig, ax):
    """ Обработчик события: движение мыши """
    x0 = event.xdata  # текущие координаты
    y0 = event.ydata  # мыши на плоскости

    if x0 is None or y0 is None:  # если мышь за пределами плоскости,
        return                    # то ничего не делать

    x = np.linspace(*XLIM, 300)  # чем больше узлов, тем точнее график, а
    y = np.linspace(*YLIM, 300)  # чем меньше узлов, тем выше отклик программы

    X, Y = np.meshgrid(x, y)  # сетка узловых точек на плоскости

    sol1 = root(f1, 2, args=(x0, y0))  # поиск произовльной постоянной k
    Z1 = f1(sol1.x[0], X, Y)  # значения функции f1 при найденном k

    sol2 = root(f2, 2, args=(x0, y0))  # поиск произовльной постоянной k
    Z2 = f2(sol2.x[0], X, Y)  # значения функции f2 при найденном k

    # plt.cla()  # очистить прошлые графики; если закомментировать эти две строчки,
    # ax.grid()  # то будет видно всё семейство кривых

    ax.contour(X, Y, Z1, [0], colors='red', alpha=0.4)  # построение графика f1 красным
    ax.contour(X, Y, Z2, [0], colors='blue', alpha=0.4)  # построение графика f2 синим

    fig.canvas.draw()  # отрисовка графиков


def main():
    fig = plt.figure()
    ax = plt.axes(xlim=XLIM, ylim=YLIM)  # инициализация системы координат
    ax.set_aspect("equal")  # соотношение масштабов по осям: 1 к 1
    ax.grid()  # отображать сетку системы координат

    # подключение обработчкика события мыши: движение мыши
    callback = lambda event: onmove(event, fig, ax)
    fig.canvas.mpl_connect("motion_notify_event", callback)

    plt.show()


if __name__ == "__main__":
    main()
