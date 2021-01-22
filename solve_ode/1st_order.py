import matplotlib
matplotlib.use('TkAgg')

from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import ode


def f(x, y):
    """ Правая часть ДУ y'=f(x, y) """
    return x/4-1/(1+y**2)


def on_move(event):
    """ Обработчик событий мыши """
    
    # начальные данные
    x0 = event.xdata
    y0 = event.ydata

    # выход курсора за пределы системы координат
    if not x0 or not y0:
        line.set_data([], [])
        fig.canvas.draw()
        return
    
    dt = 0.05  # шаг интегрирования
    sol = []  # решение

    de = ode(f)
    de.set_integrator('dop853')

    # интегрирование "вправо" от начальной точки
    de.set_initial_value(y0, x0)
    while de.successful() and de.t <= xlim.end:
        de.integrate(de.t + dt)
        sol.append((de.t, de.y[0]))

    # интегрирование "влево" от начальной точки
    de.set_initial_value(y0, x0)    
    while de.successful() and de.t >= xlim.start:
        de.integrate(de.t - dt)
        sol.append((de.t, de.y[0]))
    
    sol.sort(key=lambda x: x[0])
    sol = list(zip(*sol))
    
    if event.button:
        ax.plot(sol[0], sol[1], 'r')
    else:
        line.set_data(sol[0], sol[1])
    fig.canvas.draw()


# прямоугольная область на плоскости
Lims = namedtuple('Lims', ['start', 'end'])
xlim = Lims(-5, 5)
ylim = Lims(-5, 5)

fig = plt.figure()

# подключение обработчика событий
fig.canvas.mpl_connect('motion_notify_event', on_move)
fig.canvas.mpl_connect('button_press_event', on_move)

ax = plt.axes(xlim=xlim, ylim=ylim)
ax.set_aspect('equal')

# оси координат
ax.hlines(0, xlim.start, xlim.end, lw=0.5)
ax.vlines(0, ylim.start, ylim.end, lw=0.5)

x = np.linspace(xlim.start, xlim.end, 21)
y = np.linspace(ylim.start, ylim.end, 21)
X, Y = np.meshgrid(x, y)

# нормирующий множитель, чтобы все векторы поля
# имели одинаковую длину
norm = np.hypot(1, f(X, Y))

# поле направлений
kwargs = {'angles':'xy', 'width':0.002, 'pivot':'mid'}
ax.quiver(X, Y, 1/norm, f(X, Y)/norm, **kwargs)

# линия, которая будет отрисовывать график решения
# при движении мыши
line, = ax.plot([], [], 'm')

plt.show()
