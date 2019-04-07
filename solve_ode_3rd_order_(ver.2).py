"""
    Программа для построения интегральных кривых дифференциального уравнения 3-го порядка,
    разрешенного относительно производной y''' = f(x, y, y', y'').
    
    Левая кнопка мыши - зафиксировать начальное условие или зафиксировать интегральную кривую.
    Правая кнопка мыши - сменить началные условия.
    
"""

import matplotlib
matplotlib.use('TkAgg')

from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.integrate import ode

def f(x, y):
    """ Правая часть дифференциального уравнения y'''=(x, y, y', y'')
        Здесь y <--> y[0];  y' <--> y[1]; y'' <--> y[2]
    """
    return [y[1], y[2], x + y[0] + y[1]+ y[2]]

def on_move(event):
    global x0, y0, x1, y1, dy0

    if not event.xdata or not event.ydata:  # выход курсора за пределы области
        line.set_data([], [])
        dot.set_data([], [])
        tang.set_data([], [])
        circ.set_radius(0)
        fig.canvas.draw_idle()
        return

    if x0 is None:  # инициализация 1-го начального условия
        dot.set_data([event.xdata], [event.ydata])
        ax.set_title(f"y({event.xdata:.2f})={event.ydata:.2f}")
        if event.button == 1:  
            x0 = event.xdata
            y0 = event.ydata
    elif x1 is None:  # инициализация 2-го начального условия
        
        # восстановление доп. построений, если они были удалены при выходе
        # мыши за пределы области        
        dot.set_data([x0], [y0])
        
        tang.set_data([2*x0 - event.xdata, event.xdata], [2*y0 - event.ydata, event.ydata])
        
        delta_x = event.xdata - x0
        delta_y = event.ydata - y0

        if  delta_x == 0:  # деление на ноль запрещено
            return

        ax.set_title(f"y({x0:.2f})={y0:.2f},  y'({x0:.2f})={delta_y/delta_x:.2f}")

        if event.button == 1:
            x1 = event.xdata
            y1 = event.ydata

            dy0 = delta_y / delta_x
        
    else:  # инициализация 3-го начального условия и построение интегральной кривой
        x2 = event.xdata
        y2 = event.ydata
        
        # восстановление доп. построений, если они были удалены при выходе
        # мыши за пределы области
        dot.set_data([x0], [y0])
        tang.set_data([2*x0 - x1, x1], [2*y0 - y1, y1])
        
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

        # Начальное значение 2ой производной
        d2y0 = (1+dy0**2)**(3/2)/R*np.sign(y_c-y0)
        
        ax.set_title(f"y({x0:.2f})={y0:.2f},  y'({x0:.2f})={dy0:.2f},  y''({x0:.2f})={d2y0:.2f}")
        
        de = ode(f)
        de.set_integrator('dop853')
        # de.set_integrator('zvode', method='bdf')
        
        dt = 0.05
        sol = []
        
        de.set_initial_value([y0, dy0, d2y0], x0)
        while de.successful() and de.t <= xlim.end:
            de.integrate(de.t + dt)
            sol.append((de.t, de.y[0]))
        
        de.set_initial_value([y0, dy0, d2y0], x0)
        while de.successful() and de.t >= xlim.start:
            de.integrate(de.t - dt)
            sol.append((de.t, de.y[0]))
        
        sol.sort(key=lambda x: x[0])
        sol = list(zip(*sol))
        
        if event.button == 1:  # зафиксировать интегральную кривую
            ax.plot(sol[0], sol[1], 'r')
        elif event.button == 3:  # сменить начальную точку
            x0 = event.xdata
            y0 = event.ydata
            x1 = None
            y1 = None
            dy0 = None
            dot.set_data([x0], [y0])
            tang.set_data([], [])
            circ.set_radius(0)
            line.set_data([], [])
            ax.set_title(f"y({x0:.2f})={y0:.2f}")
        else:  # текущая интегральная кривая
            line.set_data(sol[0], sol[1])
            print(f"y''({x0:.2f})={d2y0:.2f}")
    
    fig.canvas.draw_idle()


Lim = namedtuple('Lim', ['start', 'end'])
xlim = Lim(-5, 5)
ylim = Lim(-5, 5)

x0 = None
y0 = None
x1 = None
y1 = None
dy0 = None

fig, ax = plt.subplots()

ax.grid()
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_aspect('equal')

ax.hlines(0, *xlim, lw=0.5)
ax.vlines(0, *ylim, lw=0.5)

fig.canvas.mpl_connect('button_press_event', on_move)
fig.canvas.mpl_connect('motion_notify_event', on_move)

line, = ax.plot([], [], 'm')
dot, = ax.plot([], [], '.m')
tang, = ax.plot([], [], 'g', lw=0.5)
circ = Circle((0, 0), 0, color='g', lw=0.5, fill=False)
ax.add_patch(circ)

plt.show()



