import matplotlib
matplotlib.use('TkAgg')

from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

def f(x, y):
    # y <--> y[0];  y' <--> y[1]
    return [y[1], np.sin(y[0])]

def on_move(event):
    global x0, y0

    if not event.xdata or not event.ydata:  # выход курсора за пределы области
        dot.set_data([], [])
        line.set_data([], [])
        fig.canvas.draw_idle()
        return

    if event.button == 1:  # инициализация начальной точки
        x0 = event.xdata
        y0 = event.ydata
        dot.set_data([x0], [y0])
        line.set_data([], [])
    elif x0 is not None:  # если начальная точка инициализирована
        delta_x = event.xdata - x0
        delta_y = event.ydata - y0

        if  delta_x == 0:  # деление на ноль запрещено
            return
            
        dy0 = delta_y / delta_x
        
        de = ode(f)
        de.set_integrator('dop853')
        # de.set_integrator('zvode', method='bdf')
        
        dt = 0.05
        sol = []
        
        de.set_initial_value([y0, dy0], x0)
        while de.successful() and de.t <= xlim.end:
            de.integrate(de.t + dt)
            sol.append((de.t, de.y[0]))
        
        de.set_initial_value([y0, dy0], x0)
        while de.successful() and de.t >= xlim.start:
            de.integrate(de.t - dt)
            sol.append((de.t, de.y[0]))
        
        sol.sort(key=lambda x: x[0])
        sol = list(zip(*sol))
        
        if event.button == 3:  # зафиксировать интегральную кривую
            ax.plot(sol[0], sol[1], 'r')
        else:  # текущая интегральная кривая
            dot.set_data([x0], [y0])
            line.set_data(sol[0], sol[1])
            print(f"y'({x0:.2f})={dy0:.2f}")
    
    fig.canvas.draw_idle()


Lim = namedtuple('Lim', ['start', 'end'])
xlim = Lim(-10, 10)
ylim = Lim(-10, 10)

x0 = None
y0 = None

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

plt.show()
