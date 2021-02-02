import matplotlib
matplotlib.use('TkAgg')

from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

def f(x, y):
    """ Правая часть дифференциального уравнения 
        y'' = f(x, y, y'), записанная в форме нормальной
        системы дифференциальных уравнений
        Здесь y <--> y[0];  y' <--> y[1]
    """
    return [y[1], x+y[0]+y[1]]

def on_move(event):
    x0 = event.xdata
    y0 = event.ydata
    dy0 = 0

    if not x0 or not y0:  # выход курсора за пределы области
        line.set_data([], [])
        fig.canvas.draw_idle()
        return

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
    
    if event.button:  # зафиксировать интегральную кривую
        ax.plot(sol[0], sol[1], 'r')
    else:  # текущая интегральная кривая
        line.set_data(sol[0], sol[1])
    
    fig.canvas.draw_idle()


Lim = namedtuple('Lim', ['start', 'end'])
xlim = Lim(-4, 4)
ylim = Lim(-4, 4)

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

plt.show()
