""" This program solves differential equations of 3rd order
    given in form y''' = f(x, y, y', y'').

    Define right side of equation f and region XLIM, YLIM by yourself.

    Default is y''' = 3 * y' * y''**2 / (1 + y'**2)
    in region [-4, 4] x [-4, 4] (equation of a circle).

    Start the program and define initial conditions
    using mouse clicks and mouse moving:
        1st click for freezing initial condition y(x_0)=y_0,
        2nd click for freezing initial condition y'(x_0)=y_1,
        3rd click for freezing initial condition y''(x_0)=y_2.

    Use right mouse button click to chose new initial conditions.
    Use middle mouse button click to clear all.
"""

import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import numpy as np

from scipy.integrate import ode


XLIM = (-4, 4)  # region for drawing
YLIM = (-4, 4)  # integral curves


def f(x, ys):
    """ Right part of equation y''' = f(x, y, y', y''),

        ys: (y, y', y''),
        x: independent variable.
    """
    y, dy, d2y = ys
    func = 3 * dy * d2y**2 / (1 + dy**2)  # change function in this line
                                          # x for independent variable
                                          # y for unknown function
                                          # dy for y'
                                          # and d2y for y''
    return [dy, d2y, func]


class DataManager:
    """ This class tracks mouse coordinates and computes dy and d2y
        which are 1st and 2nd derivatives respectively.
    """
    def __init__(self):
        self.data = [None for _ in range(3)]
        self.click_num = 0

    def update(self, item):
        item = self.modify(item)
        index = self.click_num if self.click_num <=2 else 2
        self.data[index] = item

    def inc(self):
        self.click_num += 1

    def modify(self, item):
        if self.click_num == 0:
            return item
        elif self.click_num == 1:
            x0, y0 = self.data[0]
            x1, y1 = item
            if x1 - x0 != 0:
                dy = (y1 - y0) / (x1 - x0)
            else:
                dy = np.nan
            return (*item, dy)
        elif self.click_num > 1:
            x0, y0 = self.data[0]
            x1, y1, dy = self.data[1]
            x2, y2 = item
            if x2 - x0 != 0:
                d2y = 2 * (y2 - y0 - dy * (x2 - x0)) / (x2 - x0)**2
            else:
                d2y = np.nan
            return (*item, d2y)

    def reset(self):
        self.data = [None for _ in range(3)]
        self.click_num = 0

    def get_click_num(self):
        return self.click_num

    def get_data(self, index):
        return self.data[index]


class InitAxes:
    """ Setup axes and event handlers """
    def __init__(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(xlim=XLIM, ylim=YLIM)
        self.tune_axes()

    def tune_axes(self):
        self.ax.set_aspect('equal')
        self.ax.grid()

    def set_events(self, callback):
        self.fig.canvas.mpl_connect('button_press_event', callback)
        self.fig.canvas.mpl_connect('motion_notify_event', callback)

    def get_ax(self):
        return self.ax


class Plots:
    """ Make all plots """
    def __init__(self, ax):
        self.ax = ax
        self.ax.hlines(0, *XLIM, lw=0.5)
        self.ax.vlines(0, *YLIM, lw=0.5)
        self.dot, = ax.plot([], [], '.g')
        self.tang, = ax.plot([], [], 'g', lw=1)
        self.parab, = ax.plot([], [], 'g', lw=1)
        self.curve, = ax.plot([], [], 'm', lw=2)
        self.added = []

    def reset(self):
        self.dot.set_data([], [])
        self.tang.set_data([], [])
        self.parab.set_data([], [])
        self.curve.set_data([], [])
        self.ax.set_title("")

    def draw_idle(self):
        self.ax.figure.canvas.draw_idle()

    def draw_dot(self, x0, y0):
        self.dot.set_data([x0], [y0])

    def draw_tang(self, x0, y0, x1, y1):
        self.tang.set_data([2*x0 - x1, x1], [2*y0 - y1, y1])

    def draw_parab(self, x0, y0, dy, d2y):
        x = np.linspace(*XLIM, 300)
        y = y0 + dy * (x - x0) + d2y * (x - x0)**2 / 2
        self.parab.set_data(x, y)

    def draw_curve(self, soln):
        self.curve.set_data(soln[0], soln[1])

    def add_new(self, soln):
        new, = self.ax.plot(soln[0], soln[1], 'r', lw=2)
        self.added.append(new)

    def set_title(self, title):
        self.ax.set_title(title)

    def clear(self):
        for item in self.added:
            item.set_data([], [])
        self.added.clear()



def dsolve(func, y0, dy0, d2y0, x0):
    """ Numerical solution with "ode" class """

    de = ode(func)
    de.set_integrator('dop853')
    # de.set_integrator('zvode', method='bdf')

    dt = 0.05
    soln = [[x0], [y0]]

    # integration to the right from start point
    de.set_initial_value([y0, dy0, d2y0], x0)
    while de.successful() and de.t <= XLIM[1]:
        de.integrate(de.t + dt)
        soln[0].append(de.t)
        soln[1].append(de.y[0])

    # integration to the left from start point
    de.set_initial_value([y0, dy0, d2y0], x0)
    while de.successful() and de.t >= XLIM[0]:
        de.integrate(de.t - dt)
        soln[0].insert(0, de.t)
        soln[1].insert(0, de.y[0])

    return soln


def on_move(event, data, plots, ax):
    """ Event handler (mouse move, mouse click) """
    xn = event.xdata
    yn = event.ydata

    if xn is None or yn is None:  # mouse is out of region XLIM, YLIM
        plots.reset()
        plots.draw_idle()
        return

    if event.button == 1:
        data.inc()
    elif event.button == 2:  # reset all
        plots.clear()
        data.reset()
    elif event.button == 3:  # reset initial conditions
        data.reset()

    data.update((xn, yn))
    title = ""

    if data.get_click_num() == 0:  # before 1st click
        x0, y0 = data.get_data(0)
        plots.reset()
        plots.draw_dot(x0, y0)
        title = f"y({xn:.2f})={yn:.2f}"
    elif data.get_click_num() == 1:  # before 2nd click
        x0, y0 = data.get_data(0)
        x1, y1, dy = data.get_data(1)
        plots.draw_dot(x0, y0)
        plots.draw_tang(x0, y0, x1, y1)
        title = f"y({x0:.2f})={y0:.2f},  y'({x0:.2f})={dy:.2f}"
    elif data.get_click_num() > 1:  # 2nd and next clicks
        x0, y0 = data.get_data(0)
        x1, y1, dy = data.get_data(1)
        x2, y2, d2y = data.get_data(2)
        plots.draw_dot(x0, y0)
        plots.draw_tang(x0, y0, x1, y1)
        plots.draw_parab(x0, y0, dy, d2y)
        title = f"y({x0:.2f})={y0:.2f},  y'({x0:.2f})={dy:.2f}, y''({x0:.2f})={d2y:.2f}"

        soln = dsolve(func=f, y0=y0, dy0=dy, d2y0=d2y, x0=x0)
        plots.draw_curve(soln)

        if event.button == 1 and data.get_click_num() > 2:
            plots.add_new(soln)  # freeze plot
            print(title)

    plots.set_title(title)
    plots.draw_idle()  # make all drawings in axes


def main():
    axes = InitAxes()
    ax = axes.get_ax()

    data = DataManager()
    plots = Plots(ax)

    callback = lambda event: on_move(event, data, plots, ax)
    axes.set_events(callback)

    plt.show()


if __name__ == "__main__":
    main()
