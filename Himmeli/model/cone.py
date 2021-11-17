from ..utils import plot_side_3d, plot_surface_3d
import numpy as np


class Cone(object):
    def __init__(self, a, b, n):
        self.a = a
        self.b = b
        self.n = n
        self.dtheta = 2 * np.pi / self.n
        self.r = self.b / (2 * np.sin(self.dtheta / 2))

        if self.a <= self.r:
            print(f"`a` must large than `r`")
            print(f"`a`={self.a}")
            print(f"`r`={self.r}")
            raise ValueError()

        self.h = np.sqrt(self.a**2 - self.r**2)

    def plot(self, ax):
        self.plot_side(ax)
        self.plot_surface(ax)
        ax.set_box_aspect((self.r, self.r, self.h / 2))

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

    def plot_side(self, ax):
        for i in range(self.n):
            theta0 = self.dtheta * i
            theta1 = self.dtheta * (i + 1)

            x0 = self.r * np.cos(theta0), self.r * np.sin(theta0), 0
            x1 = 0, 0, self.h
            plot_side_3d(ax, *x0, *x1)

            x1 = self.r * np.cos(theta1), self.r * np.sin(theta1), 0
            plot_side_3d(ax, *x0, *x1)

    def plot_surface(self, ax):
        for i in range(self.n):
            theta0 = self.dtheta * i
            theta1 = self.dtheta * (i + 1)

            x0 = self.r * np.cos(theta0), self.r * np.sin(theta0), 0
            x1 = self.r * np.cos(theta1), self.r * np.sin(theta1), 0
            x2 = 0, 0, self.h
            plot_surface_3d(ax, *x0, *x1, *x2, *x2)

    def plot_expansion(self, ax):
        for i in range(self.n):
            x0 = self.b * i
            x1 = self.b * (i + 1)
            x2 = self.b * (i + 0.5)
            y0 = 0
            y1 = 0
            y2 = np.sqrt(self.a**2 - (self.b / 2)**2)

            ax.plot([x0, x1, x2, x0], [y0, y1, y2, y0], c="C0")

        theta = np.linspace(0, 2 * np.pi, self.n + 1, endpoint=True)
        x = self.b / 2 + self.r * -np.sin(theta + self.dtheta / 2)
        y = -self.r * np.cos(self.dtheta / 2) + self.r * \
            np.cos(theta + self.dtheta / 2)
        ax.plot(x, y, c="C0", linestyle="dashed")

        theta = np.linspace(0, 2 * np.pi, 100, endpoint=True)
        x = self.b / 2 + self.r * -np.sin(theta + self.dtheta / 2)
        y = -self.r * np.cos(self.dtheta / 2) + self.r * \
            np.cos(theta + self.dtheta / 2)

        ax.plot(x, y, c="C0", linestyle="dashed")

        ax.set_aspect("equal")

        ax.set_xticks([])
        ax.set_yticks([])
