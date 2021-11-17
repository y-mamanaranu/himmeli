from ..utils import plot_side_3d, plot_surface_3d
import numpy as np


class Bicone_Lack(object):
    def __init__(self, a, b1, b2, n):
        self.a = a
        self.b1 = b1
        self.b2 = b2
        if self.b1 <= self.b2:
            print(f"`b1` must large than `b2`")
            print(f"`b1`={self.b1}")
            print(f"`b2`={self.b2}")
            raise ValueError()

        self.n = n
        self.dtheta = 2 * np.pi / self.n
        self.r1 = self.b1 / (2 * np.sin(self.dtheta / 2))
        self.r2 = self.r1 * (self.b1 - self.b2) / self.b1

        if self.a <= self.r1:
            print(f"`a` must large than `r1`")
            print(f"`a`={self.a}")
            print(f"`r1`={self.r1}")
            raise ValueError()

        self.h1 = np.sqrt(self.a**2 - self.r1**2)
        self.h2 = self.h1 * self.b2 / self.b1

    def plot(self, ax):
        self.plot_side(ax)
        self.plot_surface(ax)
        ax.set_box_aspect((self.r1, self.r1, (self.h1 + self.h2) / 2))

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

    def plot_side(self, ax):
        for i in range(self.n):
            theta0 = self.dtheta * i
            theta1 = self.dtheta * (i + 1)

            x0 = self.r1 * np.cos(theta0), self.r1 * np.sin(theta0), 0
            x1 = 0, 0, self.h1
            plot_side_3d(ax, *x0, *x1)

            x1 = self.r1 * np.cos(theta1), self.r1 * np.sin(theta1), 0
            plot_side_3d(ax, *x0, *x1)

            x1 = self.r2 * np.cos(theta0), self.r2 * np.sin(theta0), -self.h2
            plot_side_3d(ax, *x0, *x1)

            x0 = self.r2 * np.cos(theta1), self.r2 * np.sin(theta1), -self.h2
            plot_side_3d(ax, *x0, *x1)

    def plot_surface(self, ax):
        for i in range(self.n):
            theta0 = self.dtheta * i
            theta1 = self.dtheta * (i + 1)

            x0 = self.r1 * np.cos(theta0), self.r1 * np.sin(theta0), 0
            x1 = self.r1 * np.cos(theta1), self.r1 * np.sin(theta1), 0
            x2 = 0, 0, self.h1
            plot_surface_3d(ax, *x0, *x1, *x2, *x2)

            x2 = self.r2 * np.cos(theta1), self.r2 * np.sin(theta1), -self.h2
            x3 = self.r2 * np.cos(theta0), self.r2 * np.sin(theta0), -self.h2
            plot_surface_3d(ax, *x0, *x1, *x2, *x3)

    def plot_expansion(self, ax):
        for i in range(self.n):
            x0 = self.b1 * i
            x1 = self.b1 * (i + 1)
            x2 = self.b1 * (i + 0.5)
            y0 = 0
            y1 = 0
            y2 = np.sqrt(self.a**2 - (self.b1 / 2)**2)

            ax.plot([x0, x1, x2, x0], [y0, y1, y2, y0], c="C0")

        for i in range(self.n):
            x0 = self.b1 * i
            x1 = self.b1 * (i + 1)
            x2 = self.b1 * (i + 1 - 0.5 * self.b2 / self.b1)
            x3 = self.b1 * (i + 0.5 * self.b2 / self.b1)
            y0 = 0
            y1 = 0
            y2 = -np.sqrt(self.a**2 - (self.b1 / 2)**2) * self.b2 / self.b1
            y3 = -np.sqrt(self.a**2 - (self.b1 / 2)**2) * self.b2 / self.b1

            ax.plot([x0, x1, x2, x3, x0], [y0, y1, y2, y3, y0], c="C0")

        theta = np.linspace(0, 2 * np.pi, self.n + 1, endpoint=True)
        x = self.b1 / 2 + self.r1 * -np.sin(theta + self.dtheta / 2)
        y = -self.r1 * np.cos(self.dtheta / 2) + self.r1 * \
            np.cos(theta + self.dtheta / 2)

        ax.plot(x, y, c="C0", linestyle="dashed")

        theta = np.linspace(0, 2 * np.pi, 100, endpoint=True)
        x = self.b1 / 2 + self.r1 * -np.sin(theta + self.dtheta / 2)
        y = -self.r1 * np.cos(self.dtheta / 2) + self.r1 * \
            np.cos(theta + self.dtheta / 2)

        ax.plot(x, y, c="C0", linestyle="dashed")

        theta = np.linspace(0, 2 * np.pi, self.n + 1, endpoint=True)
        x = self.b1 / 2 + self.r2 * -np.sin(theta + self.dtheta / 2)
        y = -np.sqrt(self.a**2 - (self.b1 / 2)**2) * self.b2 / self.b1 - self.r2 * \
            np.cos(self.dtheta / 2) + self.r2 * np.cos(theta + self.dtheta / 2)

        ax.plot(x, y, c="C0", linestyle="dashed")

        theta = np.linspace(0, 2 * np.pi, 100, endpoint=True)
        x = self.b1 / 2 + self.r2 * -np.sin(theta + self.dtheta / 2)
        y = -np.sqrt(self.a**2 - (self.b1 / 2)**2) * self.b2 / self.b1 - self.r2 * \
            np.cos(self.dtheta / 2) + self.r2 * np.cos(theta + self.dtheta / 2)

        ax.plot(x, y, c="C0", linestyle="dashed")

        ax.set_aspect("equal")

        ax.set_xticks([])
        ax.set_yticks([])
