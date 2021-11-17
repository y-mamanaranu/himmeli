import numpy as np
from pathlib import Path

from ..himmeli import Himmeli
from ..utils import (
    plot_circle,
    plot_isosceles,
    plot_polygon,
    plot_side_3d,
    plot_surface_3d,
)


class Cone(Himmeli):
    def __init__(self, a, b, n, folder=Path(".")):
        super().__init__(folder=folder)

        self.a = a
        self.b = b
        if b is None:
            self.b = self.a
        self.n = n
        self.dtheta = 2 * np.pi / self.n
        self.r = self.b / (2 * np.sin(self.dtheta / 2))

        if self.a <= self.r:
            print(f"`a` must large than {self.r}")
            print(f"`a`={self.a}")
            print(f"In other wards,")
            b_u = self.a * (2 * np.sin(self.dtheta / 2))
            print(f"`b` must small than {b_u}")
            print(f"`b`={self.b}")
            raise ValueError()

        self.h = np.sqrt(self.a**2 - self.r**2)

    def __str__(self) -> str:
        return f"Cone: {self.a}, {self.b}, {self.n}"

    @property
    def name(self):
        return f"Cone-{self.a}-{self.b}-{self.n}"

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
        w = self.b
        h = np.sqrt(self.a**2 - (self.b / 2)**2)
        for i in range(self.n):
            x0 = self.b * i
            plot_isosceles(ax, x0, 0, w, h)

        plot_circle(ax, 0, 0, self.r, np.pi / 2 + self.dtheta / 2)
        plot_polygon(ax, 0, 0, self.r, self.n, np.pi / 2 + self.dtheta / 2)

        ax.set_aspect("equal")

        ax.set_xticks([])
        ax.set_yticks([])
