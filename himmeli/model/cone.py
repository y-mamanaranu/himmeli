from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ..himmeli import Himmeli
from ..utils import get_carr, plot_2d, plot_side_3d, plot_surface_3d


class Cone(Himmeli):
    def __init__(self,
                 r: float,
                 h: float,
                 n: int,
                 folder: Path = Path(".")):
        super().__init__(n=n, folder=folder)
        self.r = r
        self.h = h

        self.layer = [np.array([0, 0, self.h]),
                      self.radius(self.r)]

    def __str__(self) -> str:
        return f"Cone: {self.a:g}, {self.b:g}, {self.n}"

    @classmethod
    def from_in_kind(cls,
                     a: float,
                     b: float,
                     n: int,
                     folder: Path = Path(".")):
        if b is None:
            b = a

        dtheta = 2 * np.pi / n
        r = b / (2 * np.sin(dtheta / 2))

        if a <= r:
            print(f"`a` must large than {r}")
            print(f"`a`={a}")
            print(f"In other wards,")
            b_u = a * (2 * np.sin(dtheta / 2))
            print(f"`b` must small than {b_u}")
            print(f"`b`={b}")
            raise ValueError()

        h = np.sqrt(a**2 - r**2)

        return cls(r, h, n, folder=folder)

    @property
    def a(self):
        return np.sqrt(self.h**2 + self.r**2)

    @property
    def b(self):
        return 2 * self.r * np.sin(self.dtheta / 2)

    @property
    def name(self):
        return f"Cone-{self.a:g}-{self.b:g}-{self.n}"

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

    def plot_expansion(self, ax: plt.Axes):
        off = [0, 0]
        coeff = [1, 1]
        for i in range(self.n):
            p0 = self.layer[1][i]
            p1 = self.layer[1][(i + 1) % self.n]
            p2 = self.layer[0]
            parr = np.array([p0, p1, p2, p0])

            off = plot_2d(ax, parr, p0, p1, p2, off, coeff)
            off[1] = 0

        off = [0, 0]
        coeff = [1, -1]
        if True:
            p0 = self.layer[1][0]
            p1 = self.layer[1][1]
            p2 = self.layer[1][2]
            parr = np.vstack([self.layer[1], p0])
            carr = get_carr(p0, p1, p2)

            plot_2d(ax, parr, p0, p1, p2, off, coeff)
            plot_2d(ax, carr, p0, p1, p2, off, coeff, ls=":")

        off = [self.b, 0]
        coeff = [1, 1]
        if True:
            p0 = self.layer[0]
            p1 = self.layer[1][0]
            p2 = np.array([0, 0, 0])
            parr = np.array([p0, p1])

            off = plot_2d(ax, parr, p1, p2, p0, off, coeff, ls=":")
            off[1] = 0

        ax.set_aspect("equal")

        ax.set_xticks([])
        ax.set_yticks([])
