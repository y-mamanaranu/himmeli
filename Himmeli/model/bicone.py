import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from ..himmeli import Himmeli
from ..utils import (
    plot_side_3d,
    plot_surface_3d,
    plot_2d,
    get_carr,
)


class Bicone(Himmeli):
    def __init__(self,
                 r: float,
                 h1: float,
                 h2: float,
                 n: int,
                 folder: Path = Path(".")):
        super().__init__(n=n, folder=folder)
        self.r = r
        self.h1 = h1
        self.h2 = h2

        self.layer = [np.array([0, 0, self.h1]),
                      self.radius(self.r),
                      np.array([0, 0, -self.h2])]

    def __str__(self) -> str:
        return f"Bicone: {self.a1:g}, {self.a2:g}, {self.b:g}, {self.n}"

    @classmethod
    def from_in_kind(cls,
                     a1: float,
                     a2: float,
                     b: float,
                     n: int,
                     folder: Path = Path(".")):
        if a2 is None:
            a2 = a1
        b = b
        if b is None:
            b = a1

        n = n
        dtheta = 2 * np.pi / n
        r = b / (2 * np.sin(dtheta / 2))

        if a1 <= r:
            print(f"`a1` must large than {r}")
            print(f"`a1`={a1}")
            print(f"In other wards,")
            b_u = a1 * (2 * np.sin(dtheta / 2))
            print(f"`b` must small than {b_u}")
            print(f"`b`={b}")
            raise ValueError()

        if a2 <= r:
            print(f"`a2` must large than {r}")
            print(f"`a2`={a2}")
            print(f"In other wards,")
            b_u = a2 * (2 * np.sin(dtheta / 2))
            print(f"`b` must small than {b_u}")
            print(f"`b`={b}")
            raise ValueError()

        h1 = np.sqrt(a1**2 - r**2)
        h2 = np.sqrt(a2**2 - r**2)

        return cls(r, h1, h2, n, folder=folder)

    @property
    def a1(self):
        return np.sqrt(self.h1**2 + self.r**2)

    @property
    def a2(self):
        return np.sqrt(self.h2**2 + self.r**2)

    @property
    def b(self):
        return 2 * self.r * np.sin(self.dtheta / 2)

    @property
    def name(self):
        return f"Bicone-{self.a1:g}-{self.a2:g}-{self.b:g}-{self.n}"

    def plot(self, ax):
        self.plot_side(ax)
        self.plot_surface(ax)
        ax.set_box_aspect((self.r, self.r, (self.h1 + self.h2) / 2))

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

    def plot_side(self, ax):
        for i in range(self.n):
            theta0 = self.dtheta * i
            theta1 = self.dtheta * (i + 1)

            x0 = self.r * np.cos(theta0), self.r * np.sin(theta0), 0
            x1 = 0, 0, self.h1
            plot_side_3d(ax, *x0, *x1)

            x1 = self.r * np.cos(theta1), self.r * np.sin(theta1), 0
            plot_side_3d(ax, *x0, *x1)

            x1 = 0, 0, -self.h2
            plot_side_3d(ax, *x0, *x1)

    def plot_surface(self, ax):
        for i in range(self.n):
            theta0 = self.dtheta * i
            theta1 = self.dtheta * (i + 1)

            x0 = self.r * np.cos(theta0), self.r * np.sin(theta0), 0
            x1 = self.r * np.cos(theta1), self.r * np.sin(theta1), 0
            for x2 in ([
                [0, 0, self.h1],
                [0, 0, -self.h2]
            ]):
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
        for i in range(self.n):
            p0 = self.layer[1][i]
            p1 = self.layer[1][(i + 1) % self.n]
            p2 = self.layer[2]
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

            plot_2d(ax, parr, p0, p1, p2, off, coeff, ls=":")
            plot_2d(ax, carr, p0, p1, p2, off, coeff, ls=":")

        off = [self.b, 0]
        coeff = [1, 1]
        if True:
            p0 = self.layer[0]
            p1 = self.layer[1][0]
            p2 = self.layer[2]
            parr = np.array([p0, p1, p2])

            off = plot_2d(ax, parr, p1, self.C, p0, off, coeff, ls=":")
            off[1] = 0

        ax.set_aspect("equal")

        ax.set_xticks([])
        ax.set_yticks([])
