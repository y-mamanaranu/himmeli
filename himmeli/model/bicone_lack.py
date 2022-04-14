from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ..himmeli import Himmeli
from ..utils import get_carr, map_to_2d, plot_2d, plot_side_3d, plot_surface_3d


class Bicone_Lack(Himmeli):
    def __init__(self,
                 r1: float,
                 r2: float,
                 h1: float,
                 h2: float,
                 n: int,
                 folder: Path = Path(".")):
        super().__init__(n=n, folder=folder)
        self.r1 = r1
        self.r2 = r2
        self.h1 = h1
        self.h2 = h2

        self.layer = [np.array([0, 0, self.h1]),
                      self.radius(self.r1),
                      self.radius(self.r2, h=-self.h2)]

    def __str__(self) -> str:
        return f"Bicone_Lack: {self.a1:g}, {self.a2:g}, {self.b1:g}, {self.b2:g}, {self.n}"

    @classmethod
    def from_in_kind(cls,
                     a1: float,
                     a2: float,
                     b1: float,
                     b2: float,
                     n: int,
                     folder: Path = Path(".")):
        if b1 is None:
            b1 = a1
        if a2 is None and b2 is None:
            print(f"`a2` or `b2` must be not `None`")
            raise ValueError()
        if a2 is None:
            a2 = a1 * (b1 - b2) / b1
        elif b2 is None:
            b2 = b1 * (a1 - a2) / a1

        if a1 <= a2:
            print(f"`a1` must large than `a2`")
            print(f"`a1`={a1}")
            print(f"`a2`={a2}")
            raise ValueError()

        if b1 <= b2:
            print(f"`b1` must large than `b2`")
            print(f"`b1`={b1}")
            print(f"`b2`={b2}")
            raise ValueError()

        dtheta = 2 * np.pi / n
        r1 = b1 / (2 * np.sin(dtheta / 2))
        r2 = b2 / (2 * np.sin(dtheta / 2))

        if a1 <= r1:
            print(f"`a1` must large than {r1}")
            print(f"`a1`={a1}")
            print(f"In other wards,")
            b1_u = a1 * (2 * np.sin(dtheta / 2))
            print(f"`b1` must small than {b1_u}")
            print(f"`b1`={b1}")
            raise ValueError()

        if a2 <= r1 - r2:
            print(f"`a2` must large than {r1-r2}")
            print(f"`a2`={a2}")
            print(f"In other wards,")
            b2_u = (r1 - a2) * (2 * np.sin(dtheta / 2))
            print(f"`b2` must large than {b2_u}")
            print(f"`b2`={b2}")
            raise ValueError()

        h1 = np.sqrt(a1**2 - r1**2)
        h2 = np.sqrt(a2**2 - (r1 - r2)**2)

        return cls(r1, r2, h1, h2, n, folder=folder)

    @property
    def a1(self):
        return np.sqrt(self.h1**2 + self.r1**2)

    @property
    def a2(self):
        return np.sqrt(self.h2**2 + (self.r1 - self.r2)**2)

    @property
    def b1(self):
        return 2 * self.r1 * np.sin(self.dtheta / 2)

    @property
    def b2(self):
        return 2 * self.r2 * np.sin(self.dtheta / 2)

    @property
    def name(self):
        return f"Bicone_Lack-{self.a1:g}-{self.a2:g}-{self.b1:g}-{self.b2:g}-{self.n}"

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

    def plot_expansion(self, ax: plt.Axes):
        off = [0, 0]
        coeff = [1, 1]
        for i in range(self.n):
            off[1] = 0
            p0 = self.layer[1][i]
            p1 = self.layer[1][(i + 1) % self.n]
            p2 = self.layer[0]
            parr = np.array([p0, p1, p2, p0])

            off = plot_2d(ax, parr, p0, p1, p2, off, coeff)

        off = [0, 0]
        coeff = [1, -1]
        for i in range(self.n):
            off[1] = 0
            p0 = self.layer[1][i]
            p1 = self.layer[1][(i + 1) % self.n]
            p2 = self.layer[2][(i + 1) % self.n]
            p3 = self.layer[2][i]
            parr = np.array([p0, p1, p2, p3, p0])

            if i == 0:
                xarr, yarr = map_to_2d(parr, p0, p1, p2, off, coeff)
                off_ = [xarr[3], yarr[3]]
            off = plot_2d(ax, parr, p0, p1, p2, off, coeff)

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

        off = off_
        coeff = [1, -1]
        if True:
            p0 = self.layer[2][0]
            p1 = self.layer[2][1]
            p2 = self.layer[2][2]
            parr = np.vstack([self.layer[2], p0])
            carr = get_carr(p0, p1, p2)

            plot_2d(ax, parr, p0, p1, p2, off, coeff)
            plot_2d(ax, carr, p0, p1, p2, off, coeff, ls=":")

        off = [self.b1, 0]
        coeff = [1, 1]
        if True:
            p0 = self.layer[0]
            p1 = self.layer[1][0]
            p2 = self.layer[2][0]
            p3 = np.array([0, 0, -self.h2])
            parr = np.array([p0, p1, p2, p3])

            off = plot_2d(ax, parr, p1, self.C, p0, off, coeff, ls=":")
            off[1] = 0

        ax.set_aspect("equal")

        ax.set_xticks([])
        ax.set_yticks([])
