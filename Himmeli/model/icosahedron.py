import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from ..himmeli import Himmeli, a3_landscape
from ..utils import (
    plot_side_3d,
    plot_surface_3d,
    plot_2d,
    get_carr,
    map_to_2d
)


class Icosahedron(Himmeli):
    papersize = a3_landscape

    def __init__(self,
                 r1: float,
                 r2: float,
                 h1: float,
                 h2: float,
                 h3: float,
                 n: int,
                 folder: Path = Path(".")):
        super().__init__(n=n, folder=folder)
        self.r1 = r1
        self.r2 = r2
        self.h1 = h1
        self.h2 = h2
        self.h3 = h3

        self.layer = [np.array([0, 0, self.h1]),
                      self.radius(self.r1),
                      self.radius(self.r2, h=-self.h2, offset=self.dtheta / 2),
                      np.array([0, 0, -self.h2 - self.h3])]

    def __str__(self) -> str:
        return f"Icosahedron: {self.a1:g}, {self.a2:g}, {self.a3:g}, {self.b1:g}, {self.b2:g}, {self.n}"

    @classmethod
    def from_in_kind(cls,
                     a1: float,
                     a2: float,
                     a3: float,
                     b1: float,
                     b2: float,
                     n: int,
                     folder: Path = Path(".")):
        if a2 is None:
            a2 = a1
        if a3 is None:
            a3 = a1
        b1 = b1
        if b1 is None:
            b1 = a2
        b2 = b2
        if b2 is None:
            b2 = b1

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

        a2_u = np.sqrt(
            r1 ** 2
            + r2 ** 2
            - 2 * r1 * r2 * np.cos(dtheta / 2)
        )
        if a2 <= a2_u:
            print(f"`a2` must large than {a2_u}")
            print(f"`a2`={a2}")
            raise ValueError()

        if a3 <= r2:
            print(f"`a3` must large than {r2}")
            print(f"`a3`={a3}")
            print(f"In other wards,")
            b2_u = a3 * (2 * np.sin(dtheta / 2))
            print(f"`b2` must small than {b2_u}")
            print(f"`b2`={b2}")
            raise ValueError()

        h1 = np.sqrt(a1 ** 2 - r1 ** 2)
        h2 = np.sqrt(
            a2 ** 2
            - r1 ** 2
            - r2 ** 2
            + 2 * r1 * r2 * np.cos(dtheta / 2)
        )
        h3 = np.sqrt(a3 ** 2 - r2 ** 2)

        return cls(r1, r2, h1, h2, h3, n, folder=folder)

    @property
    def a1(self):
        return np.sqrt(self.h1**2 + self.r1**2)

    @property
    def a2(self):
        return np.sqrt(self.h2**2 +
                       self.r1**2 +
                       self.r2**2 -
                       2 * self.r1 * self.r2 * np.cos(self.dtheta / 2))

    @property
    def a3(self):
        return np.sqrt(self.h3**2 + self.r2**2)

    @property
    def b1(self):
        return 2 * self.r1 * np.sin(self.dtheta / 2)

    @property
    def b2(self):
        return 2 * self.r2 * np.sin(self.dtheta / 2)

    @property
    def name(self):
        return f"Icosahedron-{self.a1:g}-{self.a2:g}-{self.a3:g}-{self.b1:g}-{self.b2:g}-{self.n}"

    def plot(self, ax):
        self.plot_side(ax)
        self.plot_surface(ax)
        r = max(self.r1, self.r2)
        ax.set_box_aspect((r, r, (self.h1 + self.h2 + self.h3) / 2))

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

            theta0 = self.dtheta * i
            theta1 = self.dtheta * (i - 0.5)
            x0 = self.r1 * np.cos(theta0), self.r1 * np.sin(theta0), 0
            x1 = self.r2 * np.cos(theta1), self.r2 * np.sin(theta1), -self.h2
            plot_side_3d(ax, *x0, *x1)

            theta0 = self.dtheta * i
            theta1 = self.dtheta * (i + 0.5)
            x0 = self.r1 * np.cos(theta0), self.r1 * np.sin(theta0), 0
            x1 = self.r2 * np.cos(theta1), self.r2 * np.sin(theta1), -self.h2
            plot_side_3d(ax, *x0, *x1)

            theta0 = self.dtheta * (i - 0.5)
            theta1 = self.dtheta * (i + 0.5)

            x0 = self.r2 * np.cos(theta0), self.r2 * np.sin(theta0), -self.h2
            x1 = self.r2 * np.cos(theta1), self.r2 * np.sin(theta1), -self.h2
            plot_side_3d(ax, *x0, *x1)

            x1 = 0, 0, -self.h2 - self.h3
            plot_side_3d(ax, *x0, *x1)

    def plot_surface(self, ax):
        for i in range(self.n):
            theta0 = self.dtheta * i
            theta1 = self.dtheta * (i + 1)

            x0 = self.r1 * np.cos(theta0), self.r1 * np.sin(theta0), 0
            x1 = self.r1 * np.cos(theta1), self.r1 * np.sin(theta1), 0
            x2 = 0, 0, self.h1
            plot_surface_3d(ax, *x0, *x1, *x2, *x2)

            theta2 = self.dtheta * (i + 0.5)
            x2 = self.r2 * np.cos(theta2), self.r2 * np.sin(theta2), -self.h2
            plot_surface_3d(ax, *x0, *x1, *x2, *x2)

            theta0 = self.dtheta * (i - 0.5)
            theta1 = self.dtheta * (i + 0.5)
            x0 = self.r2 * np.cos(theta0), self.r2 * np.sin(theta0), -self.h2
            x1 = self.r2 * np.cos(theta1), self.r2 * np.sin(theta1), -self.h2
            x2 = 0, 0, -self.h2 - self.h3
            plot_surface_3d(ax, *x0, *x1, *x2, *x2)

            theta2 = self.dtheta * i
            x2 = self.r1 * np.cos(theta2), self.r1 * np.sin(theta2), 0
            plot_surface_3d(ax, *x0, *x1, *x2, *x2)

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
            p2 = self.layer[2][i]
            parr = np.array([p0, p1, p2, p0])

            off = plot_2d(ax, parr, p0, p1, p2, off, coeff)
        off_ = off

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

        off = [self.b1, 0]
        coeff = [1, 1]
        if True:
            p0 = self.layer[0]
            p1 = self.layer[1][0]
            p2 = self.layer[2][0]
            parr = np.array([p0, p1, p2])

            off = plot_2d(ax, parr, p1, (p0 + p2) / 2, p0, off, coeff, ls=":")
            off[1] = 0

        off = off_
        off[0] = 0
        coeff = [1, 1]
        for i in range(self.n):
            off[1] = off_[1]
            p0 = self.layer[2][i]
            p1 = self.layer[2][(i - 1) % self.n]
            p2 = self.layer[1][i]
            parr = np.array([p0, p1, p2, p0])

            if i == 0:
                _, yarr = map_to_2d(parr, p0, p1, p2, [0, 0], coeff)
                off_[1] -= yarr[2]
            off = plot_2d(ax, parr, p0, p1, p2, off, coeff)

        off = off_
        coeff = [1, -1]
        for i in range(self.n):
            off[1] = off_[1]
            p0 = self.layer[2][i]
            p1 = self.layer[2][(i + 1) % self.n]
            p2 = self.layer[3]
            parr = np.array([p0, p1, p2, p0])

            off = plot_2d(ax, parr, p0, p1, p2, off, coeff)

        coeff = [-1, 1]
        if True:
            off[1] = off_[1]
            p0 = self.layer[2][0]
            p1 = self.layer[2][1]
            p2 = self.layer[2][2]
            parr = np.vstack([self.layer[2], p0])
            carr = get_carr(p0, p1, p2)

            plot_2d(ax, parr, p0, p1, p2, off, coeff, ls=":")
            plot_2d(ax, carr, p0, p1, p2, off, coeff, ls=":")

        off[0] -= self.b2
        coeff = [-1, -1]
        if True:
            p0 = self.layer[3]
            p1 = self.layer[2][0]
            p2 = self.layer[1][0]
            parr = np.array([p0, p1, p2])

            off = plot_2d(ax, parr, p1, (p0 + p2) / 2, p0, off, coeff, ls=":")
            off[1] = 0

        ax.set_aspect("equal")

        ax.set_xticks([])
        ax.set_yticks([])
