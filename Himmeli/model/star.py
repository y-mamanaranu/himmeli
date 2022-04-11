import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from ..himmeli import Himmeli, a1_landscape
from ..utils import (
    plot_side_3d,
    plot_surface_3d,
    radius,
    plot_2d,
    get_carr,
    map_to_2d,
)


class Star(Himmeli):
    papersize = a1_landscape

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

        self.box = (self.r2, self.r2, (self.h1 + self.h2) / 2)

        self.layer = [np.array([0, 0, self.h1]),
                      self.radius(self.r1),
                      self.radius(self.r2, offset=self.dtheta / 2),
                      np.array([0, 0, -self.h2])]

    def __str__(self) -> str:
        return f"Star: {self.a1:g}, {self.a2:g}, {self.a3:g}, {self.b:g}, {self.n}"

    @classmethod
    def from_in_kind(cls,
                     a1: float,
                     a2: float,
                     a3: float,
                     b: float,
                     n: int,
                     folder: Path = Path(".")):
        if a2 is None:
            a2 = a1
        if a3 is None:
            a3 = a1
        if b is None:
            b = a1

        dtheta = 2 * np.pi / n
        r1 = b / (2 * np.sin(dtheta / 2))
        r2 = r1 + np.sqrt(a2 ** 2 - (b / 2) ** 2)

        h1 = np.sqrt(a1 ** 2 - r1 ** 2)
        h2 = np.sqrt(a3 ** 2 - r1 ** 2)

        return cls(r1, r2, h1, h2, n, folder=folder)

    @property
    def a1(self):
        return np.sqrt(self.h1**2 + self.r1**2)

    @property
    def a2(self):
        return np.sqrt((self.r2 - self.r1)**2 + (self.b / 2) ** 2)

    @property
    def a3(self):
        return np.sqrt(self.h2**2 + self.r1**2)

    @property
    def b(self):
        return 2 * self.r1 * np.sin(self.dtheta / 2)

    @property
    def b_star(self):
        return np.sqrt(self.r1**2 + self.r2**2 - 2 * self.r1 *
                       self.r2 * np.cos(self.dtheta / 2))

    @property
    def name(self):
        return f"Star-{self.a1:g}-{self.a2:g}-{self.a3:g}-{self.b:g}-{self.n}"

    def plot(self, ax):
        theta = np.linspace(0, 2 * np.pi, self.n, endpoint=False)
        self.p1 = [0, 0, self.h1]
        self.p2_arr = self.r1 * radius(theta)
        self.p2_brr = np.roll(self.p2_arr, 1, axis=0)
        self.p3_arr = self.r2 * radius(theta, self.dtheta / 2)
        self.p3_brr = np.roll(self.p3_arr, 1, axis=0)
        self.p4 = [0, 0, -self.h1]

        self.plot_side(ax)
        self.plot_surface(ax)
        ax.set_box_aspect(self.box)

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

    def plot_side(self, ax):
        for p2 in self.p2_arr:
            plot_side_3d(ax, *self.p1, *p2)
            plot_side_3d(ax, *self.p4, *p2)

        for p2_1, p2_2 in zip(self.p2_arr, self.p2_brr):
            plot_side_3d(ax, *p2_1, *p2_2)

        for p2, p3_1, p3_2 in zip(self.p2_arr, self.p3_arr, self.p3_brr):
            plot_side_3d(ax, *p2, *p3_1)
            plot_side_3d(ax, *p2, *p3_2)

        for p3 in self.p3_arr:
            plot_side_3d(ax, *self.p1, *p3)
            plot_side_3d(ax, *self.p4, *p3)

    def plot_surface(self, ax):
        for p2, p3_1, p3_2 in zip(self.p2_arr, self.p3_arr, self.p3_brr):
            plot_surface_3d(ax, *self.p1, *p2, *p3_1, *p3_1)
            plot_surface_3d(ax, *self.p1, *p2, *p3_2, *p3_2)

            plot_surface_3d(ax, *self.p4, *p2, *p3_1, *p3_1)
            plot_surface_3d(ax, *self.p4, *p2, *p3_2, *p3_2)

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
            p2 = self.layer[3]
            parr = np.array([p0, p1, p2, p0])

            off = plot_2d(ax, parr, p0, p1, p2, off, coeff)

        off = [0, 0]
        coeff = [1, -1]
        if True:
            p0 = self.layer[1][0]
            p1 = self.layer[1][1]
            p2 = self.layer[1][2]
            parr = np.vstack([self.layer[1], p0])
            iarr = 1 + np.arange(self.n)
            parr2 = np.insert(parr, iarr, self.layer[2], axis=0)
            carr = get_carr(p0, p1, p2)

            plot_2d(ax, parr, p0, p1, p2, off, coeff, ls=":")
            off1_ = plot_2d(ax, parr2, p0, p1, p2, off, coeff, ls=":")
            off2_ = plot_2d(ax, carr, p0, p1, p2, off, coeff, ls=":")

        off = off_ = [self.b / 2, min(off1_[1], off2_[1])]
        for i in range(self.n):
            coeff = [-1, 1]
            p0 = self.layer[2][i]
            p1 = self.layer[1][i]
            p2 = self.layer[0]
            parr = np.array([p0, p1, p2, p0])

            if i == 0:
                xarr, yarr = map_to_2d(parr, p0, p1, p2, [0, 0], coeff)
                yarr[1] = yarr[0] + yarr[2] - yarr[1]

                off_[1] -= max(yarr[1], yarr[2])

            xarr, yarr = map_to_2d(parr, p0, p1, p2, off, coeff)
            ax.plot(xarr, yarr, c="C0")

            xarr[1] = xarr[0] + xarr[2] - xarr[1]
            yarr[1] = yarr[0] + yarr[2] - yarr[1]
            ax.plot(xarr, yarr, c="C0")

            coeff = [1, -1]
            p0 = self.layer[2][i]
            p1 = self.layer[1][i]
            p2 = self.layer[3]
            parr = np.array([p0, p1, p2, p0])

            xarr, yarr = map_to_2d(parr, p0, p1, p2, off, coeff)
            ax.plot(xarr, yarr, c="C0")

            xarr[1] = xarr[0] + xarr[2] - xarr[1]
            yarr[1] = yarr[0] + yarr[2] - yarr[1]
            ax.plot(xarr, yarr, c="C0")

            off[0] += 2 * self.b_star

        ax.set_aspect("equal")

        ax.set_xticks([])
        ax.set_yticks([])
