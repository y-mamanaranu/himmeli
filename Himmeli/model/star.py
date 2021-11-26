import numpy as np
from pathlib import Path

from ..himmeli import Himmeli, a2_landscape
from ..utils import (
    plot_circle,
    plot_isosceles,
    plot_polygon,
    plot_side_3d,
    plot_surface_3d,
    plot_side,
    radius
)


class Star(Himmeli):
    papersize = a2_landscape

    def __init__(self, a1, a2, b, n, folder=Path(".")):
        super().__init__(folder=folder)
        self.a1 = a1
        self.a2 = a2
        if a2 is None:
            self.a2 = self.a1
        self.b = b
        if b is None:
            self.b = self.a1
        self.n = n
        self.dtheta = 2 * np.pi / self.n
        self.r1 = self.b / (2 * np.sin(self.dtheta / 2))
        self.r2 = self.r1 + np.sqrt(self.a2 ** 2 - (self.b / 2) ** 2)

        self.h1 = np.sqrt(self.a1 ** 2 - self.r1 ** 2)

        self.box = (self.r2, self.r2, self.h1)

    def __str__(self) -> str:
        return f"Star: {self.a1}, {self.a2}, {self.b}, {self.n}"

    @property
    def name(self):
        return f"Star-{self.a1}-{self.a2}-{self.b}-{self.n}"

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

    def plot_expansion(self, ax):
        w = self.b
        h = np.sqrt(self.a1 ** 2 - (self.b / 2) ** 2)
        for i in range(self.n):
            x0 = self.b * i
            plot_isosceles(ax, x0, 0, w, h)
            plot_isosceles(ax, x0 + w, 0, -w, -h)

        a3 = np.sqrt(self.r2 ** 2 + self.h1 ** 2)
        tmp_1 = self.a1 ** 2 + self.a2 ** 2 - a3 ** 2
        tmp_2 = 2 * self.a1 * self.a2
        cos_12 = np.abs(tmp_1 / tmp_2)
        sin_12 = np.sqrt(1 - cos_12 ** 2)
        sin_23 = self.a1 * sin_12 / a3
        cos_23 = np.sqrt(1 - sin_23 ** 2)
        sin_2_23 = 2 * sin_23 * cos_23
        cos_2_23 = cos_23 ** 2 - sin_23 ** 2

        self.qoff_arr = [[2 * self.a2 * i + w / 2, h + a3]
                         for i in range(self.n)]

        for qoff in self.qoff_arr:
            q0 = np.array([0, 0])
            q1 = np.array([-self.a2 * cos_2_23, self.a2 * sin_2_23])
            q2 = np.array([-a3 * cos_23, a3 * sin_23])
            q3 = np.array([-self.a2, 0])
            for c in [-1, 1]:
                plot_side(ax, *q0, *q1, *qoff, c=c, linestyle=":")
                plot_side(ax, *q0, *q2, *qoff, c=c)
                plot_side(ax, *q0, *q3, *qoff, c=c)
                plot_side(ax, *q1, *q2, *qoff, c=c, linestyle=":")
                plot_side(ax, *q3, *q2, *qoff, c=c, linestyle=":")

        plot_circle(ax, 0, 0, self.r1, np.pi / 2 + self.dtheta / 2)
        plot_polygon(ax, 0, 0, self.r1, self.n, np.pi / 2 + self.dtheta / 2)

        ax.set_aspect("equal")

        ax.set_xticks([])
        ax.set_yticks([])
