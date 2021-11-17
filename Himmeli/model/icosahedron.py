import numpy as np
from pathlib import Path

from ..himmeli import Himmeli, a3_landscape
from ..utils import (
    plot_circle,
    plot_isosceles,
    plot_polygon,
    plot_side_3d,
    plot_surface_3d,
)


class Icosahedron(Himmeli):
    papersize = a3_landscape

    def __init__(self, a1, a2, a3, b1, b2, n, folder=Path(".")):
        super().__init__(folder=folder)

        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        if a2 is None:
            self.a2 = self.a1
        if a3 is None:
            self.a3 = self.a1
        self.b1 = b1
        if b1 is None:
            self.b1 = self.a2
        self.b2 = b2
        if b2 is None:
            self.b2 = self.b1
        self.n = n
        self.dtheta = 2 * np.pi / self.n
        self.r1 = self.b1 / (2 * np.sin(self.dtheta / 2))
        self.r2 = self.b2 / (2 * np.sin(self.dtheta / 2))

        if self.a1 <= self.r1:
            print(f"`a1` must large than {self.r1}")
            print(f"`a1`={self.a1}")
            print(f"In other wards,")
            b1_u = self.a1 * (2 * np.sin(self.dtheta / 2))
            print(f"`b1` must small than {b_u}")
            print(f"`b1`={self.b1}")
            raise ValueError()

        a2_u = np.sqrt(
            self.r1 ** 2
            + self.r2 ** 2
            - 2 * self.r1 * self.r2 * np.cos(self.dtheta / 2)
        )
        if self.a2 <= a2_u:
            print(f"`a2` must large than {a2_u}")
            print(f"`a2`={self.a2}")
            raise ValueError()

        if self.a3 <= self.r2:
            print(f"`a2` must large than {self.r2}")
            print(f"`a2`={self.a2}")
            print(f"In other wards,")
            b2_u = self.a2 * (2 * np.sin(self.dtheta / 2))
            print(f"`b2` must small than {b2_u}")
            print(f"`b2`={self.b2}")
            raise ValueError()

        self.h1 = np.sqrt(self.a1 ** 2 - self.r1 ** 2)
        self.h2 = np.sqrt(
            self.a2 ** 2
            - self.r1 ** 2
            - self.r2 ** 2
            + 2 * self.r1 * self.r2 * np.cos(self.dtheta / 2)
        )
        self.h3 = np.sqrt(self.a3 ** 2 - self.r2 ** 2)

    def __str__(self) -> str:
        return f"Icosahedron: {self.a1}, {self.a2}, {self.a3}, {self.b1}, {self.b2}, {self.n}"

    @property
    def name(self):
        return f"Icosahedron-{self.a1}-{self.a2}-{self.a3}-{self.b1}-{self.b2}-{self.n}"

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

    def plot_expansion(self, ax):
        offset = (
            max(
                2 * self.r1 + np.sqrt(self.a2 ** 2 - (self.b2 / 2) ** 2),
                2 * self.r2 + np.sqrt(self.a2 ** 2 - (self.b1 / 2) ** 2),
            )
            / 2
        )

        w = self.b1
        h = np.sqrt(self.a1 ** 2 - (self.b1 / 2) ** 2)
        for i in range(self.n):
            x0 = self.b1 * i
            plot_isosceles(ax, x0, offset, w, h)

        plot_circle(ax, 0, offset, self.r1, np.pi / 2 + self.dtheta / 2)
        plot_polygon(
            ax,
            0,
            offset,
            self.r1,
            self.n,
            np.pi /
            2 +
            self.dtheta /
            2)

        h = np.sqrt(self.a2 ** 2 - (self.b1 / 2) ** 2)
        for i in range(self.n):
            x0 = self.b1 * i
            plot_isosceles(ax, x0, offset, w, -h)

        w = self.b2
        h = np.sqrt(self.a3 ** 2 - (self.b2 / 2) ** 2)
        for i in range(self.n):
            x0 = self.b2 * i
            plot_isosceles(ax, x0, -offset, w, -h)

        plot_circle(
            ax, self.b2 * (self.n - 1), -offset, self.r2, -
            (np.pi / 2 + self.dtheta / 2)
        )
        plot_polygon(
            ax,
            self.b2 * (self.n - 1),
            -offset,
            self.r2,
            self.n,
            -(np.pi / 2 + self.dtheta / 2),
        )

        h = np.sqrt(self.a2 ** 2 - (self.b2 / 2) ** 2)
        for i in range(self.n):
            x0 = self.b2 * i
            plot_isosceles(ax, x0, -offset, w, h)

        ax.set_aspect("equal")

        ax.set_xticks([])
        ax.set_yticks([])
