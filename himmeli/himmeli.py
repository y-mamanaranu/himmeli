import warnings
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.units import imperial

a1_landscape = [841, 594] * u.mm
a1_portrait = [594, 841] * u.mm

a2_landscape = [594, 420] * u.mm
a2_portrait = [420, 594] * u.mm

a3_landscape = [420, 297] * u.mm
a3_portrait = [297, 420] * u.mm

a4_landscape = [297, 210] * u.mm
a4_portrait = [210, 297] * u.mm


class Himmeli(object):
    papersize = a4_landscape

    def __init__(self, n, folder=Path(".")):
        self.n = n
        self.folder = folder

    @property
    def C(self):
        return np.array([0, 0, 0])

    @property
    def dtheta(self):
        return 2 * np.pi / self.n

    @property
    def theta(self):
        return np.linspace(0, 2 * np.pi, self.n, endpoint=False)

    @property
    def cos(self):
        return np.cos(self.theta)

    @property
    def sin(self):
        return np.sin(self.theta)

    @property
    def zeros(self):
        return np.zeros(self.n)

    def radius(self, r: float, h: float = 0, offset: float = 0, ):
        return np.array([r * np.cos(self.theta + offset),
                         r * np.sin(self.theta + offset),
                         self.zeros + h]).T

    @property
    def circle(self):
        return np.linspace(0, 2 * np.pi, 100, endpoint=True)

    @property
    def name(self):
        return f"Himmeli"

    @property
    def file_2d(self):
        return self.folder / f"{self.name}-2D.pdf"

    @property
    def file_3d(self):
        return self.folder / f"{self.name}-3D.png"

    def plot_expansion():
        warnings.warn("`plot_expansion` is not defined.")

    def plot_paper(
        self, papersize=None, horizontalmargin=0.02, verticalmargin=0.02
    ):
        if papersize is None:
            papersize = self.papersize
        figsize = papersize.to(imperial.inch).value
        fig, ax = plt.subplots(figsize=figsize)

        self.plot_expansion(ax)

        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()

        dx = x1 - x0
        dy = y1 - y0

        w0, h0 = papersize.to(u.mm).value
        w = w0 * (1 - 2 * horizontalmargin)
        h = h0 * (1 - 2 * verticalmargin)

        xblank = (w - dx) / 2
        yblank = (h - dy) / 2

        ax.set_xlim(x0 - xblank, x1 + xblank)
        ax.set_ylim(y0 - yblank, y1 + yblank)

        ax.text(
            x0 - xblank,
            y1 + yblank,
            str(self),
            size=30,
            horizontalalignment="left",
            verticalalignment="top",
        )

        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)

        fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98)
        return fig, ax
