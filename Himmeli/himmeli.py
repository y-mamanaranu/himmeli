import astropy.units as u
from astropy.units import imperial
import warnings
import matplotlib.pyplot as plt
from pathlib import Path

a4_landscape = [297, 210] * u.mm
a4_portrait = [210, 297] * u.mm


class Himmeli(object):
    def __init__(self, folder=Path(".")):
        self.folder = folder

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
        self, papersize=a4_landscape, horizontalmargin=0.02, verticalmargin=0.02
    ):
        figsize = w0, h0 = papersize.to(imperial.inch).value
        fig, ax = plt.subplots(figsize=figsize)

        self.plot_expansion(ax)

        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()

        dx = x1 - x0
        dy = y1 - y0

        w = 297 * (1 - 2 * horizontalmargin)
        h = 210 * (1 - 2 * verticalmargin)

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
