import numpy as np


def plot_side_3d(ax, x0, y0, z0, x1, y1, z1):
    x = np.array([x0, x1])
    y = np.array([y0, y1])
    z = np.array([z0, z1])

    ax.plot(x, y, z, c="C0")


def plot_surface_3d(ax, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3):
    x = np.array([[x0, x3], [x1, x2]])
    y = np.array([[y0, y3], [y1, y2]])
    z = np.array([[z0, z3], [z1, z2]])

    ax.plot_surface(x, y, z, alpha=0.2, color="C0")
