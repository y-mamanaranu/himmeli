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


def plot_side(ax, x0, y0, x1, y1, xoff=None, yoff=None, linestyle=None):
    x = np.array([x0, x1])
    y = np.array([y0, y1])

    ax.plot(x, y, c="C0", linestyle=linestyle)


def plot_shape(ax, *xy, linestyle=None):
    x, y = zip(*xy)
    ax.plot(x, y, c="C0", linestyle=linestyle)


def plot_isosceles(ax, x0, y0, w, h, transpose=False):
    x1 = x0 + w
    x2 = x0 + w / 2
    y1 = y0
    y2 = y0 + h
    if transpose:
        x1 = x0
        x2 = x0 - h
        y1 = y0 + w
        y2 = y0 + w / 2

    ax.plot([x0, x1, x2], [y0, y1, y2], c="C0")
    ax.plot([x2, x0], [y2, y0], c="C0", linestyle=":")


def plot_trapezoid(ax, x0, y0, w1, w2, h):
    x1 = x0 + w1
    x2 = x0 + (w1 + w2) / 2
    x3 = x0 + (w1 - w2) / 2
    y1 = y0
    y2 = y0 + h
    y3 = y0 + h

    ax.plot([x0, x1, x2, x3], [y0, y1, y2, y3], c="C0")
    ax.plot([x3, x0], [y3, y0], c="C0", linestyle=":")


def plot_circle(ax, x0, y0, r, theta0=0):
    theta = np.linspace(0, 2 * np.pi, 100, endpoint=True)
    x = x0 + r * (np.cos(theta + theta0) - np.cos(theta0))
    y = y0 + r * (np.sin(theta + theta0) - np.sin(theta0))

    ax.plot(x, y, c="C0", linestyle=":")


def plot_polygon(ax, x0, y0, r, n, theta0=0):
    theta = np.linspace(0, 2 * np.pi, n + 1, endpoint=True)
    x = x0 + r * (np.cos(theta + theta0) - np.cos(theta0))
    y = y0 + r * (np.sin(theta + theta0) - np.sin(theta0))

    ax.plot(x, y, c="C0", linestyle=":")
