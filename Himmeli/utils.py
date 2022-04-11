import numpy as np
import matplotlib.pyplot as plt


def radius(theta, offset=0, z=0):
    if hasattr(theta, "__iter__"):
        return np.array([[np.cos(t + offset), np.sin(t + offset), z]
                        for t in theta])
    else:
        return np.array([np.cos(theta + offset), np.sin(theta + offset), z])


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


def max_or_min(arr: np.array,
               coeff: float) -> float:
    """端の点を求める.

    Parameters
    ----------
    arr : np.array
        [description]
    coeff : float
        [description]

    Returns
    -------
    float
        [description]
    """
    if coeff > 0:
        return np.max(arr)
    else:
        return np.min(arr)


def get_carr(p0: np.array,
             p1: np.array,
             p2: np.array) -> np.array:
    """外接円を求める.

    Parameters
    ----------
    p0 : np.array
        1点目
    p1 : np.array
        2点目
    p2 : np.array
        3点目

    Returns
    -------
    np.array
        外接円
    """
    e_x = p0 - p1
    e_x /= np.linalg.norm(e_x)

    e_y = p2 - p1
    e_y -= np.dot(e_x, e_y) * e_x
    e_y /= np.linalg.norm(e_y)

    parr = np.array([p0, p2])
    earr = np.array([e_x, e_y])
    M = np.array(np.dot(parr - p1, earr.T))
    C = np.dot(np.linalg.inv(M), np.linalg.norm(parr - p1, axis=1)**2 / 2)

    Ooff = p1 + np.dot(C, earr)
    r = np.linalg.norm(Ooff - p1)
    s = np.linspace(0, 2 * np.pi, 100, endpoint=True)
    return Ooff + np.vstack([r * np.cos(s), r * np.sin(s), np.zeros(100)]).T


def map_to_2d(parr: np.array,
              p0: np.array,
              p1: np.array,
              p2: np.array,
              off: list,
              coeff: list):
    e_x = p1 - p0
    e_x /= np.linalg.norm(e_x)
    e_x *= coeff[0]

    e_y = p2 - p0
    e_y -= np.dot(e_x, e_y) * e_x
    e_y /= np.linalg.norm(e_y)
    e_y *= coeff[1]

    xarr = off[0] + np.dot(parr - p0, e_x)
    yarr = off[1] + np.dot(parr - p0, e_y)

    return xarr, yarr


def plot_2d(ax: plt.Axes,
            parr: np.array,
            p0: np.array,
            p1: np.array,
            p2: np.array,
            off: list,
            coeff: list,
            ls: str = None) -> list:
    """2次元にプロットする.

    Parameters
    ----------
    ax : plt.Axes
        [description]
    parr : np.array
        [description]
    p0 : np.array
        [description]
    p1 : np.array
        [description]
    p2 : np.array
        [description]
    off : list
        [description]
    coeff : list
        [description]
    ls : str, optional
        [description], by default None

    Returns
    -------
    list
        [description]
    """
    xarr, yarr = map_to_2d(parr, p0, p1, p2, off, coeff)
    ax.plot(xarr, yarr, c="C0", ls=ls)

    return [max_or_min(xarr, coeff[0]), max_or_min(yarr, coeff[1])]
