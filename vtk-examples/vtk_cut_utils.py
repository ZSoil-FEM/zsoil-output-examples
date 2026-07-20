"""
Shared helpers for scripts that cut a 3D ZSoil VTK export (from
zsoil_to_vtk.py) with a plane and plot the intersection in 2D, e.g.
3Ddeepex/plot_geology_profiles.py and 3Ddeepex/plot_result_profiles.py.
"""
from collections import namedtuple

import numpy as np


# A cut plane: `origin`/`normal` define the plane itself (vtkPlane
# convention); `axes` are the two in-plane basis vectors used to project the
# cut's 3D points onto a 2D (x1, x2) drawing plane. `axes` must be orthogonal
# to `normal` for the projection to make sense, but that isn't checked here.
# There's no sensible default for these -- origin/normal/axes are specific
# to each model's geometry and have to be picked by hand every time; see the
# "config section" in each script's __main__.
CutPlane = namedtuple("CutPlane", ["title", "origin", "normal", "axes"])


def project_on_plane(axes, origin, pt):
    """Project a 3D point onto a cut plane's 2D (x1, x2) axes, relative to origin."""
    rel = [pt[k] - origin[k] for k in range(3)]
    return (np.dot(axes[0], rel), np.dot(axes[1], rel))


def GetDiscreteColormap(minmax, colormap='ZSoil_maps', ncol=20):
    """Build a discrete (ncol, evenly spaced) colormap/norm/ticks over
    [minmax[0], minmax[1]]. 'ZSoil_maps' interpolates ZSoil's own blue-to-red
    scale for continuous fields (stress, displacement, ...); 'mat' cycles
    through a fixed 20-color palette meant for discrete material IDs."""
    import matplotlib.colors as colors

    ticks = np.linspace(minmax[0], minmax[1], ncol + 1)

    if colormap == 'ZSoil_maps':
        zsoil = np.array([[41,28,166],
                          [75,11,244],
                          [60,138,255],
                          [61,167,254],
                          [63,190,252],
                          [69,215,245],
                          [83,232,232],
                          [95,220,194],
                          [88,226,143],
                          [81,238,77],
                          [143,251,64],
                          [187,251,117],
                          [216,254,99],
                          [255,255,0],
                          [241,231,35],
                          [239,216,80],
                          [238,186,77],
                          [242,139,64],
                          [254,71,67],
                          [233,6,1],
                          [193,80,4],
                          [167,1,34]])
        colmap = np.array([[v[0]/255,v[1]/255,v[2]/255] for v in zsoil])

        new_positions = np.linspace(0, 1, ncol)
        interp_positions = np.linspace(0, 1, len(zsoil))

        new_colors = np.zeros((ncol, 3))
        for i in range(3):  # RGBA channels
            new_colors[:, i] = np.interp(new_positions, interp_positions, colmap[:, i])
        colmap = new_colors

    elif colormap == 'mat':
        c = ['#ffdcd2',
             '#ffa4a4',
             '#f98568',
             '#da180e',
             '#ffffc6',
             '#def538',
             '#b0b000',
             '#878e2b',
             '#dbfdc6',
             '#8bf391',
             '#5ac960',
             '#658750',
             '#e0e4fe',
             '#bb9af1',
             '#548bcf',
             '#fdcbfe',
             '#e75ae3',
             '#ad5ab4',
             '#abe3e7',
             '#67b1ae']
        colmap = np.array([colors.hex2color(c[k % 20]) for k in range(ncol)])

    cmap = colors.ListedColormap(colmap)
    norm = colors.BoundaryNorm(ticks, len(ticks))

    return cmap, norm, ticks
