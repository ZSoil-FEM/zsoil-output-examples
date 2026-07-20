#-------------------------------------------------------------------------------
# Name:        plot_cuts
# Purpose:     plot 2D geology cross sections through a 3D ZSoil model,
#              colored by material, by cutting a VTK volume mesh with a
#              plane and projecting the intersection onto that plane's own
#              2D axes
#
# Requires:    a VOLUMICS .vtu already exported for the requested step by
#              zsoil_to_vtk.py -- it must carry a 'Material' cell array
#              (on by default). This script reads that .vtu for geometry,
#              and the raw ZSoil project (via ZSoilPy3) only for material
#              labels used in the legend.
#
# Author:      Matthias Preisig
# Created:     2026
# Copyright:   (c) GeoDev Sarl
#-------------------------------------------------------------------------------
# -*- coding: cp1252 -*-
import os
import sys
from collections import namedtuple

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
from zsoilpy_env import ensure_zsoilpy_on_path
ensure_zsoilpy_on_path()
from zsoil_to_vtk import get_tstr

import numpy as np
import vtk
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Polygon
from matplotlib.collections import PatchCollection

from C_Mesh import Mesh


# A cut plane: `origin`/`normal` define the plane itself (vtkPlane
# convention); `axes` are the two in-plane basis vectors used to project the
# cut's 3D points onto a 2D (x1, x2) drawing plane. `axes` must be orthogonal
# to `normal` for the projection to make sense, but that isn't checked here.
# There's no sensible default for these -- origin/normal/axes are specific
# to each model's geometry and have to be picked by hand every time; see the
# "config section" in __main__ below.
CutPlane = namedtuple("CutPlane", ["title", "origin", "normal", "axes"])

# The 'mat' colormap below is a fixed 20-color palette (not a parameter --
# picking a different count wouldn't add colors, just repeat this same 20).
N_MAT_COLORS = 20


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


def project_on_plane(axes, origin, pt):
    """Project a 3D point onto a cut plane's 2D (x1, x2) axes, relative to origin."""
    rel = [pt[k] - origin[k] for k in range(3)]
    return (np.dot(axes[0], rel), np.dot(axes[1], rel))


def main(cut_planes, prob="3Ddeepex", step=6, zs_path="zsoil-files", vtk_dir="pv",
         matlist=(1, 2, 3, 4, 5, 6)):
    """
    Plot one geology cross section per entry in `cut_planes` (a list of
    CutPlane; there's no default -- see the module docstring for CutPlane),
    colored by material (restricted to `matlist`; pass an empty list/tuple
    for every material present). Reads step `step`'s VOLUMICS export from
    `vtk_dir` (see module docstring) and material labels from the ZSoil
    project at `zs_path`/`prob`.
    """
    tstr = get_tstr(step)
    vtu_path = os.path.join(vtk_dir, "{}_{}_vol.vtu".format(prob, tstr))
    if not os.path.exists(vtu_path):
        sys.exit(
            "Could not find '{}'. Run zsoil_to_vtk.py first to export step {} "
            "(at least the VOLUMICS group).".format(vtu_path, step)
        )

    mesh = Mesh(os.path.join(zs_path, prob))
    inp_material_labels = [mat.index_in_inp for mat in mesh.materials]
    matlabels = [mesh.materials[inp_material_labels.index(k)].label for k in matlist]

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_path)
    reader.Update()
    grid = reader.GetOutput()

    for plane_def in cut_planes:
        fig = plt.figure(figsize=(12, 4.4))
        fig.text(0.01, 0.01, prob, size=8)
        ax = fig.add_subplot(111)

        plane = vtk.vtkPlane()
        plane.SetOrigin(plane_def.origin)
        plane.SetNormal(plane_def.normal)
        cutter = vtk.vtkCutter()
        cutter.SetInputData(grid)
        cutter.SetCutFunction(plane)
        cutter.Update()
        cut = cutter.GetOutput()

        material = cut.GetCellData().GetArray('Material')
        if material is None:
            raise RuntimeError(
                "No 'Material' cell array in {} -- re-export it with zsoil_to_vtk.py "
                "(Material is included by default for every group).".format(vtu_path))

        patches = []
        cvect = []
        minmax = [[1e10, -1e10], [1e10, -1e10]]
        for ke in range(cut.GetNumberOfCells()):
            mat_id = material.GetTuple1(ke)
            if matlist and mat_id not in matlist:
                continue
            cell = cut.GetCell(ke)
            pts3d = [np.array(cut.GetPoint(cell.GetPointIds().GetId(kk)))
                     for kk in range(cell.GetNumberOfPoints())]
            pts2d = np.array([project_on_plane(plane_def.axes, plane_def.origin, pt) for pt in pts3d])
            for x, y in pts2d:
                minmax[0][0] = min(minmax[0][0], x)
                minmax[0][1] = max(minmax[0][1], x)
                minmax[1][0] = min(minmax[1][0], y)
                minmax[1][1] = max(minmax[1][1], y)
            patches.append(Polygon(pts2d))
            cvect.append(mat_id)

        cmap, norm, ticks = GetDiscreteColormap([1, N_MAT_COLORS], colormap='mat', ncol=N_MAT_COLORS)

        pc = PatchCollection(patches, edgecolors='none', antialiased=False)
        pc.set_array(np.array(cvect))
        pc.set_cmap(cmap)
        pc.set_norm(norm)
        ax.add_collection(pc)

        ax.set_aspect('equal')
        ax.set_xlim(minmax[0])
        ax.set_ylim(minmax[1])
        ax.set_title('Profile %s' % plane_def.title)

        handles = [Patch(color=cmap((km - 1) / N_MAT_COLORS), ec='k', label=matlabels[kk])
                   for kk, km in enumerate(matlist)]
        ax.legend(loc='lower left', handles=handles)

        fig.tight_layout()
        fig.savefig('{}_geol_{}'.format(prob, plane_def.title))
        plt.close(fig)


if __name__ == "__main__":
    # ------------------------- Config section -------------------------
    # Cut planes for THIS model -- origin/normal/axes have to be picked by
    # hand for your own geometry (e.g. from a plan view of the model).
    cut_planes = [
        CutPlane(title="East-West",
                 origin=(75, 0, -30), normal=(0, 0, 1),
                 axes=((1, 0, 0), (0, 1, 0))),
        CutPlane(title="North-South",
                 origin=(75, 0, -30), normal=(1, 0, 0),
                 axes=((0, 0, 1), (0, 1, 0))),
    ]
    # ----------------------- End of config section ----------------------

    main(cut_planes)
