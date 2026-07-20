#-------------------------------------------------------------------------------
# Name:        plot_result_profiles
# Purpose:     plot 2D cross sections through a 3D ZSoil model, colored by a
#              nodal or element (volumic) result field, using the same cut
#              planes as plot_geology_profiles.py
#
# Two coloring modes:
#   'nodal'    filled contour bands, from point data (matplotlib
#              tri.Triangulation + tricontourf). Each cut cell is
#              fan-triangulated first, since vtkCutter can produce
#              triangles, quads or larger polygons depending on which 3D
#              cell/direction it cuts through.
#   'element'  flat-shaded, one color per cut cell, from cell data
#              (matplotlib PatchCollection) -- same technique
#              plot_geology_profiles.py uses for the Material array.
#
# Requires:    a VOLUMICS .vtu already exported for the requested step by
#              zsoil_to_vtk.py, carrying both the nodal array ('nodal' mode)
#              and/or the element array ('element' mode) you want to plot.
#              Component layouts (e.g. STRESSES -> XX,YY,XY,ZZ) aren't
#              stored in the .vtu, so this script resolves component names
#              against the ZSoil project's own .rcf (via ZSoilPy3).
#
# Author:      Matthias Preisig
# Created:     2026
# Copyright:   (c) GeoDev Sarl
#-------------------------------------------------------------------------------
# -*- coding: cp1252 -*-
import os
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
from zsoilpy_env import ensure_zsoilpy_on_path
ensure_zsoilpy_on_path()
from zsoil_to_vtk import get_tstr
from vtk_cut_utils import CutPlane, project_on_plane, GetDiscreteColormap

import numpy as np
import vtk
from vtk.util import numpy_support
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from C_Rcf_info import RCF_info


# GetDiscreteColormap's 'ZSoil_maps' scale interpolates smoothly, so unlike
# N_MAT_COLORS in plot_geology_profiles.py, the exact band count here is a
# real, freely adjustable choice, not tied to a fixed palette length.
N_COLORS = 20


def resolve_component(rcf_items, field, comp):
    """
    Resolve `comp` (an int index, or a component-name string like 'XX') to a
    0-based integer index into `field`'s array, validated against
    `rcf_items` -- a list of [name, n_comp, [comp_labels]] as returned by
    RCF_info.give_ele_rcf_items_for_group() (element results) or held in
    rcf.nod_rcf (nodal results).
    """
    for name, n_comp, labels in rcf_items:
        if name.strip() != field:
            continue
        if isinstance(comp, str):
            labels = [label.strip() for label in labels]
            if comp not in labels:
                raise ValueError("'{}' has components {}, not '{}'".format(field, labels, comp))
            return labels.index(comp)
        if not (0 <= comp < n_comp):
            raise ValueError("'{}' has {} component(s) (index 0..{}), comp={} is out of range".format(
                field, n_comp, n_comp - 1, comp))
        return comp
    raise ValueError("Result array '{}' not found (available: {})".format(
        field, ", ".join(item[0].strip() for item in rcf_items)))


def get_component(vtk_array, comp_index):
    """Extract one component from a vtkDataArray as a 1D numpy array."""
    arr = numpy_support.vtk_to_numpy(vtk_array)
    return arr if arr.ndim == 1 else arr[:, comp_index]


def fan_triangles(n_pts):
    """Local (0, i, i+1) index triples fan-triangulating a convex polygon
    with n_pts vertices -- vtkCutter output isn't necessarily triangles."""
    return [(0, i, i + 1) for i in range(1, n_pts - 1)]


def main(cut_planes, prob="3Ddeepex", step=6, zs_path="zsoil-files", vtk_dir="pv",
         mode="nodal", field="DISP_TRA", comp=0, group="VOLUMICS", matlist=()):
    """
    Plot one cross section per entry in `cut_planes` (see
    plot_geology_profiles.py for the CutPlane format), colored by a result
    field read from step `step`'s exported .vtu.

    mode:     'nodal' or 'element' -- see module docstring.
    field:    the result array's name (e.g. 'DISP_TRA', 'PPRESS' for nodal;
              'STRESSES', 'STR_LEVEL', ... for element/volumic), matching
              whatever zsoil_to_vtk.py exported for `group`.
    comp:     which component of `field` to plot -- an int index, or a
              component-name string (e.g. 'XX'), resolved against the
              project's own .rcf.
    matlist:  restrict cut cells to these material IDs; empty means every
              material present.

    All cut planes share one color scale, taken from `field`'s min/max over
    the whole (uncut) model, so the profiles stay comparable to each other
    instead of each auto-scaling to its own local range.
    """
    if mode not in ("nodal", "element"):
        raise ValueError("mode must be 'nodal' or 'element', got {!r}".format(mode))

    tstr = get_tstr(step)
    group_suffix = {"VOLUMICS": "vol"}.get(group, group.lower())
    vtu_path = os.path.join(vtk_dir, "{}_{}_{}.vtu".format(prob, tstr, group_suffix))
    if not os.path.exists(vtu_path):
        sys.exit(
            "Could not find '{}'. Run zsoil_to_vtk.py first to export step {}, "
            "group {}.".format(vtu_path, step, group)
        )

    rcf = RCF_info(os.path.join(zs_path, prob))
    rcf_items = rcf.give_ele_rcf_items_for_group(group) if mode == "element" else rcf.nod_rcf
    comp_index = resolve_component(rcf_items, field, comp)

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_path)
    reader.Update()
    grid = reader.GetOutput()

    data_getter = (lambda g: g.GetPointData()) if mode == "nodal" else (lambda g: g.GetCellData())
    field_array = data_getter(grid).GetArray(field)
    if field_array is None:
        raise RuntimeError("No '{}' {} array in {}".format(
            field, "point" if mode == "nodal" else "cell", vtu_path))
    all_values = get_component(field_array, comp_index)
    vminmax = [float(np.min(all_values)), float(np.max(all_values))]
    cmap, norm, ticks = GetDiscreteColormap(vminmax, colormap='ZSoil_maps', ncol=N_COLORS)

    comp_label = comp if isinstance(comp, str) else str(comp)

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

        if mode == "nodal":
            values = get_component(cut.GetPointData().GetArray(field), comp_index)
            pts2d = np.array([project_on_plane(plane_def.axes, plane_def.origin,
                                                np.array(cut.GetPoint(pid)))
                               for pid in range(cut.GetNumberOfPoints())])

            triangles = []
            for ci in range(cut.GetNumberOfCells()):
                if matlist and material is not None and material.GetTuple1(ci) not in matlist:
                    continue
                cell = cut.GetCell(ci)
                ids = [cell.GetPointIds().GetId(kk) for kk in range(cell.GetNumberOfPoints())]
                triangles.extend([ids[a], ids[b], ids[c]] for a, b, c in fan_triangles(len(ids)))

            if not triangles:
                plt.close(fig)
                continue

            used = np.unique(np.array(triangles))
            xlim = [pts2d[used, 0].min(), pts2d[used, 0].max()]
            ylim = [pts2d[used, 1].min(), pts2d[used, 1].max()]

            triang = tri.Triangulation(pts2d[:, 0], pts2d[:, 1], triangles)
            # tricontourf linearly interpolates the per-vertex values across
            # each triangle and fills between `ticks`' band boundaries, the
            # same bands GetDiscreteColormap built cmap/norm from.
            mappable = ax.tricontourf(triang, values, levels=ticks, cmap=cmap, norm=norm)

        else:
            values = get_component(cut.GetCellData().GetArray(field), comp_index)

            patches = []
            cvect = []
            xlim = [1e10, -1e10]
            ylim = [1e10, -1e10]
            for ci in range(cut.GetNumberOfCells()):
                if matlist and material is not None and material.GetTuple1(ci) not in matlist:
                    continue
                cell = cut.GetCell(ci)
                pts3d = [np.array(cut.GetPoint(cell.GetPointIds().GetId(kk)))
                         for kk in range(cell.GetNumberOfPoints())]
                pts2d = np.array([project_on_plane(plane_def.axes, plane_def.origin, pt) for pt in pts3d])
                for x, y in pts2d:
                    xlim[0] = min(xlim[0], x)
                    xlim[1] = max(xlim[1], x)
                    ylim[0] = min(ylim[0], y)
                    ylim[1] = max(ylim[1], y)
                patches.append(Polygon(pts2d))
                cvect.append(values[ci])

            if not patches:
                plt.close(fig)
                continue

            pc = PatchCollection(patches, edgecolors='none', antialiased=False)
            pc.set_array(np.array(cvect))
            pc.set_cmap(cmap)
            pc.set_norm(norm)
            ax.add_collection(pc)
            mappable = pc

        cb = fig.colorbar(mappable, ax=ax)
        cb.set_label('{} [{}]'.format(field, comp_label))

        ax.set_aspect('equal')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_title('Profile %s' % plane_def.title)

        fig.tight_layout()
        fig.savefig('{}_{}-{}_{}_{}'.format(prob, field, comp_label, mode, plane_def.title))
        plt.close(fig)


if __name__ == "__main__":
    # ------------------------- Config section -------------------------
    # Same cut planes as plot_geology_profiles.py -- see that script for how
    # to adapt origin/normal/axes to your own model.
    cut_planes = [
        CutPlane(title="East-West",
                 origin=(75, 0, -30), normal=(0, 0, 1),
                 axes=((1, 0, 0), (0, 1, 0))),
        CutPlane(title="North-South",
                 origin=(75, 0, -30), normal=(1, 0, 0),
                 axes=((0, 0, 1), (0, 1, 0))),
    ]
    # ----------------------- End of config section ----------------------

    main(cut_planes, mode="nodal", field="DISP_TRA", comp="X")
    main(cut_planes, mode="element", field="STRESSES", comp="XX")
