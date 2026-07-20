#-------------------------------------------------------------------------------
# Name:        plot_wall_diagrams
# Purpose:     plot horizontal displacement, bending moment and shear force
#              diagrams vs. altitude for a vertical profile through a shell
#              (wall) panel, at a chosen along-wall location -- the same
#              "M/T/disp vs depth" diagram BerlinSand's plot_excavation2D.py
#              draws for a 2D wall, generalized to any wall orientation in a
#              3D model.
#
# The displacement/moment/shear components plotted are defined by the
# WALL's OWN geometry (its LocalX/LocalY/LocalZ, as exported by
# zsoil_to_vtk.py's get_element_local_axes()), not by the cut plane's
# orientation -- so the same code gives physically correct results for
# walls running in any direction, without hardcoding which raw
# SMOMENT/SQFORCE component happens to be "the right one" for a particular
# wall (that would silently break, or silently be wrong, for a wall running
# the other way). moment_and_shear_along() reproduces ZSoilPy3's own
# Element.set_transf_matrices() + Shell_EleResults transform in pure numpy,
# from the already-exported local axes and raw SMOMENT/SQFORCE -- verified
# to match the SDK-based computation exactly (0.0 difference) across every
# segment used here.
#
# This script only reads the exported .vtu -- it does not need the ZSoilPy3
# SDK (or a ZSoil license) at all, unlike the first version of this script.
#
# Requires:    a SHELLS .vtu already exported for the requested step by
#              zsoil_to_vtk.py, carrying (all on by default): Material,
#              ElementIndex, LocalX/LocalY/LocalZ, SMOMENT, SQFORCE and the
#              DISP_TRA point array.
#
# Author:      Matthias Preisig
# Created:     2026 (originally written 2018/06/22, rewritten)
# Copyright:   (c) GeoDev Sarl
#-------------------------------------------------------------------------------
# -*- coding: cp1252 -*-
import os
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
from vtk_cut_utils import CutPlane, project_on_plane, get_tstr

import numpy as np
import vtk
import matplotlib.pyplot as plt


def find_profile_segments(cut, plane_def, select_tol, wall_mat):
    """
    From a plane's vtkCutter output `cut` (of a SHELLS .vtu), return the cut
    segments belonging to material `wall_mat` that lie within `select_tol`
    of the plane's origin along `plane_def.axes[0]` (which must be a unit
    vector for `select_tol` to be a true distance) -- i.e. one
    representative vertical profile through the wall, at the along-wall
    location given by the origin. A cut plane can cross more than one wall
    (e.g. two parallel walls on opposite sides of an excavation); this is
    how a single one is picked out.

    Returns a list of dicts, one per segment, with the raw per-cell/per-point
    data needed to compute displacement/moment/shear in main().
    """
    cdata = cut.GetCellData()
    required = ["Material", "ElementIndex", "LocalX", "LocalY", "LocalZ", "SMOMENT", "SQFORCE"]
    arrays = {}
    for name in required:
        arr = cdata.GetArray(name)
        if arr is None:
            raise RuntimeError(
                "Cut is missing cell array '{}' -- re-export with zsoil_to_vtk.py; SHELLS gets "
                "all of {} by default.".format(name, required)
            )
        arrays[name] = arr
    disp_arr = cut.GetPointData().GetArray("DISP_TRA")
    if disp_arr is None:
        raise RuntimeError("Cut is missing the 'DISP_TRA' point array -- re-export with zsoil_to_vtk.py.")

    segments = []
    for ci in range(cut.GetNumberOfCells()):
        if arrays["Material"].GetTuple1(ci) != wall_mat:
            continue
        cell = cut.GetCell(ci)
        id0 = cell.GetPointIds().GetId(0)
        id1 = cell.GetPointIds().GetId(1)
        pt0 = np.array(cut.GetPoint(id0))
        pt1 = np.array(cut.GetPoint(id1))
        along0 = project_on_plane(plane_def.axes, plane_def.origin, pt0)[0]
        if abs(along0) > select_tol:
            continue
        segments.append({
            "y0": pt0[1], "y1": pt1[1],
            "disp0": np.array(disp_arr.GetTuple(id0)),
            "disp1": np.array(disp_arr.GetTuple(id1)),
            "local_x": np.array(arrays["LocalX"].GetTuple(ci)),
            "local_y": np.array(arrays["LocalY"].GetTuple(ci)),
            "local_z": np.array(arrays["LocalZ"].GetTuple(ci)),
            "M": arrays["SMOMENT"].GetTuple(ci),   # (Mxx, Myy, Mxy)
            "Q": arrays["SQFORCE"].GetTuple(ci),   # (Qx, Qy)
        })
    return segments


def moment_and_shear_along(vertical_dir, local_x, local_y, M, Q):
    """
    Resolve one element's raw local-axis SMOMENT (Mxx, Myy, Mxy) / SQFORCE
    (Qx, Qy) into the moment/shear about/across a horizontal axis aligned
    with `vertical_dir` -- the physically meaningful "cantilever bending"
    quantities for a wall, regardless of which way it faces. Reproduces
    ZSoilPy3's own Element.set_transf_matrices() + Shell_EleResults
    transform (Muser = T*M, Quser = T*Q, from C_Element_Surf_generic.py's
    get_T_user_VEC/get_T_user_TNS) directly from the exported local axes,
    without needing the SDK.

    Returns (M, T), or None if `vertical_dir` is ~parallel to this
    element's normal (can't define an in-plane "vertical" direction then).
    """
    e0 = np.dot(vertical_dir, local_x)
    e1 = np.dot(vertical_dir, local_y)
    norm = np.hypot(e0, e1)
    if norm < 1e-6:
        return None
    e0, e1 = e0 / norm, e1 / norm

    Mxx, Myy, Mxy = M
    Qx, Qy = Q
    M_v = Mxx * e0 * e0 + Myy * e1 * e1 + 2.0 * Mxy * e0 * e1
    Q_v = Qx * e0 + Qy * e1
    return M_v, Q_v


def main(cut_planes, prob, step, vtk_dir, wall_mat, vertical_dir, select_tol):
    """
    For each entry in `cut_planes` (a CutPlane; axes[0] picks the along-wall
    location, see find_profile_segments), plot a 3-panel diagram (horizontal
    displacement perpendicular to the wall, bending moment, shear force) vs.
    altitude, for the shell wall (material `wall_mat`) crossing that plane
    within `select_tol` of its origin.

    `vertical_dir` (e.g. (0, 1, 0)) is the axis moment/shear get resolved
    about/across (see moment_and_shear_along) -- the physically meaningful
    "cantilever bending" direction, independent of how the wall happens to
    face in the global XYZ frame.
    """
    tstr = get_tstr(step)
    vtu_path = os.path.join(vtk_dir, "{}_{}_shell.vtu".format(prob, tstr))
    if not os.path.exists(vtu_path):
        sys.exit(
            "Could not find '{}'. Run zsoil_to_vtk.py first to export step {}, "
            "group SHELLS.".format(vtu_path, step)
        )

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_path)
    reader.Update()
    grid = reader.GetOutput()

    for plane_def in cut_planes:
        plane = vtk.vtkPlane()
        plane.SetOrigin(plane_def.origin)
        plane.SetNormal(plane_def.normal)
        cutter = vtk.vtkCutter()
        cutter.SetInputData(grid)
        cutter.SetCutFunction(plane)
        cutter.Update()
        cut = cutter.GetOutput()

        segments = find_profile_segments(cut, plane_def, select_tol, wall_mat)
        if not segments:
            print("plot_wall_diagrams: no material-{} wall segments found near the origin "
                  "for '{}'".format(wall_mat, plane_def.title))
            continue

        fig = plt.figure(figsize=(10, 8))
        fig.text(0.01, 0.01, prob, size=8)
        AX = [fig.add_subplot(1, 3, kk + 1) for kk in range(3)]

        for seg in segments:
            result = moment_and_shear_along(
                vertical_dir, seg["local_x"], seg["local_y"], seg["M"], seg["Q"])
            if result is None:
                print("plot_wall_diagrams: skipping a near-horizontal segment "
                      "(vertical_dir ~ element normal)")
                continue
            M, T = result
            u0 = np.dot(seg["disp0"], seg["local_z"]) * 1e3
            u1 = np.dot(seg["disp1"], seg["local_z"]) * 1e3
            y0, y1 = seg["y0"], seg["y1"]

            # M/T come from one Gauss point at the element's center; plot M
            # varying linearly across the element assuming a constant shear
            # T, same convention as BerlinSand's plot_excavation2D.py.
            dy = 0.5 * (y1 - y0)
            AX[0].plot([u0, u1], [y0, y1], color='C0')
            AX[1].plot([M - dy * T, M + dy * T], [y0, y1], color='C0')
            AX[2].plot([0, T, T, 0], [y0, y0, y1, y1], color='C0')

        AX[0].set_xlabel(u'Horizontal displacement perp. to wall [mm]')
        AX[1].set_xlabel(u'Bending moment [kNm/m]')
        AX[2].set_xlabel(u'Shear force [kN/m]')
        for ax in AX:
            ax.set_ylabel(u'Altitude [m]')
            ax.grid('on')

        fig.text(0.5, 0.97, 'Profile ' + plane_def.title, size=14, ha='center')
        fig.tight_layout()
        fig.savefig('{}_wall_diagrams_{}'.format(prob, plane_def.title))
        plt.close(fig)


if __name__ == "__main__":
    # ------------------------- Config section -------------------------
    # Wall shells are material 21 (checked against the mesh). East-West
    # running walls sit at Z=-50 and Z=-10, North-South running walls at
    # X=50 and X=100. Each cut plane below picks out one wall at one
    # along-wall location via its origin + axes[0] (see
    # find_profile_segments) -- axes[0] must be a unit vector.
    cut_planes = [
        CutPlane(title="North wall, center",
                 origin=(75, 0, -50), normal=(1, 0, 0),
                 axes=((0, 0, 1), (0, 1, 0))),
        CutPlane(title="East wall, center",
                 origin=(100, 0, -30), normal=(0, 0, 1),
                 axes=((1, 0, 0), (0, 1, 0))),
    ]
    # ----------------------- End of config section ----------------------

    main(cut_planes,
         prob="3Ddeepex",
         step=6,
         vtk_dir="pv",
         wall_mat=21,
         vertical_dir=(0, 1, 0),
         select_tol=1)
