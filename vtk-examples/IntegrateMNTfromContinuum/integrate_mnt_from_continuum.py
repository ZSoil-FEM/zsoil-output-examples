#-------------------------------------------------------------------------------
# Name:        integrate_mnt_from_continuum
# Purpose:     integrate a volumic tunnel lining's 3D stress field across its
#              cross section, at points around a ring, into equivalent
#              beam-like stress resultants -- bending moment (M), normal
#              force (N) and transverse shear (T) -- the same physical
#              quantities plot_wall_diagrams.py reads directly from SMOMENT/
#              SQFORCE, derived here instead from a solid continuum lining,
#              which has no shell-style force resultants of its own.
#
# Geometry: a thin reference SHELLS surface ("neutral axis", no stiffness of
# its own -- export it with zsoil_to_vtk.py's include_inactive=True) runs
# along the lining's mid-thickness. Cutting it with one of `axial_sections`
# gives points around the ring at that axial location; cutting the VOLUMICS
# lining a second time, perpendicular to the ring at each of those points,
# gives the thin strip of solid elements spanning the lining's thickness
# there. Their STRESSES (global XX/YY/XY/ZZ/XZ/YZ) are rotated into a local
# (radial, tangential/hoop, axial) frame and integrated across the strip:
#   integral(sigma_hoop dA)                -> N
#   integral(sigma_hoop * radial_offset dA) -> M
#   integral(tau_radial-hoop dA)            -> T
#
# Requires:    a VOLUMICS .vtu (with Material, STRESSES) and a SHELLS .vtu
#              (with Material; include_inactive=True, since the neutral
#              axis carries no stiffness) already exported for the
#              requested step by zsoil_to_vtk.py.
#
# Author:      Matthias Preisig
# Created:     2026 (originally written 2018/10/24, rewritten)
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
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


def _stress_tensor(stress6):
    """Build the symmetric 3x3 stress tensor from ZSoil's 6-component
    VOLUMICS STRESSES order (XX, YY, XY, ZZ, XZ, YZ)."""
    xx, yy, xy, zz, xz, yz = stress6
    return np.array([[xx, xy, xz],
                      [xy, yy, yz],
                      [xz, yz, zz]])


def ring_centroid(axis_line, neutral_axis_mat):
    """
    Geometric center of the neutral-axis ring at this axial section (mean
    of its segment midpoints). Used as the reference point for
    integrate_ring_section()'s `radial`-outward convention -- the axial
    section's own origin (e.g. `axial_sections[k].origin`) is only a point
    ON the cutting plane, not necessarily inside the ring (e.g. it may sit
    below the tunnel to make the plane easy to define), so it can't be used
    for that directly.
    """
    material = axis_line.GetCellData().GetArray("Material")
    pts = []
    for ka in range(axis_line.GetNumberOfCells()):
        if material.GetTuple1(ka) != neutral_axis_mat:
            continue
        line = axis_line.GetCell(ka)
        pts.append(np.array(line.GetPoints().GetPoint(0)))
        pts.append(np.array(line.GetPoints().GetPoint(1)))
    if not pts:
        return None
    return np.mean(pts, axis=0)


def integrate_ring_section(axial_cut, pt0, pt1, axial_normal, axis_center, lining_mat, select_tol):
    """
    Integrate the lining's stress field across its thickness at one point
    around the ring -- the midpoint of a neutral-axis segment (pt0 -> pt1).

    `axial_cut` is the VOLUMICS grid already cut by the ring's axial plane
    (so its cells are 2D polygons lying in that plane). Cutting it again,
    through the segment midpoint and perpendicular to the ring (i.e. a
    plane containing `radial` and `axial_normal`), leaves the thin strip of
    lining cross section at this one location.

    `axis_center` (see ring_centroid() -- must be inside the ring, not just
    on the cutting plane) fixes `radial`'s sign to consistently point
    outward: pt0/pt1 come from independent line cells that vtkCutter orders
    arbitrarily, so cross(tangent, axial_normal) alone flips sign
    unpredictably from one ring segment to the next, which would turn the
    M/eccentricity diagrams into a zigzag instead of a smooth curve offset
    from the ring.

    Returns None if no `lining_mat` cells are found within `select_tol` of
    the segment midpoint (e.g. this ring point falls outside the lining's
    actual extent); otherwise a dict with the integrated M/N/T, the local
    frame, and (p0, p1, sigma_hoop) per strip cell for the caller to render
    a stress-distribution patch.
    """
    tangent = pt1 - pt0
    tangent = tangent / np.linalg.norm(tangent)
    radial = np.cross(tangent, axial_normal)
    radial = radial / np.linalg.norm(radial)
    origin = 0.5 * (pt0 + pt1)
    if np.dot(radial, origin - axis_center) < 0:
        radial = -radial
    loc_syst = np.array([radial, tangent, axial_normal])

    plane = vtk.vtkPlane()
    plane.SetOrigin(origin)
    plane.SetNormal(tangent)
    cutter = vtk.vtkCutter()
    cutter.SetInputData(axial_cut)
    cutter.SetCutFunction(plane)
    cutter.Update()
    strip = cutter.GetOutput()

    material = strip.GetCellData().GetArray("Material")
    stresses = strip.GetCellData().GetArray("STRESSES")

    M = T = N = 0.0
    cells = []
    for ci in range(strip.GetNumberOfCells()):
        if material.GetTuple1(ci) != lining_mat:
            continue
        cell = strip.GetCell(ci)
        p0 = np.array(strip.GetPoint(cell.GetPointId(0)))
        p1 = np.array(strip.GetPoint(cell.GetPointId(1)))
        ptm = 0.5 * (p0 + p1)
        r = np.dot(ptm - pt0, radial)
        if abs(r) > select_tol:
            continue

        stress_local = loc_syst.dot(_stress_tensor(stresses.GetTuple(ci))).dot(loc_syst.T)
        length = cell.GetLength2() ** 0.5
        sigma_hoop = stress_local[1, 1]

        M += length * r * sigma_hoop
        T += length * stress_local[0, 1]
        N += length * sigma_hoop
        cells.append((p0, p1, sigma_hoop))

    if not cells:
        return None
    return {"M": M, "T": T, "N": N, "dl": np.linalg.norm(pt1 - pt0),
            "origin": origin, "tangent": tangent, "radial": radial, "cells": cells}


def main(axial_sections, prob, step, step_labels, vtk_dir,
         lining_mat, neutral_axis_mat, up_dir,
         ring_select_tol, stress_diagram_scale, moment_scale, normal_scale):
    """
    For each entry in `axial_sections` (a CutPlane through the tunnel axis;
    `axes` gives the 2D drawing plane, `axes[0]` = up_dir x normal), plot a
    ring diagram: the lining cross section, a colored stress diagram showing
    the local hoop stress distribution across the lining thickness, and
    bending-moment/normal-force/eccentricity diagrams offset from the
    neutral axis.

    lining_mat:          VOLUMICS material of the tunnel lining being
                          integrated.
    neutral_axis_mat:    SHELLS material of the reference neutral-axis
                          surface (its geometry only; no stiffness/results
                          are read from it).
    ring_select_tol:      how close (in radial distance) a lining strip cell
                          must be to the neutral axis to be integrated
                          there. This makes sure that stresses from the opposite tunnel
                          wall are excluded in the integration procedure.
    stress_diagram_scale: offsets the stress-diagram patches by
                          stress_diagram_scale * sigma_hoop.
    moment_scale/normal_scale: offset the M/N diagrams from the neutral
                          axis by <scale> * M or N.
    """
    tstr = get_tstr(step)

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(os.path.join(vtk_dir, "{}_{}_vol.vtu".format(prob, tstr)))
    reader.Update()
    volumics = reader.GetOutput()

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(os.path.join(vtk_dir, "{}_{}_shell.vtu".format(prob, tstr)))
    reader.Update()
    neutral_axis_surface = reader.GetOutput()

    for ks, section in enumerate(axial_sections):
        origin0 = np.array(section.origin)
        normal0 = np.array(section.normal)
        normal0 = normal0 / np.linalg.norm(normal0)

        plane0 = vtk.vtkPlane()
        plane0.SetOrigin(origin0)
        plane0.SetNormal(normal0)

        cutter = vtk.vtkCutter()
        cutter.SetInputData(volumics)
        cutter.SetCutFunction(plane0)
        cutter.Update()
        axial_cut = cutter.GetOutput()

        cutter = vtk.vtkCutter()
        cutter.SetInputData(neutral_axis_surface)
        cutter.SetCutFunction(plane0)
        cutter.Update()
        axis_line = cutter.GetOutput()
        axis_material = axis_line.GetCellData().GetArray("Material")

        ring_center = ring_centroid(axis_line, neutral_axis_mat)
        if ring_center is None:
            print("integrate_mnt_from_continuum: no material-{} neutral axis found for '{}'".format(
                neutral_axis_mat, section.title))
            continue

        fig = plt.figure(figsize=(8, 10))
        fig.text(0.01, 0.01, prob, size=8)
        ax = fig.add_subplot(111)

        # background: the lining cross section itself, filled and semi-transparent
        xlim = [1e10, -1e10]
        ylim = [1e10, -1e10]

        def _track(pts2d):
            for x, y in pts2d:
                xlim[0] = min(xlim[0], x)
                xlim[1] = max(xlim[1], x)
                ylim[0] = min(ylim[0], y)
                ylim[1] = max(ylim[1], y)

        lining_material = axial_cut.GetCellData().GetArray("Material")
        patches = []
        cvect = []
        for ci in range(axial_cut.GetNumberOfCells()):
            if lining_material.GetTuple1(ci) != lining_mat:
                continue
            cell = axial_cut.GetCell(ci)
            pts2d = [project_on_plane(section.axes, origin0, cell.GetPoints().GetPoint(kp))
                     for kp in range(cell.GetNumberOfPoints())]
            patches.append(Polygon(pts2d))
            cvect.append(lining_material.GetTuple1(ci))
            _track(pts2d)
        pc = PatchCollection(patches, edgecolors='none', alpha=0.24, cmap=plt.get_cmap('Set1'))
        pc.set_array(np.array(cvect))
        ax.add_collection(pc)

        # walk the ring: one integration per neutral-axis segment
        MZ, QY, NX, PT, NORMALS_2D = [], [], [], [], []
        diagram_patches, diagram_colors = [], []
        for ka in range(axis_line.GetNumberOfCells()):
            if axis_material.GetTuple1(ka) != neutral_axis_mat:
                continue
            line = axis_line.GetCell(ka)
            pt0 = np.array(line.GetPoints().GetPoint(0))
            pt1 = np.array(line.GetPoints().GetPoint(1))

            result = integrate_ring_section(axial_cut, pt0, pt1, normal0, ring_center, lining_mat, ring_select_tol)
            if result is None:
                continue

            dl = result["dl"]
            M, T, N = result["M"], result["T"], result["N"]
            tangent = result["tangent"]
            radial = result["radial"]

            for p0, p1, sigma_hoop in result["cells"]:
                p0_off = p0 + tangent * stress_diagram_scale * sigma_hoop
                p1_off = p1 + tangent * stress_diagram_scale * sigma_hoop
                quad = [project_on_plane(section.axes, origin0, p) for p in (p0, p0_off, p1_off, p1)]
                diagram_patches.append(Polygon(quad))
                diagram_colors.append('r' if sigma_hoop > 0 else 'b')

            # M/T come from one integration point at the segment's center;
            # plot M varying linearly across the segment assuming a
            # constant T, same convention as BerlinSand's plot_excavation2D.py.
            MZ.append([M - 0.5 * dl * T, M + 0.5 * dl * T])
            QY.append(T)
            NX.append(N)
            PT.append([project_on_plane(section.axes, origin0, pt0),
                       project_on_plane(section.axes, origin0, pt1)])
            # M/N/eccentricity diagrams are offset perpendicular to the ring
            # (radially), not along it -- project the radial direction, not
            # the tangent.
            NORMALS_2D.append(project_on_plane(section.axes, [0, 0, 0], radial))

        if not MZ:
            print("integrate_mnt_from_continuum: no '{}' lining ring found for '{}'".format(
                lining_mat, section.title))
            plt.close(fig)
            continue

        pc = PatchCollection(diagram_patches, edgecolors='none', facecolors=diagram_colors,
                             antialiased=True, alpha=0.4)
        ax.add_collection(pc)

        # eccentricity e = M/N of the resultant axial force, offset from the
        # neutral axis; capped when N is too small for e to be meaningful.
        ecc = []
        for (m0, m1), n in zip(MZ, NX):
            if abs((m0 + m1) / (2 * n)) > 2:
                ecc.append([0, 0])
            else:
                ecc.append([m0 / n, m1 / n])

        for kp in range(len(PT)):
            m_line = [[PT[kp][kk][0] + NORMALS_2D[kp][0] * moment_scale * MZ[kp][kk] for kk in range(2)],
                      [PT[kp][kk][1] + NORMALS_2D[kp][1] * moment_scale * MZ[kp][kk] for kk in range(2)]]
            n_line = [[PT[kp][kk][0] - NORMALS_2D[kp][0] * normal_scale * NX[kp] for kk in range(2)],
                      [PT[kp][kk][1] - NORMALS_2D[kp][1] * normal_scale * NX[kp] for kk in range(2)]]
            e_line = [[PT[kp][kk][0] + NORMALS_2D[kp][0] * ecc[kp][kk] for kk in range(2)],
                      [PT[kp][kk][1] + NORMALS_2D[kp][1] * ecc[kp][kk] for kk in range(2)]]
            ax.plot(*m_line, color='b', lw=2)
            ax.plot(*n_line, color='g', lw=2)
            ax.plot([PT[kp][kk][0] for kk in range(2)],
                    [PT[kp][kk][1] for kk in range(2)], 'k--')
            ax.plot(*e_line, color='r')
            _track(zip(m_line[0], m_line[1]))
            _track(zip(n_line[0], n_line[1]))

        m_flat = [v for pair in MZ for v in pair]
        fig.text(0.95, 0.55, 'Step T=%g: %s' % (step, step_labels[step]),
                 ha='right', va='top', size=14)
        astring = (u'Bending moment (blue)\n$M_{min}=$%1.1f kNm/m, $M_{max}=$%1.1f kNm/m'
                   % (min(m_flat), max(m_flat)))
        astring += (u'\n\nNormal force (green)\n$N_{min}=$%1.1f kN/m, $N_{max}=$%1.1f kN/m'
                    % (min(NX), max(NX)))
        astring += u'\n\nPositive M for tension on the outside face\nResultant eccentricity in red'
        fig.text(0.95, 0.52, astring, ha='right', va='top', size=10)

        margin = 0.05 * max(xlim[1] - xlim[0], ylim[1] - ylim[0])
        ax.set_xlim(xlim[0] - margin, xlim[1] + margin)
        ax.set_ylim(ylim[0] - margin, ylim[1] + margin)
        ax.axis('off')
        ax.set_aspect('equal')

        fig.tight_layout()
        fig.savefig('{}_stresses_T{:g}_sect{}'.format(prob, step, ks))
        plt.close(fig)


if __name__ == "__main__":
    # ------------------------- Config section -------------------------
    # Two cross sections along the tunnel axis (Z); the drawing plane at
    # each is spanned by (up x normal, up) -- see CutPlane in vtk_cut_utils.py.
    up = (0, 1, 0)
    axial_sections = [
        CutPlane(title="Section 1", origin=(0, 0, 0.2), normal=(0, 0, 1),
                 axes=(np.cross(up, (0, 0, 1)), up)),
        CutPlane(title="Section 2", origin=(0, 0, 3.8), normal=(0, 0, 1),
                 axes=(np.cross(up, (0, 0, 1)), up)),
    ]
    # ----------------------- End of config section ----------------------

    main(axial_sections,
         prob="3Dtunnel",
         step=2,
         step_labels={2: u"End of construction"},
         vtk_dir="pv",
         lining_mat=2,
         neutral_axis_mat=4,
         up_dir=up,
         ring_select_tol=0.5,
         stress_diagram_scale=0.00005,
         moment_scale=0.01,
         normal_scale=0.001)
