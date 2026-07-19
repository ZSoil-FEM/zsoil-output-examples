#-------------------------------------------------------------------------------
# Name:        plot_piles_3D
# Purpose:     visualize a 3D pile group's axial force and shaft plasticity
#              in one combined view (an oblique 2D projection), colored to
#              flag where the shaft interface has reached full shear
#              mobilization - a combined view PostPro does not offer
#
# Author:      Matthias Preisig
# Created:     2026
# Copyright:   (c) GeoDev Sarl
# Licence:     <your licence>
#-------------------------------------------------------------------------------
# -*- coding: cp1252 -*-
import os
import sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from zsoilpy_env import ensure_zsoilpy_on_path
ensure_zsoilpy_on_path()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

import C_Mesh
import C_HistoryOfExecution as C_His
import C_Rcf_info as C_Rcf
import C_BeamResults as C_Beam_Res
import C_Contact_3D_Results as C_Cnt_Res

# ZSoil's interface PLA_CODE convention: 16 marks a Gauss point where the
# shaft has reached full shear (friction) mobilization - see the ZSoil
# interface-element manual for the complete code table.
SHEAR_FAILURE_CODE = 16


def make_projector(azim=0.3, elev=0.5, elev_correction=1.2):
    """
    Build a function that projects 3D points onto a 2D plane at the given
    azimuth/elevation (radians), for a quick oblique view of a pile group.
    `elev_correction` is an empirical adjustment to `elev` that makes the
    resulting view look more natural; it has no physical meaning - tune it
    by eye if the projection looks off for your own geometry.
    """
    elev = elev * elev_correction
    a_vec = np.array([np.cos(azim), np.sin(azim), 0])
    normal = np.cos(elev) * a_vec + np.array([0, 0, np.sin(elev)])
    z_vec = np.array([0, 1, 0])
    y_comp = z_vec - (z_vec @ normal) * normal
    y_comp = y_comp / np.linalg.norm(y_comp)
    x_comp = np.cross(y_comp, normal)
    proj_mat = np.vstack([x_comp, y_comp])

    def project_2d(point0):
        point = point0 * np.array([1, 1, -1])
        return point @ proj_mat.T

    return project_2d


project_2d = make_projector()


def main(pathname, prob, step_to_plot, pile_mat=9, shaft_mat=8, force_scale=5e-4):
    """
    For the solution step at `step_to_plot`, draw every pile as a line
    thickened by a band proportional to its axial force (`force_scale`
    controls the band width), colored red wherever the shaft interface
    (material `shaft_mat`) has reached SHEAR_FAILURE_CODE.
    """
    file_path = os.path.join(pathname, prob)

    his = C_His.HistoryOfExecution(file_path)
    rcf = C_Rcf.RCF_info(file_path)
    mesh = C_Mesh.Mesh(file_path)

    all_converged_solutions = his.give_converged_solutions()

    # select solutions without stability driver and integer time instances:
    solution_indices = []
    for sol_ind in all_converged_solutions:
        time = his.data[sol_ind][his.DATA_TIME]
        SF = his.data[sol_ind][his.DATA_SF]
        if SF == 0 and time == step_to_plot:
            solution_indices.append(sol_ind)

    # get beam and interface results:
    beam_rsl = C_Beam_Res.Beam_EleResults(mesh, his, rcf)
    cnt_rsl = C_Cnt_Res.Contact_3D_Results(mesh, his, rcf)

    # get indices for all nodes:
    all_nodes = mesh.get_list_of_nodes()

    # select pile beams and shaft interface elements:
    beams = mesh.get_list_of_elements('BEAMS', mat_filter=[pile_mat])
    shaftlist = mesh.get_list_of_elements('CONTACT', mat_filter=[shaft_mat])

    for sol in solution_indices:
        time = his.data[sol][his.DATA_TIME]

        pla_code = cnt_rsl.get_rsl_vec_for_sel_elements([ele.index for ele in shaftlist], sol, 'PLA_CODE')
        PC = [[v[0], v[1]] for v in pla_code]
        forces = beam_rsl.get_forces_vec_for_sel_elements([ele.index for ele in beams], sol)
        Fn = [v[0] for v in forces]

        fig = plt.figure(figsize=(12, 10))
        fig.text(0.01, 0.01, prob, size=8)
        ax = fig.add_subplot(111)

        for kke, beam in enumerate(beams):
            p0 = all_nodes[beam.nodes[0] - 1].get_xyz()
            p1 = all_nodes[beam.nodes[1] - 1].get_xyz()
            pF = np.array([Fn[kke] * force_scale, 0, 0])
            p0_2d = project_2d(p0)
            p1_2d = project_2d(p1)
            p0F_2d = project_2d(p0 + pF)
            p1F_2d = project_2d(p1 + pF)
            ax.plot([p0_2d[0], p1_2d[0]], [p0_2d[1], p1_2d[1]], 'k', lw=1)
            failed = min(PC[kke][kgp] for kgp in [0, 1]) == SHEAR_FAILURE_CODE
            ax.add_patch(Polygon([p0_2d, p1_2d, p1F_2d, p0F_2d],
                                  fc='r' if failed else 'k', ec=None, alpha=0.6))

        ax.annotate('T = %1.1f' % (time), xy=(0.05, 0.15), xycoords='axes fraction', size=14)
        ax.axis('off')

        fig.tight_layout()
        fig.savefig(prob + '_FN_T%1.0f_%1.0f' % (int(time), (time % 1) * 10))
        plt.close(fig)


if __name__ == "__main__":
    for t in np.linspace(3, 7, 9):
        main('res-files', 'pile_raft', t)
