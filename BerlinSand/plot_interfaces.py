#-------------------------------------------------------------------------------
# Name:        plot_interfaces
# Purpose:     extract and plot the normal earth pressure behind the
#              retaining wall, and the mobilized friction along each anchor
#              row, vs. depth
#
# Author:      Matthias Preisig
# Created:     2026
# Copyright:   (c) GeoDev Sarl
#-------------------------------------------------------------------------------
# -*- coding: cp1252 -*-
import os
import sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from zsoilpy_env import ensure_zsoilpy_on_path
ensure_zsoilpy_on_path()

import numpy as np
import matplotlib.pyplot as plt

import C_Mesh
import C_HistoryOfExecution as C_His
import C_Rcf_info as C_Rcf
import C_Contact_2D_Results as C_Cnt_Res

from result_utils import select_static_steps

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


def _anchor_row(y, depth_bounds):
    """
    Classify an anchor interface element into a row index from its
    y-coordinate. `depth_bounds` is a descending sequence of elevations,
    e.g. (-15, -20), which splits anchors into 3 rows: y > -15,
    -20 < y <= -15, and y <= -20.
    """
    for k, bound in enumerate(depth_bounds):
        if y > bound:
            return k
    return len(depth_bounds)


def main(pathname, prob, steps_to_plot, step_labels,
         wall_interface_mat=7, anchor_interface_mat=9, anchor_depth_bounds=(-15, -20)):
    """
    Plot the normal interface stress behind the wall (material
    `wall_interface_mat`) and the mobilized shear stress along the anchors'
    fixed-length zone (material `anchor_interface_mat`), vs. depth, at each
    of `steps_to_plot`.

    Unlike the anchor trusses (one material per row), the anchor interface
    elements in this model all share a single material, so rows are told
    apart by elevation instead, using `anchor_depth_bounds`.
    """
    file_path = os.path.join(pathname, prob)

    his = C_His.HistoryOfExecution(file_path)
    rcf = C_Rcf.RCF_info(file_path)
    mesh = C_Mesh.Mesh(file_path)

    solution_indices, _ = select_static_steps(his, steps=steps_to_plot)

    # get interface results:
    cnt_rsl = C_Cnt_Res.Contact_2D_Results(mesh, his, rcf)

    # select interface elements of retaining wall:
    wilist = mesh.get_list_of_elements('CONTACT', mat_filter=[wall_interface_mat])

    # select interface elements of anchors:
    ancilist = mesh.get_list_of_elements('CONTACT', mat_filter=[anchor_interface_mat])

    nAnchors = len(anchor_depth_bounds) + 1

    fig_wall = plt.figure(figsize=(8, 8))
    ax_wall = fig_wall.add_subplot(1, 1, 1)

    fig_anc = plt.figure(figsize=(8, 10))
    AX_anc = [fig_anc.add_subplot(nAnchors, 1, k + 1) for k in range(nAnchors)]

    h_wall = []
    h_anc = []

    ylim = [1e10, -1e10]
    for kkt, sol in enumerate(solution_indices):
        tx = his.data[sol][his.DATA_TIME]

        # select normal stress behind the wall for current time step:
        sigN = cnt_rsl.get_rsl_vec_for_sel_elements([ele.index for ele in wilist], sol, 'STRESSES', comp="YY")
        # select mobilized shear (friction) in the anchors for current time step:
        tau = cnt_rsl.get_rsl_vec_for_sel_elements([ele.index for ele in ancilist], sol, 'STRESSES', comp="XY")

        h_wall.append(0)
        h_anc.append(0)

        for kk, wall_element in enumerate(wilist):
            crds = wall_element.get_ele_coord()
            ylim[0] = min(ylim[0], min(crds[:, 1]))
            ylim[1] = max(ylim[1], max(crds[:, 1]))
            n = np.array([-crds[1, 1] + crds[0, 1], crds[1, 0] - crds[0, 0]])
            l = np.linalg.norm(n)
            n /= l
            v = n[0] * np.array(sigN[kk])
            h_wall[-1] = ax_wall.plot(v, crds[:, 1], color=colors[kkt],
                                       label='T=%1.0f: %s' % (tx, step_labels[tx]))[0]

        for kk, anc_interf in enumerate(ancilist):
            crds = anc_interf.get_ele_coord()
            p0 = crds[0, :]
            p1 = crds[1, :]
            ax_anc = AX_anc[_anchor_row(p0[1], anchor_depth_bounds)]
            h_anc[-1] = ax_anc.plot([p0[0], p1[0]], tau[kk], color=colors[kkt],
                                     label='T=%1.0f: %s' % (tx, step_labels[tx]))[0]

    ax_wall.set_ylabel(u'Depth [m]')
    ax_wall.grid('on')
    ax_wall.set_ylim(ylim)
    ax_wall.legend(handles=h_wall, labels=[v.get_label() for v in h_wall])

    for ax_anc in AX_anc:
        ax_anc.grid('on')
        ax_anc.set_xlabel('x-coordinate [m]')
        ax_anc.set_ylabel('Mobilized friction [kPa]')
        ax_anc.legend(handles=h_anc, labels=[v.get_label() for v in h_anc])

    fig_wall.tight_layout()
    fig_wall.savefig('earth_pressures_' + prob)
    plt.close(fig_wall)

    fig_anc.tight_layout()
    fig_anc.savefig('anchors_' + prob)
    plt.close(fig_anc)


if __name__ == "__main__":
    from demo_config import STEPS_TO_PLOT, STEP_LABELS
    main(sys.argv[1], sys.argv[2], STEPS_TO_PLOT, STEP_LABELS)
