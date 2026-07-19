#-------------------------------------------------------------------------------
# Name:        plot_excavation2D
# Purpose:     extract and plot the retaining wall's bending moment, shear
#              force and horizontal displacement profiles vs. depth, one
#              curve per construction stage
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
import C_BeamResults as C_Beam_Res
import C_NodalResults as C_NodRes

from result_utils import select_static_steps

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


def main(pathname, prob, steps_to_plot, step_labels,
         wall_mat=3, anchor_mats=(4, 5, 6), anchor_wall_x=30.001):
    """
    Plot the wall's bending moment, shear force and horizontal displacement
    vs. depth at each of `steps_to_plot`, one colored curve per stage
    (labelled from `step_labels`).

    `anchor_wall_x` picks out the wall-side node of each anchor row (as
    opposed to its far, fixed-length end) to place a depth marker on the
    bending-moment plot; it only needs changing if the wall isn't at x=0
    in your model.
    """
    file_path = os.path.join(pathname, prob)

    his = C_His.HistoryOfExecution(file_path)
    rcf = C_Rcf.RCF_info(file_path)
    mesh = C_Mesh.Mesh(file_path)

    solution_indices, _ = select_static_steps(his, steps=steps_to_plot)

    # get nodal and beam results:
    nd_rsl = C_NodRes.NodalResults(mesh, his, rcf)
    beam_rsl = C_Beam_Res.Beam_EleResults(mesh, his, rcf)

    # select beam elements of retaining wall:
    blist = mesh.get_list_of_elements('BEAMS', mat_filter=[wall_mat])

    # get indices for all nodes in model:
    node_indices = [node.index for node in mesh.nodes]

    fig = plt.figure(figsize=(12, 8))
    AX = [fig.add_subplot(1, 3, k + 1) for k in range(3)]
    handles = [[] for _ in range(3)]

    ylim = [1e10, -1e10]
    Mlim = [1e10, -1e10]
    Tlim = 0
    ulim = 0
    for kkt, sol in enumerate(solution_indices):
        tx = his.data[sol][his.DATA_TIME]

        # select MZ for current time step:
        moment = beam_rsl.get_moment_for_sel_elements([beam.index for beam in blist], sol, "Z")
        # select TY for current time step:
        shear = beam_rsl.get_force_for_sel_elements([beam.index for beam in blist], sol, "Y")
        # select ux for current time step:
        disp = nd_rsl.get_nodes_displacements(node_indices, sol, comp="X")

        for kk in range(3):
            handles[kk].append(0)

        for ke, beam in enumerate(blist):
            inel = [nd.index for nd in beam.get_nodes()]
            crds = beam.get_ele_coord()
            ylim[0] = min(ylim[0], min(crds[:, 1]))
            ylim[1] = max(ylim[1], max(crds[:, 1]))
            n = np.array([-crds[1, 1] + crds[0, 1], crds[1, 0] - crds[0, 0]])
            l = np.linalg.norm(n)
            n /= l
            v = n[0] * np.array([moment[ke] - 0.5 * l * shear[ke],
                                  moment[ke] + 0.5 * l * shear[ke]])
            Mlim[0] = min(Mlim[0], min(v))
            Mlim[1] = max(Mlim[1], max(v))
            handles[0][-1] = AX[0].plot(v, crds[:, 1], color=colors[kkt],
                                         label='T=%1.0f: %s' % (tx, step_labels[tx]))[0]

            v = n[0] * np.array([shear[ke], shear[ke]])
            Tlim = max([Tlim, max(v, key=abs)], key=abs)
            handles[1][-1] = AX[1].plot(v, crds[:, 1], color=colors[kkt],
                                         label='T=%1.0f: %s' % (tx, step_labels[tx]))[0]

            v = 1e3 * np.array([disp[inel[0] - 1], disp[inel[1] - 1]])
            ulim = max([ulim, max(v, key=abs)], key=abs)
            handles[2][-1] = AX[2].plot(v, crds[:, 1], color=colors[kkt],
                                         label='T=%1.0f: %s' % (tx, step_labels[tx]))[0]

    # mark the elevation of each anchor row on the bending-moment plot:
    tlist = [truss for mat in anchor_mats for truss in mesh.get_list_of_elements('TRUSSES', mat_filter=[mat])]
    for truss in tlist:
        crds = truss.get_ele_coord()
        if crds[0][0] < anchor_wall_x:
            AX[0].annotate('', xy=(0, crds[0][1]), xytext=(-20, 0),
                           textcoords='offset points', ha='right', va='center',
                           arrowprops=dict(facecolor='black', headwidth=6, width=2))

    for kax, ax in enumerate(AX):
        ax.set_ylabel(u'Depth [m]')
        ax.grid('on')
        ax.set_ylim(ylim)
        ax.legend(handles=handles[kax], labels=[v.get_label() for v in handles[kax]], loc='lower right')

    AX[0].set_xlabel('Bending moment [kNm/m]')
    AX[1].set_xlabel('Shear force [kN/m]')
    AX[2].set_xlabel('Horizontal displacement [mm]')

    AX[0].annotate('$M_{min}=%1.0f$ kNm/m\n$M_{max}=%1.0f$ kNm/m' % (Mlim[0], Mlim[1]),
                   xy=(0.05, 0.02), xycoords='axes fraction', rotation=90, va='bottom')
    AX[1].annotate('$|T|_{max}=%1.0f$ kN/m' % (Tlim),
                   xy=(0.05, 0.02), xycoords='axes fraction', rotation=90, va='bottom')
    AX[2].annotate('$u_{h,min}=%1.0f$ mm' % (ulim),
                   xy=(0.05, 0.02), xycoords='axes fraction', rotation=90, va='bottom')

    fig.tight_layout()
    fig.savefig('wall_' + prob)

    plt.close(fig)


if __name__ == "__main__":
    from demo_config import STEPS_TO_PLOT, STEP_LABELS
    main(sys.argv[1], sys.argv[2], STEPS_TO_PLOT, STEP_LABELS)
