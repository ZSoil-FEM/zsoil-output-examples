#-------------------------------------------------------------------------------
# Name:        plot_nodal_maps
# Purpose:     plot the deformed mesh, colored by resultant nodal
#              displacement, with the retaining wall and anchors overlaid
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
import matplotlib.tri as tri

import C_Mesh
import C_HistoryOfExecution as C_His
import C_Rcf_info as C_Rcf
import C_NodalResults as C_NodRes

from result_utils import select_static_steps


def main(pathname, prob, step_to_plot, wall_mat=3, anchor_mats=(4, 5, 6), scale=100):
    """
    Plot the deformed mesh at `step_to_plot`, colored by resultant nodal
    displacement magnitude, with the wall (material `wall_mat`) and anchors
    (materials `anchor_mats`) overlaid in their deformed position.
    """
    file_path = os.path.join(pathname, prob)

    his = C_His.HistoryOfExecution(file_path)
    rcf = C_Rcf.RCF_info(file_path)
    mesh = C_Mesh.Mesh(file_path)

    solution_indices, _ = select_static_steps(his, steps=[step_to_plot])

    # get indices for all nodes:
    all_nodes = mesh.get_list_of_nodes()

    # select beam elements of retaining wall:
    beams = mesh.get_list_of_elements('BEAMS', mat_filter=[wall_mat])

    # select truss elements sorted by material:
    trusses = [mesh.get_list_of_elements('TRUSSES', mat_filter=[mat]) for mat in anchor_mats]

    # get all volume elements active at the plotted step:
    vollist = mesh.get_list_of_elements('VOLUMICS', time=step_to_plot, only_active=True)

    # get nodal results:
    nd_rsl = C_NodRes.NodalResults(mesh, his, rcf)

    triangles = []
    for ele in vollist:
        inel = ele.nodes
        triangles.append([inel[kv] - 1 for kv in [0, 1, 2]])
        triangles.append([inel[kv] - 1 for kv in [0, 2, 3]])

    all_node_indices = [node.index for node in all_nodes]

    XY = [mesh.get_node(index).get_xyz()[:2] for index in all_node_indices]
    res = [nd_rsl.get_nodes_displacements(all_node_indices, solution_indices[-1], comp="X"),
           nd_rsl.get_nodes_displacements(all_node_indices, solution_indices[-1], comp="Y")]

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(1, 1, 1)

    x = np.array([XY[kn][0] + res[0][kn] * scale for kn in range(len(res[0]))])
    y = np.array([XY[kn][1] + res[1][kn] * scale for kn in range(len(res[0]))])
    triang = tri.Triangulation(x, y, triangles)

    val = [(res[0][kn] ** 2 + res[1][kn] ** 2) ** 0.5 * 1e3 for kn in range(len(res[0]))]
    a = ax.tripcolor(triang, val)
    cb = fig.colorbar(a)
    cb.set_label('Absolute displacement [mm]')

    for beam in beams:
        nodes = beam.get_nodes()
        xb = [nd.get_xyz()[0] + res[0][nd.index - 1] * scale for nd in nodes]
        yb = [nd.get_xyz()[1] + res[1][nd.index - 1] * scale for nd in nodes]
        ax.plot(xb, yb, color='orange', lw=2)

    for trusses0 in trusses:
        for truss in trusses0:
            nodes = truss.get_nodes()
            xb = [nd.get_xyz()[0] + res[0][nd.index - 1] * scale for nd in nodes]
            yb = [nd.get_xyz()[1] + res[1][nd.index - 1] * scale for nd in nodes]
            ax.plot(xb, yb, color='orange', lw=2)

    ax.axis('off')
    ax.set_aspect('equal')

    fig.tight_layout()
    fig.savefig('disp_map_' + prob + '_T=%1.0f' % (step_to_plot))

    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]))
