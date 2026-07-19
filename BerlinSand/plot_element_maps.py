#-------------------------------------------------------------------------------
# Name:        plot_element_maps
# Purpose:     plot the deformed mesh, colored by element horizontal stress
#              (sigma_xx), with the retaining wall and anchors overlaid
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

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

import C_Mesh
import C_HistoryOfExecution as C_His
import C_Rcf_info as C_Rcf
import C_NodalResults as C_NodRes
import C_ContinuumResults

from result_utils import select_static_steps


def main(pathname, prob, step_to_plot, wall_mat=3, anchor_mats=(4, 5, 6), scale=100):
    """
    Plot the deformed mesh at `step_to_plot`, colored by element horizontal
    stress (sigma_xx), with the wall (material `wall_mat`) and anchors
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

    # get nodal and element results:
    nd_rsl = C_NodRes.NodalResults(mesh, his, rcf)
    vol_rsl = C_ContinuumResults.Continuum_EleResults(mesh, his, rcf)

    all_node_indices = [node.index for node in all_nodes]
    XY = [mesh.get_node(index).get_xyz()[:2] for index in all_node_indices]

    res = [nd_rsl.get_nodes_displacements(all_node_indices, solution_indices[-1], comp="X"),
           nd_rsl.get_nodes_displacements(all_node_indices, solution_indices[-1], comp="Y")]
    XYdisp = [(xy[0] + scale * res[0][index], xy[1] + scale * res[1][index]) for index, xy in enumerate(XY)]

    all_vol_indices = [element.index for element in vollist]

    # STRESSES result per element: one stress vector per Gauss point; [0][0]
    # takes Gauss point 0's first component, i.e. sigma_xx (matching the
    # "Horizontal stress" label below).
    stresses = vol_rsl.get_rsl_for_sel_elements(all_vol_indices, solution_indices[-1], 'STRESSES')

    patches = []
    cvect = []
    for ke, ele in enumerate(vollist):
        inel = ele.nodes
        crds = [XYdisp[index - 1] for index in inel]
        patches.append(Polygon(crds))
        cvect.append(stresses[ke][0][0])

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(1, 1, 1)

    pc = PatchCollection(patches, cmap='viridis', edgecolors='none', antialiased=False)
    pc.set_array(cvect)
    ax.add_collection(pc)

    cb = fig.colorbar(pc)
    cb.set_label('Horizontal stress [kPa]')

    for beam in beams:
        nodes = beam.get_nodes()
        xb = [XYdisp[nd.index - 1][0] for nd in nodes]
        yb = [XYdisp[nd.index - 1][1] for nd in nodes]
        ax.plot(xb, yb, color='orange', lw=2)

    for trusses0 in trusses:
        for truss in trusses0:
            nodes = truss.get_nodes()
            xb = [XYdisp[nd.index - 1][0] for nd in nodes]
            yb = [XYdisp[nd.index - 1][1] for nd in nodes]
            ax.plot(xb, yb, color='orange', lw=2)

    ax.axis('off')
    ax.set_aspect('equal')

    fig.tight_layout()
    fig.savefig('Sxx_map_' + prob + '_T=%1.0f' % (step_to_plot))

    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]))
