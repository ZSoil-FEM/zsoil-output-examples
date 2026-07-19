#-------------------------------------------------------------------------------
# Name:        plot_piles_diagrams
# Purpose:     extract and plot normal force and mobilized friction along
#              each pile, including the *available* (maximum) shaft
#              friction recovered from the interface stress level, and the
#              mobilized tip resistance - quantities PostPro does not plot
#              directly
#
# Author:      Matthias Preisig
# Created:     2026
# Copyright:   (c) GeoDev Sarl
# Licence:     <your licence>
#-------------------------------------------------------------------------------
# -*- coding: cp1252 -*-
import os
import sys
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from zsoilpy_env import ensure_zsoilpy_on_path
ensure_zsoilpy_on_path()

import matplotlib.pyplot as plt

import C_Mesh
import C_HistoryOfExecution as C_His
import C_Rcf_info as C_Rcf
import C_BeamResults as C_Beam_Res
import C_Contact_3D_Results as C_Cnt_Res


def _group_by_pile(elements, all_nodes, precision):
    """
    Group elements into piles by rounding the (x, z) footprint of their
    first node. There is no explicit "pile" grouping in the mesh data, so
    this is how beams and their shaft/tip interfaces belonging to the same
    physical pile get matched up; `precision` is the rounding tolerance (in
    model length units) and may need loosening if piles aren't perfectly
    vertical/aligned in your mesh.
    """
    grouped = defaultdict(list)
    for ele in elements:
        pt = all_nodes[ele.nodes[0] - 1].get_xyz()
        pile_key = (round(pt[0], precision), round(pt[2], precision))
        grouped[pile_key].append(ele)
    return grouped


def main(pathname, prob, steps=None, pile_mat=9, shaft_mat=8, tip_mat=10, precision=3):
    """
    Plot, for each pile and each requested step (every converged static
    step if `steps` is None): the axial force N(z) from the pile beams
    (material `pile_mat`), and the shaft friction from the shaft interface
    (material `shaft_mat`) - both mobilized (q_s,mob, directly from
    STRESSES) and available (q_s,max, recovered as q_s,mob / STR_LEVEL,
    since STR_LEVEL is the ratio of mobilized to available shear
    resistance). The mobilized tip resistance from the tip interface
    (material `tip_mat`) is annotated on each friction subplot.
    """
    file_path = os.path.join(pathname, prob)

    his = C_His.HistoryOfExecution(file_path)
    rcf = C_Rcf.RCF_info(file_path)
    mesh = C_Mesh.Mesh(file_path)

    # select solutions without stability driver, optionally restricted to
    # specific step times:
    solution_indices = []
    for sol_ind in his.give_converged_solutions():
        time = his.data[sol_ind][his.DATA_TIME]
        SF = his.data[sol_ind][his.DATA_SF]
        if SF == 0 and (steps is None or time in steps):
            solution_indices.append(sol_ind)

    # get beam and interface results:
    beam_rsl = C_Beam_Res.Beam_EleResults(mesh, his, rcf)
    cnt_rsl = C_Cnt_Res.Contact_3D_Results(mesh, his, rcf)

    # get indices for all nodes:
    all_nodes = mesh.get_list_of_nodes()

    # select pile beams and shaft/tip interface elements:
    beams = mesh.get_list_of_elements('BEAMS', mat_filter=[pile_mat])
    beaminds = [ele.index for ele in beams]

    shaftlist = mesh.get_list_of_elements('CONTACT', mat_filter=[shaft_mat])
    shaftinds = [ele.index for ele in shaftlist]

    tiplist = mesh.get_list_of_elements('CONTACT', mat_filter=[tip_mat])
    tipinds = [ele.index for ele in tiplist]

    # identify piles and match beams and contact elements with piles:
    pile_to_beams = _group_by_pile(beams, all_nodes, precision)
    pile_to_shafts = _group_by_pile(shaftlist, all_nodes, precision)
    pile_to_tips = _group_by_pile(tiplist, all_nodes, precision)
    piles = sorted(pile_to_beams.keys())
    nPiles = len(piles)

    # create a figure for each solution
    for sol in solution_indices:
        time = his.data[sol][his.DATA_TIME]

        shaft_str = cnt_rsl.get_rsl_vec_for_sel_elements(shaftinds, sol, 'STRESSES')
        str_level = cnt_rsl.get_rsl_vec_for_sel_elements(shaftinds, sol, 'STR_LEVEL')
        tip_str = cnt_rsl.get_rsl_vec_for_sel_elements(tipinds, sol, 'STRESSES')
        forces = beam_rsl.get_forces_vec_for_sel_elements(beaminds, sol)

        fig = plt.figure(figsize=(20, 12))
        fig.text(0.01, 0.01, prob, size=8)
        AX = [[], []]
        handles = [0, 0]

        # create a subplot pair for each pile
        for kpile in range(nPiles):
            # one subplot for normal force diagram:
            AX[0].append(fig.add_subplot(2, nPiles, kpile + 1))
            # second subplot for mobilized/available friction along piles:
            AX[1].append(fig.add_subplot(2, nPiles, nPiles + kpile + 1))

            pile_key = piles[kpile]

            ax0 = AX[0][kpile]
            pile_beams = pile_to_beams[pile_key]
            for beam in pile_beams:
                kke = beaminds.index(beam.index)
                pt0 = all_nodes[beam.nodes[0] - 1].get_xyz()
                pt1 = all_nodes[beam.nodes[1] - 1].get_xyz()

                # NX is the 1st component, value is constant over one beam
                ax0.plot([forces[kke][0], forces[kke][0]],
                         [pt0[1], pt1[1]], 'k')

            ax1 = AX[1][kpile]
            pile_shafts = pile_to_shafts[pile_key]
            for cnt in pile_shafts:
                kke = shaftinds.index(cnt.index)
                pt0 = all_nodes[cnt.nodes[0] - 1].get_xyz()
                pt1 = all_nodes[cnt.nodes[1] - 1].get_xyz()

                # plot the mobilized shaft friction:
                handles[0] = ax1.plot([shaft_str[kke][0][0], shaft_str[kke][1][0]],
                                      [pt0[1], pt1[1]], 'b', label=u'$q_{s,mob}$')[0]
                # STR_LEVEL is the ratio of mobilized to available shear
                # stress (0-1); dividing it out recovers the available
                # (maximum) shaft friction - not a quantity PostPro plots
                # directly. A stress level of 0 means nothing is mobilized
                # yet, so the available friction is just the mobilized
                # (zero) stress itself.
                if min(str_level[kke]) > 0:
                    handles[1] = ax1.plot([abs(shaft_str[kke][0][0]) / str_level[kke][0],
                                           abs(shaft_str[kke][1][0]) / str_level[kke][1]],
                                          [pt0[1], pt1[1]], 'k--', label=u'$q_{s,max}$')[0]
                else:
                    handles[1] = ax1.plot([abs(shaft_str[kke][0][0]),
                                           abs(shaft_str[kke][1][0])],
                                          [pt0[1], pt1[1]], 'k--', label=u'$q_{s,max}$')[0]

            # write the mobilized tip resistance in the qs-mob-plots:
            pile_tips = pile_to_tips[pile_key]
            ax1.annotate('$\\sigma_p=$%1.0f kPa' % (tip_str[tipinds.index(pile_tips[0].index)][0][2]),
                         xy=(0.5, 0.06), xycoords='axes fraction', ha='center', rotation=90)
            # write the x-z coordinates in the NX-plots
            ax0.annotate('(%1.1f,%1.1f)' % (pile_key[0], pile_key[1]),
                         xy=(0.1, 0.92), xycoords='axes fraction')

        for ax in AX[0]:
            ax.grid(visible=True, which='both')
            ax.set_xlabel(u'$F_N$ [kN]')
        for ax in AX[1]:
            ax.grid(visible=True, which='both')
            ax.set_xlabel(u'$q_{s,mob}$ [kPa]')
        AX[1][-1].legend(handles=handles)

        fig.tight_layout()
        fig.subplots_adjust(wspace=0)
        fig.savefig(prob + '_frott_T%1.0f_%1.0f' % (int(time), (time % 1) * 10))
        plt.close(fig)


if __name__ == "__main__":
    main('res-files', 'pile_raft')
