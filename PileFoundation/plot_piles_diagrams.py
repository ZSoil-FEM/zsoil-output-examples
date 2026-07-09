#-------------------------------------------------------------------------------
# Name:        plot_piles_diagrams
# Purpose:     extract and plot normal force and mobilized friction along piles,
#               plot diagrams in separate subplots, indicate also mobilized
#               tip resistance
#               a 2D projection with coloring according to plastic codes of
#               interface elements
#
# Author:      Matthias Preisig
#
# Created:     2026
# Copyright:   (c) GeoDev Sàrl
# Licence:     <your licence>
#-------------------------------------------------------------------------------
# -*- coding: cp1252 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.patches import Polygon
import os,sys
from collections import defaultdict

__temp_path = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), r".."))
__temp_path = os.path.normpath('C:\Program Files\ZSoil\Tools v2026\ZSoilPy3')
sys.path.append(__temp_path) if __temp_path not in sys.path else None

import C_Mesh
import C_HistoryOfExecution as C_His
import C_Rcf_info as C_Rcf
import C_EleResults as C_Ele_Res
import C_BeamResults as C_Beam_Res
import C_NodalResults as C_NodRes
import C_ContinuumResults
import C_Contact_3D_Results as C_Cnt_Res

def main(pathname, prob, step_to_plot):

    file_path = os.path.join(pathname, prob)

    his = C_His.HistoryOfExecution(file_path)
    rcf = C_Rcf.RCF_info(file_path)
    mesh = C_Mesh.Mesh(file_path)

    all_converged_solutions = his.give_converged_solutions()

    # select solutions without stability driver and integer time instances:
    solution_indices = []
    all_solution_indices = []
    for sol_ind in all_converged_solutions:
        time = his.data[sol_ind][his.DATA_TIME]
        SF = his.data[sol_ind][his.DATA_SF]
        if SF==0:
            all_solution_indices.append(sol_ind)
##            if time == step_to_plot:
            solution_indices.append(sol_ind)

    # get nodal, beam and interface results:
    nd_rsl = C_NodRes.NodalResults(mesh, his, rcf)
    beam_rsl = C_Beam_Res.Beam_EleResults(mesh, his, rcf)
    cnt_rsl = C_Cnt_Res.Contact_3D_Results(mesh, his, rcf)

    # get indices for all nodes:
    all_nodes = mesh.get_list_of_nodes()

    # select beam elements:
    beams = mesh.get_list_of_elements('BEAMS',mat_filter=[9])
    beaminds = [ele.index for ele in beams]

    # select interface elements of pile shafts:
    shaftlist = mesh.get_list_of_elements('CONTACT',mat_filter=[8])
    shaftinds = [ele.index for ele in shaftlist]

    # select interface elements of pile tips:
    tiplist = mesh.get_list_of_elements('CONTACT',mat_filter=[10])
    tipinds = [ele.index for ele in tiplist]

    # identify piles and match beams and contact elements with piles:
    piles = set()
    pile_to_beams = defaultdict(list)
    for beam in beams:
        pt = all_nodes[beam.nodes[0] - 1].get_xyz()

        precision = 3
        x_coord = round(pt[0], precision)
        z_coord = round(pt[2], precision)

        pile_key = (x_coord, z_coord)
        pile_to_beams[pile_key].append(beam)
        piles.add(pile_key)
    pile_to_shafts = defaultdict(list)
    for shaft in shaftlist:
        pt = all_nodes[shaft.nodes[0] - 1].get_xyz()

        precision = 3
        x_coord = round(pt[0], precision)
        z_coord = round(pt[2], precision)

        pile_key = (x_coord, z_coord)
        pile_to_shafts[pile_key].append(shaft)
    pile_to_tips = defaultdict(list)
    for tip in tiplist:
        pt = all_nodes[tip.nodes[0] - 1].get_xyz()

        precision = 3
        x_coord = round(pt[0], precision)
        z_coord = round(pt[2], precision)

        pile_key = (x_coord, z_coord)
        pile_to_tips[pile_key].append(tip)
    piles = sorted(list(piles))
    nPiles = len(piles)

    # create a figure for each solution
    for sol in solution_indices:
        time = his.data[sol][his.DATA_TIME]

        shaft_str = cnt_rsl.get_rsl_vec_for_sel_elements(shaftinds, sol, 'STRESSES')
##        pla_code = cnt_rsl.get_rsl_vec_for_sel_elements(shaftinds, sol, 'PLA_CODE')
        str_level = cnt_rsl.get_rsl_vec_for_sel_elements(shaftinds, sol, 'STR_LEVEL')
        tip_str = cnt_rsl.get_rsl_vec_for_sel_elements(tipinds, sol, 'STRESSES')
        forces = beam_rsl.get_forces_vec_for_sel_elements(beaminds, sol)
    

        fig = plt.figure(figsize=(20,12))
        fig.text(0.01,0.01,prob,size=8)
        AX = [[],[]]
        handles = [0,0]

        # create a subplot for each
        for kpile in range(nPiles):
            # one subplot for normal force diagram:
            AX[0].append(fig.add_subplot(2,nPiles,kpile+1))
            # second subplot for mobilized friction along piles:
            AX[1].append(fig.add_subplot(2,nPiles,nPiles+kpile+1))

            pile_key = piles[kpile]

            ax0 = AX[0][kpile]
            pile_beams = pile_to_beams[pile_key]
            for kk,beam in enumerate(pile_beams):
                kke = beaminds.index(beam.index)
                pt0 = all_nodes[beam.nodes[0] - 1].get_xyz()
                pt1 = all_nodes[beam.nodes[1] - 1].get_xyz()

                # NX is the 1st component, value is constant over one beam
                ax0.plot([forces[kke][0],forces[kke][0]],
                         [pt0[1],pt1[1]],'k')

            ax1 = AX[1][kpile]
            pile_shafts = pile_to_shafts[pile_key]
            for kk,cnt in enumerate(pile_shafts):
                kke = shaftinds.index(cnt.index)
                pt0 = all_nodes[cnt.nodes[0] - 1].get_xyz()
                pt1 = all_nodes[cnt.nodes[1] - 1].get_xyz()

                # plot the mobilized stress
                handles[0] = ax1.plot([shaft_str[kke][0][0],shaft_str[kke][1][0]],
                                      [pt0[1],pt1[1]],'b',label=u'$q_{s,mob}$')[0]
                # plot the mobilized stress divided by the stress level, representing the
                # shaft resistance
                if min(str_level[kke])>0:
                    handles[1] = ax1.plot([abs(shaft_str[kke][0][0])/str_level[kke][0],
                                           abs(shaft_str[kke][1][0])/str_level[kke][1]],
                                          [pt0[1],pt1[1]],'k--',label=u'$q_{s,max}$')[0]
                else:
                    handles[1] = ax1.plot([abs(shaft_str[kke][0][0]),
                                           abs(shaft_str[kke][1][0])],
                                          [pt0[1],pt1[1]],'k--',label=u'$q_{s,max}$')[0]

            # write the mobilized tip resistance in the qs-mob-plots:
            pile_tips = pile_to_tips[pile_key]
            ax1.annotate('$\sigma_p=$%1.0f kPa'%(tip_str[tipinds.index(pile_tips[0].index)][0][2]),
                         xy=(0.5,0.06),xycoords='axes fraction',ha='center',rotation=90)
            # write the x-z coordinates in the NX-plots
            ax0.annotate('(%1.1f,%1.1f)'%(pile_key[0],pile_key[1]),
                         xy=(0.1,0.92),xycoords='axes fraction')

        for ax in AX[0]:
            ax.grid(visible=True,which='both')
            ax.set_xlabel(u'$F_N$ [kN]')
        for ax in AX[1]:
            ax.grid(visible=True,which='both')
            ax.set_xlabel(u'$q_{s,mob}$ [kPa]')
        AX[1][-1].legend(handles=handles)

        fig.tight_layout()
        fig.subplots_adjust(wspace=0)
        fig.savefig(prob+'_frott_T%1.0f_%1.0f'%(int(time),(time%1)*10))
        plt.close(fig)



if __name__ == "__main__":
##    for t in np.linspace(3,7,9):
    main('res-files','pile_raft',7)
##    main(sys.argv[1], sys.argv[2])

