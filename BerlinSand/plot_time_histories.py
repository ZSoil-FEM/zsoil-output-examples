# -*- coding: cp1252 -*-
import numpy as np
import matplotlib.pyplot as plt
import os,sys

__temp_path = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), r".."))
__temp_path = os.path.normpath('C:\Program Files\ZSoil\Tools v2026\ZSoilPy3')
sys.path.append(__temp_path) if __temp_path not in sys.path else None

import C_Mesh
import C_HistoryOfExecution as C_His
import C_Rcf_info as C_Rcf
import C_EleResults as C_Ele_Res
import C_TrussResults as C_Truss_Res
import C_BeamResults as C_Beam_Res
import C_NodalResults as C_NodRes

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

def main(pathname, prob):
    
    file_path = os.path.join(pathname, prob)

    his = C_His.HistoryOfExecution(file_path)
    rcf = C_Rcf.RCF_info(file_path)
    mesh = C_Mesh.Mesh(file_path)

    all_solution_indices = his.give_converged_solutions()

    # select solutions without stability driver and integer time instances:
    solution_indices = []
    for sol_ind in all_solution_indices:
        time = his.data[sol_ind][his.DATA_TIME]
        SF = his.data[sol_ind][his.DATA_SF]
        if SF==0:# and abs(time%1)<1e-3:
            solution_indices.append(sol_ind)

    # get nodal, beam and truss results:
    nd_rsl = C_NodRes.NodalResults(mesh, his, rcf)
    beam_rsl = C_Beam_Res.Beam_EleResults(mesh, his, rcf)
    truss_rsl = C_Truss_Res.Truss_EleResults(mesh, his, rcf)

    # select beam elements of retaining wall:
    blist = mesh.get_list_of_elements('BEAMS',mat_filter=[3])
    beam_indices = [beam.index for beam in blist]

    # select truss elements sorted by material:
    tlists = [mesh.get_list_of_elements('TRUSSES',mat_filter=[mat]) for mat in [4,5,6]]

    # get indices for all nodes in retaining wall:
    wall_node_indices = [node.index for node in mesh.get_list_of_nodes(ele_filter=blist)]

    # get indices for all nodes in retaining wall:
    surface_node_indices = [node.index for node in mesh.get_list_of_nodes(zoom_filter=[0,[[30,200],[-0.1,0.1]]])]

    # time instances:
    TX = []

    # maximum horizontal displacement of wall:
    UX_max = []
    # maximum settlement behind wall:
    UY_back_max = []
    # maximum positive bending moment:
    MZ_pos_max = []
    # maximum negative bending moment:
    MZ_neg_max = []
    # anchor force for all anchors:
    nAnchors = 3
    N_anc = [[] for _ in range(nAnchors)]

    FIG = [plt.figure(figsize=(6,3)) for _ in range(2)]
    AX = [fig.add_subplot(1,1,1) for fig in FIG]
    fig_M = plt.figure(figsize=(6,6))
    AX_M = [fig_M.add_subplot(2,1,k+1) for k in range(2)]
    fig_anc = plt.figure(figsize=(6,8))
    AX_anc = [fig_anc.add_subplot(nAnchors,1,k+1) for k in range(nAnchors)]
    
    for sol in solution_indices:
        tx = his.data[sol][his.DATA_TIME]
        TX.append(tx)

        # select MZ for current time step:
        moment = beam_rsl.get_moment_for_sel_elements(beam_indices, sol, "Z")
        MZ_pos_max.append(max(moment))
        MZ_neg_max.append(min(moment))
        # select minimum wall displacement for current time step:
        UX_max.append(1e3*min(nd_rsl.get_nodes_displacements(wall_node_indices, sol, comp="X")))
        # select minimum settlement for current time step:
        UY_back_max.append(1e3*min(nd_rsl.get_nodes_displacements(surface_node_indices, sol, comp="Y")))
        # select maximum anchor force for current time step:
        for ka in range(nAnchors):
            indices = [truss.index for truss in tlists[ka]]
            N_anc[ka].append(max(truss_rsl.get_force_for_sel_elements(indices, sol)))

    AX[0].plot(TX,UX_max)
    AX[1].plot(TX,UY_back_max)
    AX_M[0].plot(TX,MZ_pos_max)
    AX_M[1].plot(TX,MZ_neg_max)
    for ka in range(nAnchors):
        AX_anc[ka].plot(TX,N_anc[ka])

    AXES = AX
    AXES.extend(AX_M)
    AXES.extend(AX_anc)
    for kax,ax in enumerate(AXES):
        ax.set_xlabel(u'Time step [-]')
        ax.grid('on')
        ax.set_xlim([min(TX),max(TX)])

    AX[0].set_ylabel('Maximum wall displacement [mm]')
    AX[1].set_ylabel('Maximum settlement behind wall [mm]')
    AX_M[0].set_xlabel('Maximum positive bending moment [kNm/m]')
    AX_M[1].set_xlabel('Maximum negative bending moment [kNm/m]')
    AX_anc[0].set_xlabel('Force in anchor row 1 [kN]')
    AX_anc[1].set_xlabel('Force in anchor row 2 [kN]')
    AX_anc[2].set_xlabel('Force in anchor row 3 [kN]')

    FIGS = FIG
    FIG.append(fig_M)
    FIG.append(fig_anc)
    for fig in FIGS:
        fig.tight_layout()
    FIG[0].savefig('wall_disp_evol_'+prob)
    FIG[1].savefig('settlement_evol_'+prob)
    fig_M.savefig('bending_moment_evol_'+prob)
    fig_anc.savefig('anchor_force_evol_'+prob)

    for fig in FIGS:
        plt.close(fig)
        


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])





