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
import C_Contact_2D_Results as C_Cnt_Res
import C_NodalResults as C_NodRes

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


def main(pathname, prob, steps_to_plot, step_labels):
    file_path = os.path.join(pathname, prob)

    his = C_His.HistoryOfExecution(file_path)
    rcf = C_Rcf.RCF_info(file_path)
    mesh = C_Mesh.Mesh(file_path)

    all_solution_indices = his.give_converged_solutions()

    # select solutions without stability driver and specific steps:
    solution_indices = []
    for sol_ind in all_solution_indices:
        time = his.data[sol_ind][his.DATA_TIME]
        SF = his.data[sol_ind][his.DATA_SF]
        if SF==0 and time in steps_to_plot:
            solution_indices.append(sol_ind)

    # get nodal, beam and truss results:
    nd_rsl = C_NodRes.NodalResults(mesh, his, rcf)
    beam_rsl = C_Beam_Res.Beam_EleResults(mesh, his, rcf)
    truss_rsl = C_Truss_Res.Truss_EleResults(mesh, his, rcf)
    cnt_rsl = C_Cnt_Res.Contact_2D_Results(mesh, his, rcf)

    # select interface elements of retaining wall:
    wilist = mesh.get_list_of_elements('CONTACT',mat_filter=[7])

    # select interface elements of anchors:
    ancilist = mesh.get_list_of_elements('CONTACT',mat_filter=[9])

    # get indices for all nodes in model:
    node_indices = [node.index for node in mesh.nodes]

    fig_wall = plt.figure(figsize=(8,8))
    ax_wall = fig_wall.add_subplot(1,1,1)

    fig_anc = plt.figure(figsize=(8,10))
    nAnchors = 3
    AX_anc = [fig_anc.add_subplot(nAnchors,1,k+1) for k in range(nAnchors)]

    h_wall = []
    h_anc = []
    
    ylim = [1e10,-1e10]
    kkt = -1
    for sol in solution_indices:
        tx = his.data[sol][his.DATA_TIME]
        if abs(tx%1)<1e-3 and tx>2:
            kkt += 1
            # select sig_N for current time step:
            sigN = cnt_rsl.get_rsl_vec_for_sel_elements([ele.index for ele in wilist],sol,'STRESSES',comp="YY")
            wSL = cnt_rsl.get_rsl_vec_for_sel_elements([ele.index for ele in wilist],sol,'STR_LEVEL',comp="")
            wPC = cnt_rsl.get_rsl_vec_for_sel_elements([ele.index for ele in wilist],sol,'PLA_CODE',comp="")

            # select mobilized traction in anchors for current time step:
            tau = cnt_rsl.get_rsl_vec_for_sel_elements([ele.index for ele in ancilist],sol,'STRESSES',comp="XY")
            SL = cnt_rsl.get_rsl_vec_for_sel_elements([ele.index for ele in ancilist],sol,'STR_LEVEL',comp="")
            PC = cnt_rsl.get_rsl_vec_for_sel_elements([ele.index for ele in ancilist],sol,'PLA_CODE',comp="")

            h_wall.append(0)
            h_anc.append(0)

            for kk,wall_element in enumerate(wilist):
                inel = [nd.index for nd in wall_element.get_nodes()]
                crds = wall_element.get_ele_coord()
                ylim[0] = min(ylim[0],min(crds[:,1]))
                ylim[1] = max(ylim[1],max(crds[:,1]))
                n = np.array([-crds[1,1]+crds[0,1],crds[1,0]-crds[0,0]])
                l = np.linalg.norm(n)
                n /= l
                v = n[0]*np.array(sigN[kk])
                h_wall[-1] = ax_wall.plot(v,crds[:,1],color=colors[kkt],label='T=%1.0f: %s'%(tx,step_labels[tx]))[0]

            sc = 5e-3
            for kk,anc_interf in enumerate(ancilist):
                inel = [nd.index for nd in anc_interf.get_nodes()]
                crds = anc_interf.get_ele_coord()
                p0 = crds[0,:]
                p1 = crds[1,:]
                n = np.cross(p1-p0,np.array([0,0,1]))
                l = np.linalg.norm(n)
                n /= l
                v0 = p0+n*sc*tau[kk][0]
                v1 = p1+n*sc*tau[kk][1]
                if p0[1]>-15:
                    ax_anc = AX_anc[0]
                elif p0[1]>-20:
                    ax_anc = AX_anc[1]
                else:
                    ax_anc = AX_anc[2]
                h_anc[-1] = ax_anc.plot([p0[0],p1[0]],tau[kk],color=colors[kkt],label='T=%1.0f: %s'%(tx,step_labels[tx]))[0]                


    if True:
        ax_wall.set_ylabel(u'Depth [m]')
        ax_wall.grid('on')
        ax_wall.set_ylim(ylim)
        ax_wall.legend(handles=h_wall,labels=[v.get_label() for v in h_wall])
    if True:
        for ax_anc in AX_anc:
            ax_anc.grid('on')
            ax_anc.set_xlabel('x-coordinate [m]')
            ax_anc.set_ylabel('Mobilized friction [kPa]')
            ax_anc.legend(handles=h_anc,labels=[v.get_label() for v in h_anc])

    fig_wall.tight_layout()
    fig_wall.savefig('earth_pressures_'+prob)

    plt.close(fig_wall)

    fig_anc.tight_layout()
    fig_anc.savefig('anchors_'+prob)

    plt.close(fig_anc)
        


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])





