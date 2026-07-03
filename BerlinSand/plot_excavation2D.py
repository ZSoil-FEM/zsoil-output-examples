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

    # select beam elements of retaining wall:
    blist = mesh.get_list_of_elements('BEAMS',mat_filter=[3])
    beam_indices = [beam.index for beam in blist]

    # get indices for all nodes in model:
    node_indices = [node.index for node in mesh.nodes]

    fig = plt.figure(figsize=(12,8))
    AX = [fig.add_subplot(1,3,k+1) for k in range(3)]
    handles = [[] for _ in range(3)]
    
    ylim = [1e10,-1e10]
    Mlim = [1e10,-1e10]
    Tlim = 0
    ulim = 0
    kkt = -1
    for sol in solution_indices:
        tx = his.data[sol][his.DATA_TIME]
        if abs(tx%1)<1e-3 and tx>2:
            kkt += 1
            # select MZ for current time step:
            moment = beam_rsl.get_moment_for_sel_elements(beam_indices, sol, "Z")
            # select TY for current time step:
            shear = beam_rsl.get_force_for_sel_elements(beam_indices, sol, "Y")
            # select ux for current time step:
            disp = nd_rsl.get_nodes_displacements(node_indices, sol, comp="X")

            for kk in range(3):
                handles[kk].append(0)

            for ke,beam_index in enumerate(beam_indices):
                beam = blist[ke]
                inel = [nd.index for nd in beam.get_nodes()]
                crds = beam.get_ele_coord()
                ylim[0] = min(ylim[0],min(crds[:,1]))
                ylim[1] = max(ylim[1],max(crds[:,1]))
                n = np.array([-crds[1,1]+crds[0,1],crds[1,0]-crds[0,0]])
                l = np.linalg.norm(n)
                n /= l
                v = n[0]*np.array([moment[ke]-0.5*l*shear[ke],
                                   moment[ke]+0.5*l*shear[ke]])
                Mlim[0] = min(Mlim[0],min(v))
                Mlim[1] = max(Mlim[1],max(v))
                handles[0][-1] = AX[0].plot(v,crds[:,1],color=colors[kkt],label='T=%1.0f: %s'%(tx,step_labels[tx]))[0]

                v = n[0]*np.array([shear[ke],
                                   shear[ke]])
                Tlim = max([Tlim,max(v,key=abs)],key=abs)
                handles[1][-1] = AX[1].plot(v,crds[:,1],color=colors[kkt],label='T=%1.0f: %s'%(tx,step_labels[tx]))[0]

                v = 1e3*np.array([disp[inel[0]-1],
                                  disp[inel[1]-1]])
                ulim = max([ulim,max(v,key=abs)],key=abs)
                handles[2][-1] = AX[2].plot(v,crds[:,1],color=colors[kkt],label='T=%1.0f: %s'%(tx,step_labels[tx]))[0]


    tlist = mesh.get_list_of_elements('TRUSSES',mat_filter=[4,5,6])
    for truss in tlist:
        crds = truss.get_ele_coord()
        if crds[0][0]<30.001:
            AX[0].annotate('',xy=(0,crds[0][1]),xytext=(-20,0),
                           textcoords='offset points',ha='right',va='center',
                           arrowprops=dict(facecolor='black',headwidth=6,width=2))
                

    for kax,ax in enumerate(AX):
        ax.set_ylabel(u'Depth [m]')
        ax.grid('on')
        ax.set_ylim(ylim)
        ax.legend(handles=handles[kax],labels=[v.get_label() for v in handles[kax]],loc='lower right')

    AX[0].set_xlabel('Bending moment [kNm/m]')
    AX[1].set_xlabel('Shear force [kN/m]')
    AX[2].set_xlabel('Horizontal displacement [mm]')

    AX[0].annotate('$M_{min}=%1.0f$ kNm/m\n$M_{max}=%1.0f$ kNm/m'%(Mlim[0],Mlim[1]),
                   xy=(0.05,0.02),xycoords='axes fraction',rotation=90,va='bottom')
    AX[1].annotate('$|T|_{max}=%1.0f$ kN/m'%(Tlim),
                   xy=(0.05,0.02),xycoords='axes fraction',rotation=90,va='bottom')
    AX[2].annotate('$u_{h,min}=%1.0f$ mm'%(ulim),
                   xy=(0.05,0.02),xycoords='axes fraction',rotation=90,va='bottom')

    fig.tight_layout()
    fig.savefig('wall_'+prob)

    plt.close(fig)
        

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])





