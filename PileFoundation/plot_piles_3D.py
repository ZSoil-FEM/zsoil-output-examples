#-------------------------------------------------------------------------------
# Name:        plot_piles_3D
# Purpose:     extract and plot normal force diagrams along piles, plot diagrams
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


# Projection of 3D view onto 2D:
azim = 0.3#ax.azim*pi/180
elev = 0.5#ax.elev*pi/180
elev *= 1.2          # this seems to improve the outcome

a_vec = np.array([np.cos(azim),np.sin(azim),0])
normal = np.cos(elev)*a_vec + np.array([0,0,np.sin(elev)])
z_vec = np.array([0,1,0])
y_comp = z_vec - (z_vec@normal)*normal
y_comp = y_comp/np.linalg.norm(y_comp)
x_comp = np.cross(y_comp,normal)

def project_2d(point0):
    point = point0*np.array([1,1,-1])
    proj_mat = np.vstack([x_comp,y_comp]) # build projection matrix
    points_2D = point @ proj_mat.T

    return points_2D

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
            if time == step_to_plot:
                solution_indices.append(sol_ind)

    # get nodal, beam and interface results:
    nd_rsl = C_NodRes.NodalResults(mesh, his, rcf)
    beam_rsl = C_Beam_Res.Beam_EleResults(mesh, his, rcf)
    cnt_rsl = C_Cnt_Res.Contact_3D_Results(mesh, his, rcf)

    # get indices for all nodes:
    all_nodes = mesh.get_list_of_nodes()

    # select beam elements:
    beams = mesh.get_list_of_elements('BEAMS',mat_filter=[9])

    # select interface elements of pile shafts:
    shaftlist = mesh.get_list_of_elements('CONTACT',mat_filter=[8])

    for sol in solution_indices:
        time = his.data[sol][his.DATA_TIME]

        pla_code = cnt_rsl.get_rsl_vec_for_sel_elements([ele.index for ele in shaftlist], sol, 'PLA_CODE')
        PC = [[v[0],v[1]] for v in pla_code]
        forces = beam_rsl.get_forces_vec_for_sel_elements([ele.index for ele in beams], sol)
        Fn = [v[0] for v in forces]


        scF = 5e-4

        fig = plt.figure(figsize=(12,10))
        fig.text(0.01,0.01,prob,size=8)
        ax = fig.add_subplot(111)

        for kke,beam in enumerate(beams):
            p0 = all_nodes[beam.nodes[0]-1].get_xyz()
            p1 = all_nodes[beam.nodes[1]-1].get_xyz()
            pF = np.array([Fn[kke]*scF,0,0])
            p0_2d = project_2d(p0)
            p1_2d = project_2d(p1)
            p0F_2d = project_2d(p0+pF)
            p1F_2d = project_2d(p1+pF)
            ax.plot([p0_2d[0],p1_2d[0]],[p0_2d[1],p1_2d[1]],'k',lw=1)
            if min([PC[kke][kgp] for kgp in [0,1]])==16:
                ax.add_patch(Polygon([p0_2d,p1_2d,p1F_2d,p0F_2d],
                                     fc='r',ec=None,alpha=0.6))
            else:
                ax.add_patch(Polygon([p0_2d,p1_2d,p1F_2d,p0F_2d],
                                     fc='k',ec=None,alpha=0.6))


        ax.annotate('T = %1.1f'%(time),
                    xy=(0.05,0.15),
                    xycoords='axes fraction',
                    size=14)
        for ka,ax in enumerate([ax]):
            ax.axis('off')

        fig.tight_layout()

        fig.savefig(prob+'_FN_T%1.0f_%1.0f'%(int(time),(time%1)*10))
        plt.close(fig)



if __name__ == "__main__":
    for t in np.linspace(3,7,9):
        main('res-files','pile_raft',t)
##    main(sys.argv[1], sys.argv[2])


