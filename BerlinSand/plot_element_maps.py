# -*- coding: cp1252 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
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
import C_ContinuumResults

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


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

    # get indices for all nodes:
    all_nodes = mesh.get_list_of_nodes()

    # select beam elements of retaining wall:
    beams = mesh.get_list_of_elements('BEAMS',mat_filter=[3])

    # select truss elements sorted by material:
    trusses = [mesh.get_list_of_elements('TRUSSES',mat_filter=[mat]) for mat in [4,5,6]]

    # get all volume elements:
    vollist = mesh.get_list_of_elements('VOLUMICS', time=9, only_active=True)

    # get nodal and element results:
    nd_rsl = C_NodRes.NodalResults(mesh, his, rcf)
    vol_rsl = C_ContinuumResults.Continuum_EleResults(mesh, his, rcf)

    all_node_indices = [node.index for node in all_nodes]
    XY = [mesh.get_node(index).get_xyz()[:2] for index in all_node_indices]

    res = [nd_rsl.get_nodes_displacements(all_node_indices, solution_indices[-1], comp="X"),
           nd_rsl.get_nodes_displacements(all_node_indices, solution_indices[-1], comp="Y")]
    scale = 100
    XYdisp = [(xy[0]+scale*res[0][index],
               xy[1]+scale*res[1][index]) for index, xy in enumerate(XY)]

    all_vol_indices = [element.index for element in vollist]

    res = vol_rsl.get_rsl_for_sel_elements(all_vol_indices, solution_indices[-1],
                                           'STRESSES')

    patches = []
    cvect = []
    for ke, ele in enumerate(vollist):
        inel = ele.nodes
        crds = [XYdisp[index-1] for index in inel]
        patches.append(Polygon(crds))
        cvect.append(res[ke][0][0])

    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(1,1,1)

    pc = PatchCollection(patches, cmap='viridis',
                         edgecolors='none', antialiased=False)
    pc.set_array(cvect)
    ax.add_collection(pc)    
    
    cb = fig.colorbar(pc)
    cb.set_label('Horizontal stress [kPa]')

    for beam in beams:
        nodes = beam.get_nodes()
        xb = [XYdisp[nd.index-1][0] for nd in nodes]
        yb = [XYdisp[nd.index-1][1] for nd in nodes]
        ax.plot(xb,yb,color='orange',lw=2)

    for trusses0 in trusses:
        for truss in trusses0:
            nodes = truss.get_nodes()
            xb = [XYdisp[nd.index-1][0] for nd in nodes]
            yb = [XYdisp[nd.index-1][1] for nd in nodes]
            ax.plot(xb,yb,color='orange',lw=2)

    ax.axis('off')
    ax.set_aspect('equal')

    fig.tight_layout()
    fig.savefig('Sxx_map_'+prob+'_T=%1.0f'%(step_to_plot))

    plt.close(fig)
        

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])





