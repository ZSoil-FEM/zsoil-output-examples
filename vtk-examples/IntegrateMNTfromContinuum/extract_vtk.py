import os
import sys

# vtk-examples/, where zsoil_to_vtk.py and zsoilpy_env.py live:
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

from zsoil_to_vtk import export_to_vtk

export_to_vtk(
    "zsoil-files/3Dtunnel",
    outdir="pv",
    steps=[2],                    # default: every converged step
    groups=["VOLUMICS"],        # default: every group present
    # nodal_arrays=["DISP_TRA", "PPRESS"],   # default: DISP_TRA, PPRESS, PRES_HEAD
    # element_arrays={"VOLUMICS": ["STRESSES", "PLA_CODE"]},  # default: see zsoil_to_vtk.DEFAULT_ELEMENT_ARRAYS
)
export_to_vtk(
    "zsoil-files/3Dtunnel",
    outdir="pv",
    steps=[2],                    # default: every converged step
    groups=["SHELLS"],        # default: every group present
    include_inactive=True,
    # nodal_arrays=["DISP_TRA", "PPRESS"],   # default: DISP_TRA, PPRESS, PRES_HEAD
    # element_arrays={"VOLUMICS": ["STRESSES", "PLA_CODE"]},  # default: see zsoil_to_vtk.DEFAULT_ELEMENT_ARRAYS
)
