#-------------------------------------------------------------------------------
# Name:        run_all
# Purpose:     drive a ZSoil analysis end to end: run the solver on the
#              Berlin Sand excavation benchmark, then call each plot_*.py
#              module to turn the result files into figures
#
# Author:      Matthias Preisig
# Created:     2026
# Copyright:   (c) GeoDev Sarl
#-------------------------------------------------------------------------------
import os
import shutil
import subprocess

import plot_excavation2D
import plot_interfaces
import plot_time_histories
import plot_nodal_maps
import plot_element_maps

from demo_config import PROBLEM, STEPS_TO_PLOT, STEP_LABELS


# ------------------------- Config section -------------------------
# Edit these two for your own machine/problem. STEPS_TO_PLOT / STEP_LABELS
# come from demo_config.py, which the plot_*.py scripts also read as their
# own command-line defaults.
ZSoil_exe = r'C:/Program Files/ZSoil/ZSoil 2026 v26.03/Z_Soil.exe'
working_dir = os.path.realpath('./res-files')
# ----------------------- End of config section ----------------------

if not os.path.isfile(ZSoil_exe):
    raise FileNotFoundError(
        "ZSoil executable not found at '%s' - edit ZSoil_exe at the top of "
        "run_all.py to point at your installation." % ZSoil_exe
    )

os.makedirs(working_dir, exist_ok=True)
shutil.copy('./inp-files/%s.inp' % PROBLEM, '%s/%s.inp' % (working_dir, PROBLEM))

try:
    subprocess.check_call(
        '%s "%s\\%s" /S' % (ZSoil_exe, working_dir, PROBLEM + '.inp'),
        cwd=working_dir,
    )

    plot_excavation2D.main(working_dir, PROBLEM, STEPS_TO_PLOT, STEP_LABELS)
    plot_interfaces.main(working_dir, PROBLEM, STEPS_TO_PLOT, STEP_LABELS)
    plot_time_histories.main(working_dir, PROBLEM)
    plot_nodal_maps.main(working_dir, PROBLEM, STEPS_TO_PLOT[-1])
    plot_element_maps.main(working_dir, PROBLEM, STEPS_TO_PLOT[-1])
except Exception:
    print("run_all.py: an error occurred, leaving '%s' in place for inspection." % working_dir)
    raise
else:
    shutil.rmtree(working_dir)
