import os
import shutil
import subprocess
import numpy as np


import plot_excavation2D
import plot_interfaces
import plot_time_histories
import plot_nodal_maps
import plot_element_maps


# Config section

prob = r'HS-Brick-Exc-Berlin-Sand-2phase_2026'

ZSoil_exe = r'C:/Program Files/ZSoil/ZSoil 2026 v26.03/Z_Soil.exe'
working_dir = os.path.realpath('./res-files')

# End of config section

os.mkdir(working_dir)
shutil.copy('./inp-files/%s.inp'%(prob),'%s/%s.inp'%(working_dir,prob))

subprocess.check_call(
        '%s "%s\\%s" /S' % (ZSoil_exe, working_dir, prob+'.inp'),
        cwd=working_dir,
    )


step_labels = {
    3:'Exc. at -4.8 m',
    4:'Anchor level 1',
    5:'Exc. at -9.3 m',
    6:'Anchor level 2',
    7:'Exc. at -14.35 m',
    8:'Anchor level 3',
    9:'Exc. at -16.8 m'
    }
steps_to_plot = [3,5,7,9]

plot_excavation2D.main(working_dir, prob, steps_to_plot, step_labels)
plot_interfaces.main(working_dir, prob, steps_to_plot, step_labels)
plot_time_histories.main(working_dir, prob)
plot_nodal_maps.main(working_dir, prob, steps_to_plot[-1])
plot_element_maps.main(working_dir, prob, steps_to_plot[-1])

shutil.rmtree(working_dir)
