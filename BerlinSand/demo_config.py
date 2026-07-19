"""
Config for the Berlin Sand demo problem, shared by run_all.py and by the
individual plot_*.py scripts' command-line entry points. Edit this file to
point the demo at a different problem or set of construction stages - every
script that reads it will pick up the change.
"""

PROBLEM = r'HS-Brick-Exc-Berlin-Sand-2phase_2026'

# Construction stages to plot, and a human-readable label for each.
STEPS_TO_PLOT = [3, 5, 7, 9]
STEP_LABELS = {
    3: 'Exc. at -4.8 m',
    4: 'Anchor level 1',
    5: 'Exc. at -9.3 m',
    6: 'Anchor level 2',
    7: 'Exc. at -14.35 m',
    8: 'Anchor level 3',
    9: 'Exc. at -16.8 m',
}
