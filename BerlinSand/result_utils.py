"""
Helper shared by the Berlin Sand plotting scripts for selecting which solver
steps to read results from.
"""


def select_static_steps(his, steps=None):
    """
    Return the solution indices for converged, non-stability-driver steps
    (SF == 0) recorded in a HistoryOfExecution.

    If `steps` is given, only the entries whose time matches one of the
    values in `steps` are kept; otherwise every converged static step is
    returned.

    Returns (selected_indices, all_static_indices) - all_static_indices
    lists every converged static step regardless of `steps`, which is handy
    when a caller wants both a subset to plot and the full run for context.
    """
    selected = []
    all_static = []
    for sol_ind in his.give_converged_solutions():
        time = his.data[sol_ind][his.DATA_TIME]
        sf = his.data[sol_ind][his.DATA_SF]
        if sf == 0:
            all_static.append(sol_ind)
            if steps is None or time in steps:
                selected.append(sol_ind)
    return selected, all_static
