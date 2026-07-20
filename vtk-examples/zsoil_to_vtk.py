#-------------------------------------------------------------------------------
# Name:        zsoil_to_vtk
# Purpose:     Read a ZSoil project (mesh + full result history) via the ZSoilPy3
#              SDK and export it as a ParaView-readable VTK time series: one .vtu
#              per (converged step, element group) plus a .pvd collection. This
#              mirrors the file-naming and per-group-file conventions of
#              M. Preisig's zsoil_tools/vtktools.py (write_vtu / get_tstr).
#
# Usage:       python zsoil_to_vtk.py <project_path_without_extension> [--outdir DIR]
#                  [--steps 3,5,7,9 | --last-only] [--groups VOLUMICS,CONTACT]
#                  [--nodal-arrays DISP_TRA,PPRESS | all]
#                  [--element-arrays GROUP=NAME1,NAME2,... | GROUP=all] (repeatable)
#              Also usable programmatically via export_to_vtk(project, ...) - see
#              its docstring for the same options as real Python arguments.
#
# Requires:    numpy, vtk  (pip install numpy vtk)
#              ZSoilPy3 SDK on disk (see zsoilpy_env.py for how to point at yours)
#-------------------------------------------------------------------------------

import os
import sys
import argparse
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from zsoilpy_env import ensure_zsoilpy_on_path
ensure_zsoilpy_on_path()

from C_Mesh import *  # noqa: F401,F403  (brings in Mesh, Element, Material, ...)
from C_HistoryOfExecution import HistoryOfExecution
from C_Rcf_info import RCF_info
from C_NodalResults import NodalResults
from C_EleResults import EleResults

import numpy as np
import vtk
from vtk.util import numpy_support


# =====================================================================
# ZSoil element label -> (VTK cell type, expected node count)
# expected node count == None means "use all nodes of the element" (poly-line)
# Node ordering is assumed to follow the SDK's own local numbering, which for
# the continuum elements (Q4, T3, B8, W6, TH4) matches the standard VTK
# convention (bottom face first, then top face, both counter-clockwise).
# =====================================================================
VTK_CELL_MAP = {
    "Q4":    (vtk.VTK_QUAD,       4),
    "Q4ES":  (vtk.VTK_QUAD,       4),
    "T3":    (vtk.VTK_TRIANGLE,   3),
    "B8":    (vtk.VTK_HEXAHEDRON, 8),
    "B8ES":  (vtk.VTK_HEXAHEDRON, 8),
    "W6":    (vtk.VTK_WEDGE,      6),
    "TH4":   (vtk.VTK_TETRA,      4),
    "SXQ4":  (vtk.VTK_QUAD,       4),
    "SHQ4":  (vtk.VTK_QUAD,       4),
    "TRS2":  (vtk.VTK_LINE,       2),
    "BEL2":  (vtk.VTK_LINE,       2),
    "SBEL2": (vtk.VTK_LINE,       2),
    "M_L2":  (vtk.VTK_LINE,       2),
    "M_Q4":  (vtk.VTK_QUAD,       4),
    "M_T3":  (vtk.VTK_TRIANGLE,   3),
    # C_L2, C_Q4, C_T3, CB_L2, CNNPL are interface elements handled separately
    # (see SPLIT_LABELS below): they carry two coincident geometries (one per
    # side of the interface) that must become two separate VTK cells.
    "HEXL2": (vtk.VTK_POLY_LINE,  None),
    "BHPPE": (vtk.VTK_POLY_LINE,  None),
    "BHE1U": (vtk.VTK_POLY_LINE,  None),
    "BHE2U": (vtk.VTK_POLY_LINE,  None),
    "BHCXA": (vtk.VTK_POLY_LINE,  None),
    "BHCXC": (vtk.VTK_POLY_LINE,  None),
}

# Interface/contact element labels: each one models two coincident geometries
# (e.g. the soil side and the wall side of a diaphragm-wall interface) that can
# separate/slide relative to each other during the analysis. In the .dat file
# they are stored as ONE element whose node list is the concatenation of both
# sides' nodes: [side_A_node_1..N, side_B_node_1..N]. This was confirmed
# empirically for C_L2 (a 2D line-to-line contact) on a real project
# (tests/diaphragm-wall): a 'C_L2' element there has 4 nodes, nodes[0]/nodes[2]
# and nodes[1]/nodes[3] are coincident pairs, i.e. first half = side A (a line),
# second half = side B (a line). The same halving convention is assumed here for
# C_Q4 (2 quads), C_T3 (2 triangles), CB_L2 (2 lines) and CNNPL (2 points), by
# analogy, but is UNVERIFIED for those since no sample project with those labels
# was available to check against.
# value: (single-side VTK cell type, single-side node count)
SPLIT_LABELS = {
    "C_L2":  (vtk.VTK_LINE,     2),
    "C_Q4":  (vtk.VTK_QUAD,     4),
    "C_T3":  (vtk.VTK_TRIANGLE, 3),
    "CB_L2": (vtk.VTK_LINE,     2),
    "CNNPL": (vtk.VTK_VERTEX,   1),
}

# Element groups written out as separate .vtu files, with the short filename
# suffix used in zsoil_tools/vtktools.py (write_vtu: '_vol', '_shell', '_truss',
# '_beam', '_cnt', '_mem'). HEAT_EXCH has no precedent there; 'heat' is our own
# extension. The dict order also fixes each group's stable .pvd "part" index.
GROUP_SHORT_NAMES = {
    "VOLUMICS":  "vol",
    "SHELLS":    "shell",
    "TRUSSES":   "truss",
    "BEAMS":     "beam",
    "CONTACT":   "cnt",
    "MEMBRANE":  "mem",
    "HEAT_EXCH": "heat",
}

# Mesh only pre-computes per-element .sNN seek offsets for groups in
# Mesh.groups_with_rsl, which still excludes MEMBRANE as of this SDK version;
# ensure_seek_positions() patches that in for MEMBRANE (and any other group
# missing from that list) on first use.
# (Previously Element.group used "MEMBRANES" (plural) while the .rcf file and
# Mesh.groups_with_rsl used "MEMBRANE" (singular), which broke rcf lookups for
# this group entirely -- this has since been fixed directly in the ZSoilPy3 SDK
# (C_Element_memb_L2/Q4/T3.py, C_Element.py) to consistently use "MEMBRANE".)

# Nodal result blocks: (dict_key for NodalResults.get_nodal_results, rcf attribute holding item layout)
NODAL_RESULT_BLOCKS = [
    ("NODAL_RSL", "nod_rcf"),
    ("NODAL_VEL", "nod_rcf_v"),
    ("NODAL_ACC", "nod_rcf_a"),
]

# "Material", "EF" and "LF" aren't .rcf result items -- they're per-element
# material/exf/unlf properties written by add_element_identity_arrays() from
# the mesh itself, not read via rcf.give_ele_rcf_items_for_group(). They're
# still filtered against each group's DEFAULT_ELEMENT_ARRAYS entry like any
# other array, but add_element_result_arrays() (which only knows about real
# .rcf items) needs to ignore them when checking for "requested but missing"
# names, since it will never find them there.
IDENTITY_ARRAY_NAMES = {"Material", "EF", "LF"}

# Which named result arrays to export by default, instead of every array the
# .rcf happens to contain. `None` (as a dict value) means "export everything"
# for that group -- Material/EF/LF included, since add_element_identity_arrays()
# treats None the same way.
DEFAULT_NODAL_ARRAYS = ["DISP_TRA", "PPRESS", "PRES_HEAD"]
DEFAULT_ELEMENT_ARRAYS = {
    "VOLUMICS":  ["Material", "EF", "LF", "STRESSES", "STRAINS", "STR_LEVEL", "PLA_CODE", "SAT_EFF"],
    "SHELLS":    ["Material", "EF", "LF", "THICK", "SMFORCE", "SMOMENT", "SQFORCE"],
    "TRUSSES":   ["FORCE"],
    "BEAMS":     ["FORCE", "MOMENT"],
    "CONTACT":   ["STRESSES", "STRAINS", "PLA_CODE", "STR_LEVEL"],
    "MEMBRANE":  None,
    "HEAT_EXCH": None,
}


def get_tstr(t):
    """Format a time value the same way as zsoil_tools.vtktools.get_tstr():
    T=1.2 -> '001_20' (integer part zero-padded to 3 digits, fractional part
    as hundredths zero-padded to 2 digits)."""
    intpart = int(float("%1.2f" % t))
    frac = int(round(100 * (t - intpart)))
    if frac >= 100:
        intpart += 1
        frac -= 100
    return "{:03d}_{:02d}".format(intpart, frac)


def get_vtk_cell_type(element, warned_labels):
    n_nodes = len(element.nodes)
    if n_nodes == 0:
        return None

    entry = VTK_CELL_MAP.get(element.label)
    if entry is not None:
        vtk_type, expected_n = entry
        if expected_n is None or expected_n == n_nodes:
            return vtk_type
        if element.label not in warned_labels:
            print("  warning: element label '{}' expected {} nodes but has {}; "
                  "using a generic fallback shape".format(element.label, expected_n, n_nodes))
            warned_labels.add(element.label)

    fallback_by_count = {1: vtk.VTK_VERTEX, 2: vtk.VTK_LINE, 3: vtk.VTK_TRIANGLE, 4: vtk.VTK_QUAD}
    if element.label not in VTK_CELL_MAP:
        if element.label not in warned_labels:
            print("  warning: unrecognized element label '{}' ({} nodes); "
                  "using a generic fallback shape".format(element.label, n_nodes))
            warned_labels.add(element.label)
        if n_nodes in fallback_by_count:
            return fallback_by_count[n_nodes]
    if n_nodes >= 2:
        return vtk.VTK_POLY_LINE
    return vtk.VTK_POLY_VERTEX


def expand_element_to_cells(element, warned_labels):
    """Return a list of (vtk_cell_type, node_id_list, side_label) for one
    element. Normally a single entry with side_label=None. For interface
    elements (SPLIT_LABELS) whose node list holds both sides, this returns two
    entries (side_label 'A'/'B'), one per physical geometry, so the two sides
    can move independently in the VTK output instead of being drawn as one
    self-crossing cell."""
    split_entry = SPLIT_LABELS.get(element.label)
    if split_entry is not None:
        vtk_type, side_n = split_entry
        n_nodes = len(element.nodes)
        if n_nodes == 2 * side_n:
            return [
                (vtk_type, element.nodes[:side_n], "A"),
                (vtk_type, element.nodes[side_n:], "B"),
            ]
        elif n_nodes == side_n:
            return [(vtk_type, element.nodes, None)]
        else:
            if element.label not in warned_labels:
                print("  warning: interface element label '{}' expected {} or {} nodes but has {}; "
                      "using a generic fallback shape".format(element.label, side_n, 2 * side_n, n_nodes))
                warned_labels.add(element.label)
            vtk_type = get_vtk_cell_type(element, warned_labels)
            if vtk_type is None:
                return []
            return [(vtk_type, element.nodes, None)]

    vtk_type = get_vtk_cell_type(element, warned_labels)
    if vtk_type is None:
        return []
    return [(vtk_type, element.nodes, None)]


def ensure_seek_positions(mesh, rcf, ele_group):
    """Mesh only precomputes .sNN seek offsets for groups in mesh.groups_with_rsl.
    Patch in the same computation for any other group (e.g. MEMBRANE) on first use."""
    if ele_group in mesh.groups_with_rsl:
        return
    flag_name = "_seek_done_" + ele_group
    if getattr(mesh, flag_name, False):
        return
    gp_size = rcf.give_one_gp_size(ele_group)
    if gp_size <= 0:
        raise RuntimeError("no .rcf entry found for group '{}'".format(ele_group))
    last_pos = 0
    for e in mesh.elements:
        if e.group == ele_group:
            e.seek_in_rsl = last_pos
            last_pos += len(e.xsiGP) * gp_size * 4
    setattr(mesh, flag_name, True)


def is_element_active(element, time):
    exf = element.give_exf()
    if exf is None:
        return True
    return exf.is_ON(time)


def build_mesh_index(mesh):
    points = vtk.vtkPoints()
    node_id_of = {}
    for i, node in enumerate(mesh.nodes):
        node_id_of[node.index] = i
        xyz = node.xyz
        points.InsertNextPoint(xyz[0], xyz[1], xyz[2] if len(xyz) > 2 else 0.0)

    group_elements = defaultdict(list)
    for e in mesh.elements:
        if e.group in GROUP_SHORT_NAMES:
            group_elements[e.group].append(e)

    return points, node_id_of, group_elements


def build_group_grid_for_step(points, node_id_of, elements_all, time, warned_labels):
    """Build a fresh grid containing only the elements of one group that are
    active at the given time. Interface elements (SPLIT_LABELS) become two
    cells, one per side. Returns (grid, cell_records, active_elements):
      - cell_records[i] = (element, side_label) for cell i ('A'/'B'/None)
      - active_elements = unique active elements (for result lookups)."""
    active_elements = [e for e in elements_all if is_element_active(e, time)]

    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)

    cell_records = []
    for e in active_elements:
        for vtk_type, node_list, side in expand_element_to_cells(e, warned_labels):
            id_list = vtk.vtkIdList()
            for n in node_list:
                id_list.InsertNextId(node_id_of[n])
            grid.InsertNextCell(vtk_type, id_list)
            cell_records.append((e, side))

    return grid, cell_records, active_elements


def add_element_identity_arrays(grid, mesh, cell_records, wanted_names=None):
    """
    Always writes ElementIndex/ElementLabel/Side (needed to identify cells).
    Material/EF/LF are written too unless `wanted_names` is given and
    excludes them (see IDENTITY_ARRAY_NAMES).
    """
    n = len(cell_records)
    ele_index = np.empty(n, dtype=np.int32)
    material = np.empty(n, dtype=np.int32)
    ef = np.empty(n, dtype=np.int32)
    lf = np.empty(n, dtype=np.int32)

    label_arr = vtk.vtkStringArray()
    label_arr.SetName("ElementLabel")
    label_arr.SetNumberOfTuples(n)
    side_arr = vtk.vtkStringArray()
    side_arr.SetName("Side")
    side_arr.SetNumberOfTuples(n)

    for i, (e, side) in enumerate(cell_records):
        mat = mesh.materials[e.material_index - 1]
        ele_index[i] = e.index
        material[i] = mat.index_in_inp
        ef[i] = mat.exf_index
        lf[i] = mat.unlf_index
        label_arr.SetValue(i, e.label)
        side_arr.SetValue(i, side or "")

    cdata = grid.GetCellData()
    ele_index_arr = numpy_support.numpy_to_vtk(ele_index, deep=True)
    ele_index_arr.SetName("ElementIndex")
    cdata.AddArray(ele_index_arr)

    for name, arr in [("Material", material), ("EF", ef), ("LF", lf)]:
        if wanted_names is not None and name not in wanted_names:
            continue
        vtk_arr = numpy_support.numpy_to_vtk(arr, deep=True)
        vtk_arr.SetName(name)
        cdata.AddArray(vtk_arr)
    cdata.AddArray(label_arr)
    cdata.AddArray(side_arr)


def add_nodal_result_arrays(grid, mesh, rcf, nodal_results, sol_index, wanted_names=None, warned_labels=None):
    """
    `wanted_names`, if given, restricts output to nodal result items whose
    name matches (e.g. "DISP_TRA"); a wanted item can live in any of the
    NODAL_RESULT_BLOCKS. None exports every item found.
    """
    remaining = set(wanted_names) if wanted_names is not None else None
    for dict_key, rcf_attr in NODAL_RESULT_BLOCKS:
        rcf_items = getattr(rcf, rcf_attr)
        if not rcf_items:
            continue
        try:
            rows = nodal_results.get_nodal_results(mesh.nodes, [sol_index], dict_key)
        except (IOError, OSError) as ex:
            print("  warning: could not read {} results: {}".format(dict_key, ex))
            continue
        if not rows:
            continue
        data = np.array(rows, dtype=np.float32)  # columns: [node_index, v0, v1, ...]
        offset = 0
        for name, n_comp, comps in rcf_items:
            stripped = name.strip()
            if wanted_names is not None and stripped not in wanted_names:
                offset += n_comp
                continue
            if remaining is not None:
                remaining.discard(stripped)
            block = np.ascontiguousarray(data[:, 1 + offset: 1 + offset + n_comp])
            arr = numpy_support.numpy_to_vtk(block, deep=True)
            arr.SetName(stripped)
            grid.GetPointData().AddArray(arr)
            offset += n_comp

    for name in (remaining or ()):
        key = "missing-nodal-array:{}".format(name)
        if warned_labels is not None and key in warned_labels:
            continue
        print("  warning: requested nodal array '{}' not found".format(name))
        if warned_labels is not None:
            warned_labels.add(key)


def add_element_result_arrays(grid, rcf, ele_results, group, cell_records, active_elements, sol_index,
                               wanted_names=None, warned_labels=None):
    """
    `wanted_names`, if given, restricts output to element result items of
    `group` whose name matches (e.g. "STRESSES" for VOLUMICS). None exports
    every item available for that group.
    """
    if not active_elements:
        return
    items = rcf.give_ele_rcf_items_for_group(group)
    if not items:
        return

    try:
        ensure_seek_positions(ele_results.mesh, rcf, group)
        rows = ele_results.get_element_results_ex(active_elements, [sol_index])
        gp_size = rcf.give_one_gp_size(group)
    except (IOError, OSError, RuntimeError) as ex:
        print("  warning: could not read element results for group '{}': {}".format(group, ex))
        return
    if not rows or gp_size <= 0:
        return

    # average all Gauss points of each (unique) element into a single value.
    # Interface elements split into two cells (see SPLIT_LABELS) share the
    # same underlying Gauss-point result -- it describes the interface as a
    # whole, not one side -- so both of their cells get the same row here.
    sums = {}
    counts = {}
    for row in rows:
        ele_idx = row[0]
        vals = np.asarray(row[2:2 + gp_size], dtype=np.float64)
        if ele_idx not in sums:
            sums[ele_idx] = np.zeros(gp_size)
            counts[ele_idx] = 0
        sums[ele_idx] += vals
        counts[ele_idx] += 1

    avg_by_ele_index = {idx: sums[idx] / counts[idx] for idx in sums}

    avg_matrix = np.zeros((len(cell_records), gp_size), dtype=np.float64)
    for i, (e, side) in enumerate(cell_records):
        if e.index in avg_by_ele_index:
            avg_matrix[i] = avg_by_ele_index[e.index]

    cdata = grid.GetCellData()
    offset = 0
    found = set()
    for name, n_comp, comps in items:
        stripped = name.strip()
        if wanted_names is not None and stripped not in wanted_names:
            offset += n_comp
            continue
        found.add(stripped)
        block = np.ascontiguousarray(avg_matrix[:, offset:offset + n_comp]).astype(np.float32)
        arr = numpy_support.numpy_to_vtk(block, deep=True)
        arr.SetName(stripped)
        cdata.AddArray(arr)
        offset += n_comp

    # Material/EF/LF are handled by add_element_identity_arrays(), not here --
    # exclude them so they don't get reported as missing .rcf items.
    missing = (set(wanted_names) - found - IDENTITY_ARRAY_NAMES) if wanted_names is not None else ()
    for name in missing:
        key = "missing-element-array:{}:{}".format(group, name)
        if warned_labels is not None and key in warned_labels:
            continue
        print("  warning: requested element array '{}' not found for group '{}'".format(name, group))
        if warned_labels is not None:
            warned_labels.add(key)


def write_pvd(pvd_path, entries):
    lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
        '  <Collection>',
    ]
    for fname, t, part in entries:
        lines.append('    <DataSet timestep="{:.10g}" group="" part="{}" file="{}"/>'.format(t, part, fname))
    lines.append('  </Collection>')
    lines.append('</VTKFile>')
    with open(pvd_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def export_to_vtk(project, outdir=None, steps=None, last_only=False, groups=None,
                   nodal_arrays=DEFAULT_NODAL_ARRAYS, element_arrays=DEFAULT_ELEMENT_ARRAYS):
    """
    Export one project to a VTK/ParaView time series, one .vtu per (step,
    element group) plus a .pvd collection.

    steps:          list of step times to export (exact match against the
                     project's own time values, e.g. [3, 5, 7, 9]), or None
                     for every converged step. Ignored if `last_only` is True.
    last_only:       export only the last converged step.
    groups:          list of element groups to export (see GROUP_SHORT_NAMES
                     for the choices), or None for every group present in
                     the mesh.
    nodal_arrays:    list of nodal result names to export (e.g. "DISP_TRA"),
                     or None to export every nodal array available.
    element_arrays:  {group: list_of_names} restricting element result
                     arrays per group; a group missing from this dict (or
                     mapped to None) exports every array available for it.
    """
    if not os.path.exists(project + ".dat"):
        sys.exit("Could not find '{}.dat'".format(project))

    outdir = outdir or os.path.join(os.path.dirname(os.path.abspath(project)), "vtk_output")
    os.makedirs(outdir, exist_ok=True)
    basename = os.path.basename(project)

    if groups is not None:
        unknown = set(groups) - set(GROUP_SHORT_NAMES)
        if unknown:
            sys.exit("Unknown element group(s): {}. Choices: {}".format(
                ", ".join(sorted(unknown)), ", ".join(GROUP_SHORT_NAMES)))
    group_names = groups if groups is not None else list(GROUP_SHORT_NAMES)
    element_arrays = element_arrays or {}

    print("Reading history of execution...")
    his = HistoryOfExecution(project)
    print("Reading mesh...")
    mesh = Mesh(project)
    print("Reading result configuration (.rcf)...")
    rcf = RCF_info(project)

    solution_indices = his.give_converged_time_solutions()
    if not solution_indices:
        print("No converged time steps found; falling back to all stored steps.")
        solution_indices = his.give_solutions_with_plot_status()
    if last_only:
        solution_indices = solution_indices[-1:]
    elif steps is not None:
        wanted_steps = set(steps)
        solution_indices = [si for si in solution_indices if his.data[si][his.DATA_TIME] in wanted_steps]
        if not solution_indices:
            sys.exit("None of the requested --steps {} matched a converged step in this project.".format(
                sorted(wanted_steps)))
    if not solution_indices:
        sys.exit("No stored solution steps found in this project.")

    print("Indexing mesh ({} nodes, {} elements)...".format(len(mesh.nodes), len(mesh.elements)))
    points, node_id_of, group_elements = build_mesh_index(mesh)
    part_index = {group: i for i, group in enumerate(GROUP_SHORT_NAMES)}
    warned_labels = set()

    nodal_results = NodalResults(mesh, his, rcf)
    ele_results = EleResults(mesh, his, rcf)

    entries = []
    for step_num, sol_index in enumerate(solution_indices):
        time = his.data[sol_index][his.DATA_TIME]
        tstr = get_tstr(time)
        print("Exporting step {}/{} (solution index {}, time={:.6g})...".format(
            step_num + 1, len(solution_indices), sol_index, time))

        for group in group_names:
            short_name = GROUP_SHORT_NAMES[group]
            elements_all = group_elements.get(group, [])
            if not elements_all:
                continue

            grid, cell_records, active_elements = build_group_grid_for_step(
                points, node_id_of, elements_all, time, warned_labels)
            if not cell_records:
                continue

            add_element_identity_arrays(grid, mesh, cell_records, wanted_names=element_arrays.get(group))
            add_nodal_result_arrays(grid, mesh, rcf, nodal_results, sol_index,
                                     wanted_names=nodal_arrays, warned_labels=warned_labels)
            add_element_result_arrays(grid, rcf, ele_results, group, cell_records, active_elements, sol_index,
                                       wanted_names=element_arrays.get(group), warned_labels=warned_labels)

            fname = "{}_{}_{}.vtu".format(basename, tstr, short_name)
            out_path = os.path.join(outdir, fname)
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetDataModeToBinary()
            writer.SetFileName(out_path)
            writer.SetInputData(grid)
            writer.Write()

            entries.append((fname, time, part_index[group]))

    pvd_path = os.path.join(outdir, basename + ".pvd")
    write_pvd(pvd_path, entries)
    print("Done. Wrote {} file(s) across {} step(s).".format(len(entries), len(solution_indices)))
    print("Open in ParaView: {}".format(pvd_path))


def _parse_csv(value):
    return [v.strip() for v in value.split(",") if v.strip()]


def main():
    parser = argparse.ArgumentParser(
        description="Export a ZSoil project (mesh + full result history) to a VTK/ParaView time series, "
                     "one .vtu per (step, element group)."
    )
    parser.add_argument("project", help="ZSoil project path without file extension, e.g. .../my_model")
    parser.add_argument("--outdir", default=None, help="Output directory (default: <project_dir>/vtk_output)")

    step_group = parser.add_mutually_exclusive_group()
    step_group.add_argument("--steps", default=None,
                             help="Comma-separated step times to export, e.g. '3,5,7,9' "
                                  "(default: every converged step)")
    step_group.add_argument("--last-only", action="store_true",
                             help="Export only the last converged step")

    parser.add_argument("--groups", default=None,
                         help="Comma-separated element groups to export (default: all present). "
                              "Choices: {}".format(", ".join(GROUP_SHORT_NAMES)))
    parser.add_argument("--nodal-arrays", default=None,
                         help="Comma-separated nodal result names to export, or 'all'. "
                              "Default: {}".format(",".join(DEFAULT_NODAL_ARRAYS)))
    parser.add_argument("--element-arrays", default=[], action="append", metavar="GROUP=NAME1,NAME2,...",
                         help="Restrict one group's element result arrays; NAMEs may be 'all'. Repeatable. "
                              "Groups not given here keep DEFAULT_ELEMENT_ARRAYS' own default for that "
                              "group (every group defaults to at least Material,EF,LF; "
                              "VOLUMICS: {}; SHELLS adds THICK).".format(
                                  ",".join(DEFAULT_ELEMENT_ARRAYS["VOLUMICS"])))
    args = parser.parse_args()

    steps = [float(s) for s in _parse_csv(args.steps)] if args.steps else None

    groups = _parse_csv(args.groups) if args.groups else None

    nodal_arrays = DEFAULT_NODAL_ARRAYS
    if args.nodal_arrays:
        nodal_arrays = None if args.nodal_arrays.strip().lower() == "all" else _parse_csv(args.nodal_arrays)

    element_arrays = dict(DEFAULT_ELEMENT_ARRAYS)
    for spec in args.element_arrays:
        if "=" not in spec:
            sys.exit("--element-arrays expects GROUP=NAME1,NAME2,...; got '{}'".format(spec))
        group, names = spec.split("=", 1)
        group = group.strip()
        if group not in GROUP_SHORT_NAMES:
            sys.exit("Unknown element group '{}'. Choices: {}".format(group, ", ".join(GROUP_SHORT_NAMES)))
        element_arrays[group] = None if names.strip().lower() == "all" else _parse_csv(names)

    export_to_vtk(args.project, outdir=args.outdir, steps=steps, last_only=args.last_only,
                  groups=groups, nodal_arrays=nodal_arrays, element_arrays=element_arrays)


if __name__ == "__main__":
    main()
