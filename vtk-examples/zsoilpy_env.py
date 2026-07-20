"""
Locate the ZSoilPy3 package shipped with a ZSoil installation and add it to
sys.path so the C_* result-reading modules (C_Mesh, C_HistoryOfExecution,
C_NodalResults, ...) can be imported.

Every example script in this folder calls ensure_zsoilpy_on_path() before
importing anything from ZSoilPy3. To point these examples at your own ZSoil
installation, either:

  * set the ZSOILPY3_PATH environment variable to your ZSoilPy3 folder, e.g.
        setx ZSOILPY3_PATH "C:\\Program Files\\ZSoil\\Tools v2026\\ZSoilPy3"
  * or edit DEFAULT_ZSOILPY3_PATH below.
"""
import os
import sys

DEFAULT_ZSOILPY3_PATH = r"C:\Program Files\ZSoil\Tools v2026\ZSoilPy3"


def ensure_zsoilpy_on_path():
    path = os.path.normpath(os.environ.get("ZSOILPY3_PATH", DEFAULT_ZSOILPY3_PATH))
    if not os.path.isdir(path):
        raise FileNotFoundError(
            "ZSoilPy3 not found at '%s'. Set the ZSOILPY3_PATH environment "
            "variable, or edit DEFAULT_ZSOILPY3_PATH in zsoilpy_env.py, to "
            "point at your ZSoil installation's Tools/ZSoilPy3 folder." % path
        )
    if path not in sys.path:
        sys.path.append(path)
