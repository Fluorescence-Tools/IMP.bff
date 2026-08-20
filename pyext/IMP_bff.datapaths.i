/*
 * Where the shipped cgdye data lives, and how to reach it.
 *
 * Path arithmetic over `IMP.bff.get_data_path`, which is a generated Python
 * function of this module -- so this is one of the few things that is more
 * natural in `%pythoncode` than in C++: the surface is `pathlib.Path` and the
 * parts are variadic, and a std::string port would be a downgrade at both ends.
 *
 * It was `cgdye/utils.py`, then `tools.py`, and the module it lived in never
 * had anything else in it.
 */

%pythoncode %{
import pathlib as _pathlib


def _cgdye_data_root():
    """The cgdye data directory.

    Templates, input structures and restraint files are IMP module *data*, not
    package sources: they live in ``imp.bff/data/cgdye`` and are reached through
    ``IMP.bff.get_data_path``. Deriving them from ``__file__`` instead would tie
    them to where the package happens to sit, which is exactly what broke when
    cgdye moved out of imp-tricks.
    """
    return _pathlib.Path(get_data_path("cgdye"))


def _join_parts(parts):
    """Join path parts, handling nested path strings."""
    result = _pathlib.Path()
    for part in parts:
        result = result / _pathlib.Path(part)
    return result


def get_template_dir(*parts):
    """The template directory, optionally with subpath parts."""
    return _cgdye_data_root() / "templates" / _join_parts(parts)


def get_structure_dir(*parts):
    """The input-structure directory, optionally with subpath parts."""
    return _cgdye_data_root() / "inputs" / "structures" / _join_parts(parts)


def get_output_dir(*parts):
    """The output directory, optionally with subpath parts.

    Relative to the working directory, not to the package: module data is
    read-only and installed, so nothing may be written beside it.
    """
    return _pathlib.Path.cwd() / "output" / _join_parts(parts)


def ensure_dir(path):
    """Ensure a directory exists and return it."""
    path = _pathlib.Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path
%}
