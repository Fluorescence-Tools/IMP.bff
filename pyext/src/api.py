"""The public surface of ``IMP.bff``: domain-scoped, with a flat view over it.

The code is organised by **domain** -- ``dye``, ``label``, ``representation``,
``scoring``, ``observables``,
``io``, ``restraints`` -- and that organisation is what a reader should learn.
A domain leaves this map when its last Python module does: the names it carried
are then attributes of ``IMP.bff`` itself, resolved by SWIG rather than lazily
here. ``observables`` was the first to go that way.
A domain is a *name*, not necessarily a directory: most are one flat module,
and the handful that are still packages earned it (see ``pyext/src/README.md``). ``IMP.bff.representation.rotamer.RotamerEnsemble``
says where a name lives and what it is about; ``IMP.bff.RotamerEnsemble`` is a
convenience over it.

So the domains are the authored thing here (:data:`BY_DOMAIN`) and the flat map
(:data:`EXPORTS`) is **derived** from them. That inversion is the point of
PRD-113 stage 7. Before it, the flat map was authored and the grouping existed
only as comments, which meant nothing could check that a name was filed
correctly -- and the test that stood in for it matched export *names* against a
regex of naming families, a list that had to grow with every feature and could
only ever assert that a name looked plausible.

What is checked now instead, in ``test/test_public_api_names.py``:

* every export resolves and carries a docstring;
* every export's module actually lies inside the domain it is filed under, so a
  misfiled name is an error rather than a comment that drifted;
* no name appears in two domains -- the flat view would silently keep one.

Names resolve **lazily** through the ``__getattr__`` that ``pyext/swig.i-in``
installs on the package, so ``import IMP.bff`` imports almost nothing and
nothing here needs ``click``. Retired names are removed, never aliased.

``IMP.bff.cgdye`` is **not** a domain and has no names here. It is explicit
all-atom dye modelling under a force field -- molecular mechanics, which is
IMP's own territory rather than a spectroscopy library's -- and as of
2026-08-18 it sits outside the layout, reachable only by module path. Its 38
flat names were removed. What was *not* molecular mechanics came out of it
first: the rotamer library into ``representation``, attachment into ``label``,
its readers into ``io``. See that package's docstring for the table.
"""

from __future__ import annotations

import importlib

#: domain -> {public name -> the module that defines it}. **The authored map.**
BY_DOMAIN = {
    # -- io -- the formats: fps.json, its legacy ancestors, structures, templates
    "io": {
        # IMP.bff.io.structure
        "read_dcd": "IMP.bff.io.structure",
        "read_trajectory": "IMP.bff.io.structure",
        # IMP.bff.io.cif
        # IMP.bff.io.structure
        "read_rotamer_library_rmf": "IMP.bff.io.structure",
        "write_rotamer_library_rmf": "IMP.bff.io.structure",
    },
    # -- label -- a dye attached at a site, its linker and frame, and the quenchers it sees
    "label": {
        # IMP.bff.label
        "SITE_KEEP_ATOM_NAMES": "IMP.bff.label",
        "align_hierarchies": "IMP.bff.label",
        "attach_dyes": "IMP.bff.label",
        "place_dye_from_coords": "IMP.bff.label",
        "place_dye_from_rotamer_cb": "IMP.bff.label",
        "resolve_dye_site": "IMP.bff.label",
        "strip_sidechain_at_site": "IMP.bff.label",
        # IMP.bff.label
        "backbone_frame": "IMP.bff.label",
        "backbone_frame_from_coords": "IMP.bff.label",
        "backbone_transformation": "IMP.bff.label",
        "backbone_transformation_from_coords": "IMP.bff.label",
        # IMP.bff.label
        # IMP.bff.label
        "FLUOROPHORE_TYPES": "IMP.bff.label",
        "Label": "IMP.bff.label",
        # IMP.bff.label
        "select_atoms": "IMP.bff.label",
        "strip_hierarchy": "IMP.bff.label",
        "strip_obstacles": "IMP.bff.label",
    },
    # -- representation -- stage 1 -- where the dye can be: accessible volume, rotamer library, distributions
    "representation": {
        # IMP.bff.representation.rotamer
        "RotamerEnsemble": "IMP.bff.representation.rotamer",
        "SIMULATION_TYPE_R1": "IMP.bff.representation.rotamer",
        "rotamer_ensembles_from_fps": "IMP.bff.representation.rotamer",
        # IMP.bff.representation.rotamer
        "RotamerDistance": "IMP.bff.representation.rotamer",
        "RotamerPosition": "IMP.bff.representation.rotamer",
        "distances_from_ensembles": "IMP.bff.representation.rotamer",
        "read_rotamer_fps": "IMP.bff.representation.rotamer",
        "rotamer_ensemble_payload": "IMP.bff.representation.rotamer",
        "rotamer_fret_from_fps": "IMP.bff.representation.rotamer",
        "rotamer_position_payload": "IMP.bff.representation.rotamer",
        "write_rotamer_fps": "IMP.bff.representation.rotamer",
        # IMP.bff.representation.rotamer
        "RotamerFRET": "IMP.bff.representation.rotamer",
        # IMP.bff.representation.rotamer
        "load_protein_frames": "IMP.bff.representation.rotamer",
        "load_rotamer_library": "IMP.bff.representation.rotamer",
        "resolve_rotamer_library_path": "IMP.bff.representation.rotamer",
        "rotamer_library_metadata": "IMP.bff.representation.rotamer",
        "rotamer_library_registry": "IMP.bff.representation.rotamer",
    },
    # -- scoring -- stage 2 -- is this configuration allowed, and how heavily weighted
    "scoring": {
        # IMP.bff.scoring
        "compute_rotamer_score": "IMP.bff.scoring",
    },
}

#: The flat view: ``IMP.bff.<Name>``. Derived from :data:`BY_DOMAIN`, never
#: edited directly -- a name added to two domains must be an error, not a
#: silently surviving duplicate.
EXPORTS = {}
for _domain, _members in BY_DOMAIN.items():
    for _name, _module in _members.items():
        if _name in EXPORTS:
            raise RuntimeError(
                f"{_name!r} is exported by two domains: "
                f"{EXPORTS[_name]} and {_module}")
        EXPORTS[_name] = _module
del _domain, _members, _name, _module

#: The two exports whose module name differs from their public one.
#:
#: There were ten. Eight were legacy spellings -- ``kappasq_dwt`` for
#: ``kappa2_distribution_diffusion_with_traps``, ``s2delta`` for
#: ``s2_delta_from_anisotropy`` -- which meant the package's own code never read
#: the way its API did, and a reader following a public name landed on a
#: different word. Those implementations were renamed to match, and the aliases
#: deleted.
#:
#: These two remain because they are not a spelling difference: they are two
#: genuinely different functions that a name collision forced apart, and the
#: flat surface keeps the meaning it always had.
#: Names whose module spelling differs from their public one. Empty: the two
#: that were here were the two `compute_av`s, an *array* door and a *structure*
#: door that a name collision had forced apart. Both are C++ now, under the two
#: names they always meant.
_SOURCE_NAME = {}


def domain_of(name: str) -> str:
    """Which domain a public name belongs to (raises ``KeyError`` if none)."""
    for domain, members in BY_DOMAIN.items():
        if name in members:
            return domain
    raise KeyError(name)


def resolve(name: str):
    """Import and return the object behind a public name (raises AttributeError)."""
    target = EXPORTS.get(name)
    if target is None:
        raise AttributeError(f"module 'IMP.bff' has no attribute {name!r}")
    module = importlib.import_module(target)
    return getattr(module, _SOURCE_NAME.get(name, name))


def module_getattr(name: str):
    """PEP 562 hook installed on ``IMP.bff``: resolve lazily, then cache."""
    import IMP.bff
    value = resolve(name)
    setattr(IMP.bff, name, value)
    return value


def module_dir():
    """``dir(IMP.bff)`` including the lazy exports."""
    import IMP.bff
    return sorted(set(vars(IMP.bff)) | set(EXPORTS))


def public_names():
    """The sorted list of flat public names."""
    return sorted(EXPORTS)
