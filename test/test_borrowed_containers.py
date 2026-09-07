"""A container read off a temporary must not be a view of freed memory.

An accessor that returns `const std::map&` or `const std::vector&` hands SWIG
a **borrowed** pointer into the owner, and the proxy Python gets does not keep
that owner alive. So this, which is how anyone would write it:

    IMP.bff.read_forcefield_cif(path).components.values()

read freed memory: the system is a temporary, it dies as the expression
unwinds, and the map proxy outlives it. It returned a path of twelve zero
bytes on a good day and a segmentation fault on a bad one -- and *silently*,
which is worse than the crash.

The accessors cannot simply return by value: `Scoring.cpp` calls
`system.get_bonds()` in the condition of a loop over the bonds, and a copy per
iteration is quadratic. The copy is made at the language boundary instead --
`%owned_container_out` in `IMP_bff.types.i` -- so it costs one copy per Python
access, and the lifetime becomes Python's.
"""

import tempfile

import pytest

import IMP.bff

#: `create_forcefield_system` is overloaded, so SWIG allows no keywords: every
#: parameter up to `relative_to` is positional, and these are its own defaults.
DEFAULTS = (2000.0, 400.0, 12.0, 1.5, 40.0, 180.0, 120.0, 220.0,
            20000, 100, 1.7, 12.0, 5.0, 6.0, 200)


@pytest.fixture(scope="module")
def system_cif():
    tmp = tempfile.mkdtemp()
    mol2 = str(IMP.bff.get_structure_dir("atto655.mol2"))
    system = IMP.bff.create_forcefield_system(
        [IMP.bff.FFComponentSpec("dye", mol2, "", "mobile")], *DEFAULTS, "")
    path = f"{tmp}/system.cif"
    IMP.bff.write_probe_forcefield_cif(path, system)
    return path, mol2


def test_a_map_of_values_survives_its_owner(system_cif):
    path, mol2 = system_cif
    # the system is a temporary and is gone by the time `values()` is read
    components = IMP.bff.read_forcefield_cif(path).components.values()
    assert [c.mol2_path for c in components] == [mol2]


def test_a_map_of_strings_survives_its_owner(system_cif):
    path, _ = system_cif
    groups = IMP.bff.read_forcefield_cif(path).groups
    assert "dye_all" in list(groups)


def test_a_vector_of_values_survives_its_owner(system_cif):
    path, _ = system_cif
    sites = IMP.bff.read_forcefield_cif(path).sites
    assert len(sites) == 70
    assert all(s.id for s in sites)


def test_a_vector_of_strings_survives_its_owner(system_cif):
    path, _ = system_cif
    assert list(IMP.bff.read_forcefield_cif(path).fixed_groups) == []


def test_a_map_of_doubles_survives_its_owner(system_cif):
    path, _ = system_cif
    bond_types = IMP.bff.read_forcefield_cif(path).bond_types
    assert "B1" in list(bond_types)
    assert bond_types["B1"] > 0.0


def test_the_copy_is_a_copy(system_cif):
    """Each access hands out its own; writing to one does not reach the C++."""
    path, _ = system_cif
    system = IMP.bff.read_forcefield_cif(path)
    first = system.groups
    first["invented"] = ["nothing"]
    assert "invented" not in list(system.groups)
