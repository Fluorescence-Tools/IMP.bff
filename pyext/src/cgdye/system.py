"""DyeForceFieldSystem-level helpers for cgdye N-component pipelines."""

import os
from pathlib import Path
from IMP.bff.io.cif import read_dye_forcefield_cif, write_dye_forcefield_cif


class DyeForceFieldSystem:
    """Wrapper for the cgdye system dictionary (mmCIF backed)."""

    def __init__(self, data=None):
        self.data = data or {}

    @classmethod
    def from_cif(cls, path):
        """Load a system from an mmCIF file."""
        data = read_dye_forcefield_cif(path)
        return cls(data)

    def write_cif(self, path):
        """Write the system to an mmCIF file."""
        write_dye_forcefield_cif(path, self.data)

    def get_fixed_components(self):
        """Return list of component names with role 'fixed'."""
        comps = self.data.get("components", {})
        return [name for name, spec in comps.items() if spec.role == "fixed"]

    def get_mobile_components(self):
        """Return list of component names with role 'mobile'."""
        comps = self.data.get("components", {})
        return [name for name, spec in comps.items() if spec.role == "mobile"]

    def get_component_role(self, name):
        """Return role ('fixed' or 'mobile') for a component, or None if not found."""
        spec = self.data.get("components", {}).get(name, {})
        return spec.role

    def get_fixed_component(self):
        """Return the single fixed component name, or None if not uniquely defined."""
        fixed = self.get_fixed_components()
        if len(fixed) == 1:
            return fixed[0]
        return None

    def find_group_for_component(self, component_name, group_suffix):
        """Find a group name like '{component_name}_{group_suffix}' in the system.

        Returns the group name if found, None otherwise.
        """
        groups = self.data.get("groups", {})
        target = f"{component_name}_{group_suffix}"
        return target if target in groups else None

    @property
    def name(self):
        return self.data.get("name", "unknown")

    @property
    def sites(self):
        return self.data.get("sites", [])

    @property
    def bonds(self):
        return self.data.get("bonds", [])


def fixed_components(system):
    """Return list of component names with role 'fixed'."""
    from IMP.bff.io.cif import as_forcefield_system
    system = as_forcefield_system(system)
    if isinstance(system, DyeForceFieldSystem):
        return system.get_fixed_components()
    comps = system.components
    return [name for name, spec in comps.items() if spec.role == "fixed"]


def mobile_components(system):
    """Return list of component names with role 'mobile'."""
    from IMP.bff.io.cif import as_forcefield_system
    system = as_forcefield_system(system)
    if isinstance(system, DyeForceFieldSystem):
        return system.get_mobile_components()
    comps = system.components
    return [name for name, spec in comps.items() if spec.role == "mobile"]


def component_role(system, name):
    """Return role ('fixed' or 'mobile') for a component, or None if not found."""
    from IMP.bff.io.cif import as_forcefield_system
    system = as_forcefield_system(system)
    if isinstance(system, DyeForceFieldSystem):
        return system.get_component_role(name)
    spec = system.components.get(name, {})
    return spec.role


def derive_system_name(*mol2_paths):
    """Derive system name from MOL2 file stems.

    Example: ('cx4.mol2', 'atto655.mol2') -> 'cx4_atto655'
    """
    stems = [Path(p).stem for p in mol2_paths]
    return "_".join(stems)


def fixed_component(system):
    """Return the single fixed component name, or None if not uniquely defined."""
    from IMP.bff.io.cif import as_forcefield_system
    system = as_forcefield_system(system)
    if isinstance(system, DyeForceFieldSystem):
        return system.get_fixed_component()
    fixed = fixed_components(system)
    if len(fixed) == 1:
        return fixed[0]
    return None


def find_group_for_component(system, component_name, group_suffix):
    """Find a group name like '{component_name}_{group_suffix}' in the system.

    Returns the group name if found, None otherwise.
    """
    from IMP.bff.io.cif import as_forcefield_system
    system = as_forcefield_system(system)
    if isinstance(system, DyeForceFieldSystem):
        return system.find_group_for_component(component_name, group_suffix)
    groups = system.groups
    target = f"{component_name}_{group_suffix}"
    return target if target in groups else None
