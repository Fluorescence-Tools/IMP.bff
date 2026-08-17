#!/usr/bin/env python
"""Backwards-compatible shim for cgdye mmCIF IO."""

import sys
from pathlib import Path



from IMP.bff.cgdye.io.cif import read_dye_forcefield_cif, write_dye_forcefield_cif

__all__ = ["read_dye_forcefield_cif", "write_dye_forcefield_cif"]
