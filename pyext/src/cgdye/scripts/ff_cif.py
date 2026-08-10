#!/usr/bin/env python
"""Backwards-compatible shim for cgdye mmCIF IO."""

import sys
from pathlib import Path



from IMP.bff.cgdye.io.cif import read_ff_system, write_ff_system

__all__ = ["read_ff_system", "write_ff_system"]
