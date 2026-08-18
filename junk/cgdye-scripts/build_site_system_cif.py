#!/usr/bin/env python
"""Backwards-compatible CLI wrapper for cgdye topology builder."""

import sys
from pathlib import Path



from IMP.bff.cgdye.topology.builder import main


if __name__ == "__main__":
    main()
