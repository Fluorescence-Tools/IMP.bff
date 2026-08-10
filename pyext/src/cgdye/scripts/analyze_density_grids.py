#!/usr/bin/env python
"""Backwards-compatible CLI wrapper for cgdye density analysis."""

import sys
from pathlib import Path



from IMP.bff.cgdye.analysis.density import main


if __name__ == "__main__":
    main()
