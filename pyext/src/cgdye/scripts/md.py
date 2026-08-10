#!/usr/bin/env python
"""Entry-point for cgdye-backed IMP MD runner."""

import sys
from pathlib import Path



from IMP.bff.cgdye.sim.runner import main


if __name__ == "__main__":
    main()
