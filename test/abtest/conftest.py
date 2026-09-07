"""Put the frozen Python chinet reference on the path for the A/B suite only."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
