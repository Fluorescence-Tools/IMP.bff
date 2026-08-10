"""Fluorescent Protein (FP) library management."""

import json
from pathlib import Path
from typing import Dict, List, Tuple

def load_fp_library() -> Dict[str, Dict]:
    """Load the FP library from the bundled JSON file."""
    path = Path(__file__).parent / "fp_library.json"
    with open(path, "r") as f:
        return json.load(f)

def get_fp_names() -> List[str]:
    """Return a list of all supported FP names."""
    return list(load_fp_library().keys())

def get_fp_data(name: str) -> Dict:
    """Return the data for a specific FP."""
    lib = load_fp_library()
    if name not in lib:
        raise ValueError(f"FP '{name}' not found in library. Available: {', '.join(lib.keys())}")
    return lib[name]
