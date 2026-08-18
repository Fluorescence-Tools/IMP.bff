"""The legacy C# FPS ``.txt`` formats — **read only**.

A format nobody can still read is data that has been lost, and a decade of
measurements live in these files. So they are readable, and deliberately not
writable: nothing should produce another one.

Split out of ``IMP.bff.fret.io`` by PRD-113 stage 7.
"""

from __future__ import annotations

import json
import os
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from IMP.bff.io import fps_schema

__all__ = ["read_old_lps_txt", "read_old_distances_txt"]


def read_old_lps_txt(
    path: str | os.PathLike,
    pdb_paths: Optional[List[str]] = None,
) -> Tuple[Dict, List[str]]:
    """Read a C# FPS format labeling positions (.txt) file.

    Parameters
    ----------
    path : str or PathLike
        Path to the labeling positions text file.
    pdb_paths : list of str, optional
        Paths to PDB files to resolve atom IDs to chain/residue info.

    Returns
    -------
    positions : dict
        Standard positions dict.
    molecules : list of str
        The unique list of molecule names in order of appearance.
    """
    if pdb_paths is None:
        pdb_paths = []

    if not pdb_paths:
        dir_name = os.path.dirname(str(path))
        try:
            if dir_name:
                pdb_paths = [os.path.join(dir_name, f) for f in os.listdir(dir_name) if f.endswith(".pdb")]
            else:
                pdb_paths = [f for f in os.listdir(".") if f.endswith(".pdb")]
        except Exception:
            pass

    # Map molecule name (basename without ext) to its full PDB path
    pdb_map = {}
    for p in pdb_paths:
        base = os.path.splitext(os.path.basename(p))[0].lower()
        pdb_map[base] = p

    positions = {}
    molecules = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue

            name = parts[0]
            mol = parts[1]
            dye = parts[2]
            av_type = parts[3]

            if mol not in molecules:
                molecules.append(mol)
            body_id = molecules.index(mol)

            # Resolve atom ID to chain, residue sequence number, and atom name
            chain = ""
            resseq = 1
            aname = "CA"
            atom_id = None

            # Look up PDB file if available
            pdb_file = pdb_map.get(mol.lower())

            # Parse parameters
            r2 = 0.0
            r3 = 0.0
            if av_type == "AV1" and len(parts) >= 8:
                ll = float(parts[4])
                lw = float(parts[5])
                r1 = float(parts[6])
                atom_id = int(parts[7])
            elif av_type == "AV3" and len(parts) >= 10:
                ll = float(parts[4])
                lw = float(parts[5])
                r1 = float(parts[6])
                r2 = float(parts[7])
                r3 = float(parts[8])
                atom_id = int(parts[9])
            elif av_type == "XYZ" and len(parts) >= 7:
                # Fixed coordinates
                x = float(parts[4])
                y = float(parts[5])
                z = float(parts[6])
                positions[name] = {
                    "simulation_type": "XYZ",
                    "x": x,
                    "y": y,
                    "z": z,
                    "body_id": body_id,
                }
                continue
            else:
                continue

            if atom_id is not None and pdb_file and os.path.exists(pdb_file):
                try:
                    with open(pdb_file, "r") as pf:
                        for pline in pf:
                            if pline.startswith(("ATOM  ", "HETATM")):
                                try:
                                    cur_id = int(pline[6:11].strip())
                                except ValueError:
                                    continue
                                if cur_id == atom_id:
                                    chain = pline[21].strip()
                                    try:
                                        resseq = int(pline[22:26].strip())
                                    except ValueError:
                                        resseq = 1
                                    aname = pline[12:16].strip()
                                    break
                except Exception:
                    pass
            elif atom_id is not None:
                # Fallback proxy if PDB not parsed
                resseq = atom_id

            positions[name] = {
                "chain_identifier": chain,
                "residue_seq_number": resseq,
                "atom_name": aname,
                "linker_length": ll,
                "linker_width": lw,
                "radius1": r1,
                "radius2": r2,
                "radius3": r3,
                "simulation_grid_resolution": 1.5,
                "simulation_type": av_type,
                "body_id": body_id,
            }

    return positions, molecules


def read_old_distances_txt(
    path: str | os.PathLike,
) -> Dict:
    """Read a C# FPS format experimental distances (.txt) file.

    Parameters
    ----------
    path : str or PathLike

    Returns
    -------
    distances : dict
        Standard distances dict.
    """
    distances = {}
    distance_type = "RDAMean"

    with open(path, "r") as f:
        # Check first line for distance type
        first_line = f.readline().strip()
        if first_line and not first_line.split()[0].isalnum():
            # e.g. "RDAMeanE"
            pass
        elif first_line:
            parts = first_line.split()
            if len(parts) == 1:
                distance_type = parts[0]
            else:
                # Re-wind or process as first distance
                f.seek(0)

        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue

            pos1 = parts[0]
            pos2 = parts[1]
            dist = float(parts[2])
            err_neg = float(parts[3])
            err_pos = float(parts[4])
            forster = float(parts[5]) if len(parts) >= 6 else 52.0

            dname = f"{pos1}_{pos2}"
            distances[dname] = {
                "position1_name": pos1,
                "position2_name": pos2,
                "distance": dist,
                "error_neg": err_neg,
                "error_pos": err_pos,
                "distance_type": distance_type,
                "Forster_radius": forster,
            }

    return distances
