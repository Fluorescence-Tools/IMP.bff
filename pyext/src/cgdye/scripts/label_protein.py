#!/usr/bin/env python
"""Label a protein (optionally from PDB) with a dye and linker."""

import os
import sys
import urllib.request
from pathlib import Path

import click
import IMP
import IMP.atom
import numpy as np

# Ensure local package is importable


from IMP.bff.cgdye.labeling.attachment import attach_dyes, place_dye_from_coords
from IMP.bff.cgdye.io.rotamer_rmf import read_rotamer_library_rmf
from IMP.bff.cgdye.sampling.rotamer import apply_rotamer_coords
from IMP.bff.cgdye.utils import get_structure_dir, get_template_dir, _data_root

def _rotamer_library_dir():
    """Directory of the bundled rotamer library.

    It ships as IMP.bff module data rather than as a vendored package,
    so it is reached through get_data_path and not by walking up from
    __file__ -- which broke the moment the module moved.
    """
    from pathlib import Path as _Path
    import IMP.bff
    return _Path(IMP.bff.get_data_path("rotamer_library"))



def download_pdb(pdb_id, output_path):
    """Download a PDB file from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    click.echo(f"Downloading {pdb_id} from {url}...")
    try:
        urllib.request.urlretrieve(url, output_path)
        return True
    except Exception as e:
        click.echo(f"Error downloading PDB: {e}")
        return False


def find_dye_structure(dye_name, linker_type=None):
    """Search for a dye structure (MOL2, PDB or RMF) in common locations."""
    # 1. Search in cgdye/templates/rotamer/ for RMF libraries
    templates_dir = get_template_dir("rotamer")
    
    # Map common names to identifiers
    mapping = {
        "alexa488": "A48",
        "alexa350": "A35",
        "alexa532": "A53",
        "alexa568": "A56",
        "alexa594": "A59",
        "alexa647": "A64",
        "atto390": "T39",
        "atto425": "T42",
        "atto465": "T46",
        "atto488": "T48",
        "atto495": "T49",
        "atto520": "T52",
        "atto610": "T61",
        "atto655": "T65",
    }
    
    clean_dye = dye_name.lower().replace(" ", "").replace("-", "")
    dye_id = mapping.get(clean_dye, dye_name)
    
    patterns = []
    if linker_type:
        patterns.append(f"{dye_id}_{linker_type}.rmf3")
    patterns.append(f"{dye_id}_C1R.rmf3")
    patterns.append(f"{dye_id}_C2R.rmf3")
    patterns.append(f"{dye_id}_L1R.rmf3")
    
    for p in patterns:
        f = templates_dir / p
        if f.exists():
            return f, "rmf"

    # 2. Search the bundled rotamer library (IMP.bff module data) for PDB/DCD bases
    lib_dir = _rotamer_library_dir()
    if lib_dir.exists():
        pdb_patterns = []
        if linker_type:
            pdb_patterns.append(f"{dye_id}_{linker_type}.pdb")
        pdb_patterns.append(f"{dye_id}_C1R.pdb")
        pdb_patterns.append(f"{dye_id}_C2R.pdb")
        pdb_patterns.append(f"{dye_id}_L1R.pdb")
        for p in pdb_patterns:
            f = lib_dir / p
            if f.exists():
                return f, "pdb"

    # 3. Search in inputs/structures/ for MOL2 files
    structures_dir = get_structure_dir()
    mol2_patterns = []
    if linker_type:
        mol2_patterns.append(f"{dye_name}_{linker_type}.mol2")
    mol2_patterns.append(f"{dye_name}.mol2")
    
    # Common fallbacks
    if "488" in dye_name.lower():
        mol2_patterns.append("alexa488_r48.mol2")
    if "655" in dye_name.lower() or "atto655" in dye_name.lower():
        mol2_patterns.append("atto655.mol2")
        
    for p in mol2_patterns:
        f = structures_dir / p
        if f.exists():
            return f, "mol2"
            
    return None, None


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("pdb_id_or_path")
@click.option("--chain", default="E", help="Attachment chain ID")
@click.option("--residue", required=True, type=int, help="Attachment residue number")
@click.option("--dye", required=True, help="Dye name (e.g., Alexa488, Alexa647, A64)")
@click.option("--linker", help="Linker type (e.g., C1R, C2R, L1R)")
@click.option("--output", help="Output PDB path (default: {pdb}_{dye}_{residue}.pdb)")
@click.option("--strip-sidechain", is_flag=True, default=True, help="Remove site side-chain")
def main(pdb_id_or_path, chain, residue, dye, linker, output, strip_sidechain):
    """Label a protein with a dye.
    
    Example: python scripts/label_protein.py 148L --residue 132 --dye Alexa647 --linker C2R
    """
    model = IMP.Model()
    
    # Handle PDB input
    if os.path.exists(pdb_id_or_path):
        pdb_path = pdb_id_or_path
        pdb_name = Path(pdb_path).stem
    else:
        pdb_name = pdb_id_or_path.upper()
        pdb_path = f"inputs/structures/{pdb_name}.pdb"
        if not os.path.exists(pdb_path):
            os.makedirs("inputs/structures", exist_ok=True)
            if not download_pdb(pdb_id_or_path, pdb_path):
                sys.exit(1)
                
    protein = IMP.atom.read_pdb(pdb_path, model, IMP.atom.NonWaterPDBSelector())
    
    # Handle Dye input
    struct_path, struct_type = find_dye_structure(dye, linker)
    if not struct_path:
        click.echo(f"Error: Could not find structure for dye {dye} and linker {linker}")
        sys.exit(1)
        
    click.echo(f"Using dye structure: {struct_path} ({struct_type})")
    
    if struct_type == "rmf":
        lib = read_rotamer_library_rmf(str(struct_path))
        # Find the matching PDB for the RMF to use as hierarchy template
        dye_id = struct_path.stem.split("_")[0]
        base_pdb, _ = find_dye_structure(dye_id, linker)
        if not base_pdb or not str(base_pdb).endswith(".pdb"):
             # Fallback
             base_pdb = _rotamer_library_dir() / f"{struct_path.stem}.pdb"
             
        if base_pdb.exists():
            dye_hier = IMP.atom.read_pdb(str(base_pdb), model, IMP.atom.AllPDBSelector())
        else:
            dye_hier = IMP.atom.read_mol2("inputs/structures/alexa488_r48.mol2", model)
            
        try:
            apply_rotamer_coords(dye_hier, lib["coords"][1])
        except Exception as e:
            click.echo(f"Warning: Could not apply rotamer coords ({e}). Using base structure.")
    elif struct_type == "pdb":
        dye_hier = IMP.atom.read_pdb(str(struct_path), model, IMP.atom.AllPDBSelector())
    else:
        dye_hier = IMP.atom.read_mol2(str(struct_path), model)
    
    # Attachment
    try:
        attach_dyes(protein, [(dye_hier, chain, residue)], strip_sidechain)
    except Exception as e:
        click.echo(f"Error during attachment: {e}")
        sys.exit(1)
        
    # Build hierarchy for output
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "root"))
    root.add_child(protein)
    root.add_child(dye_hier)
    
    if not output:
        output = f"output/test_systems/{pdb_name}_{dye}_{residue}.pdb"
    
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    IMP.atom.write_pdb(root, output)
    click.echo(f"Wrote labeled system to {output}")


if __name__ == "__main__":
    main()
