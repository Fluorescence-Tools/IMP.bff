"""Unified CLI for cgdye."""

import os
import sys
import random
import urllib.request
from pathlib import Path

import click
import IMP
import IMP.atom
import IMP.core
import IMP.rmf
import RMF
import numpy as np

from .labeling.attachment import attach_dyes, place_dye_from_coords, resolve_dye_site
from .io.rotamer_rmf import read_rotamer_library_rmf, write_rotamer_library_rmf
from .sampling.rotamer import apply_rotamer_coordinates, sample_rotamer_index
from .sampling.kinetic import reconstruct_rotamer_trajectory, rotamer_rotational_correlation_time, rotamer_correlation_times
from .sampling.library_gen import generate_linker_rotamers
from .utils import get_template_dir, get_structure_dir, get_output_dir, ensure_dir

def _chain_sequence(hierarchy, chain_id):
    """One-letter sequence of a chain in an IMP hierarchy.

    Returns None when the chain is not present, matching the caller's
    "sequence not found" branch.
    """
    for ch in IMP.atom.get_by_type(hierarchy, IMP.atom.CHAIN_TYPE):
        if IMP.atom.Chain(ch).get_id() != chain_id:
            continue
        letters = []
        for res in IMP.atom.get_by_type(ch, IMP.atom.RESIDUE_TYPE):
            letters.append(IMP.atom.get_one_letter_code(
                IMP.atom.Residue(res).get_residue_type()))
        return "".join(letters) or None
    return None


def _rotamer_library_dir():
    """Directory of the bundled rotamer library.

    It ships as IMP.bff module data rather than as a vendored package,
    so it is reached through get_data_path and not by walking up from
    __file__ -- which broke the moment the module moved.
    """
    from pathlib import Path as _Path
    import IMP.bff
    return _Path(IMP.bff.get_data_path("rotamer_library"))



@click.group()
def dye():
    """cgdye: Coarse-Grained Dye Simulations in IMP."""
    pass


# --- label ---

def download_pdb(pdb_id, output_path):
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        urllib.request.urlretrieve(url, output_path)
        return True
    except Exception as e:
        click.echo(f"Error downloading PDB: {e}")
        return False


def find_dye_structure(dye_name, linker_type=None):
    templates_dir = get_template_dir("rotamer")
    mapping = {
        "alexa488": "A48", "alexa350": "A35", "alexa532": "A53",
        "alexa568": "A56", "alexa594": "A59", "alexa647": "A64",
        "atto390": "T39", "atto425": "T42", "atto465": "T46",
        "atto488": "T48", "atto495": "T49", "atto520": "T52",
        "atto610": "T61", "atto655": "T65",
    }
    clean_dye = dye_name.lower().replace(" ", "").replace("-", "")
    dye_id = mapping.get(clean_dye, dye_name)
    
    patterns = []
    if linker_type: patterns.append(f"{dye_id}_{linker_type}.rmf3")
    patterns.extend([f"{dye_id}_C1R.rmf3", f"{dye_id}_C2R.rmf3", f"{dye_id}_L1R.rmf3"])
    
    for p in patterns:
        f = templates_dir / p
        if f.exists(): return f, "rmf"

    lib_dir = _rotamer_library_dir()
    if lib_dir.exists():
        pdb_patterns = []
        if linker_type: pdb_patterns.append(f"{dye_id}_{linker_type}.pdb")
        pdb_patterns.extend([f"{dye_id}_C1R.pdb", f"{dye_id}_C2R.pdb", f"{dye_id}_L1R.pdb"])
        for p in pdb_patterns:
            f = lib_dir / p
            if f.exists(): return f, "pdb"

    structures_dir = get_structure_dir()
    mol2_patterns = []
    if linker_type: mol2_patterns.append(f"{dye_name}_{linker_type}.mol2")
    mol2_patterns.append(f"{dye_name}.mol2")
    if "488" in dye_name.lower(): mol2_patterns.append("alexa488_r48.mol2")
    if "655" in dye_name.lower() or "atto655" in dye_name.lower(): mol2_patterns.append("atto655.mol2")
    for p in mol2_patterns:
        f = structures_dir / Path(p)
        if f.exists(): return f, "mol2"
            
    return None, None


def resolve_protein_pdb(pdb_id_or_path):
    """Resolve a protein PDB ID or path to an existing PDB file."""
    if os.path.exists(pdb_id_or_path):
        return pdb_id_or_path
    pdb_name = pdb_id_or_path.upper()
    pdb_path = str(get_structure_dir(f"{pdb_name}.pdb"))
    if os.path.exists(pdb_path):
        return pdb_path
    ensure_dir(get_structure_dir())
    if download_pdb(pdb_id_or_path, pdb_path):
        return pdb_path
    raise FileNotFoundError(f"Unable to find or download PDB file for {pdb_id_or_path}")


@dye.command()
@click.argument("pdb_id_or_path")
@click.option("--chain", default="A", help="Chain ID")
@click.option("--residue", required=True, type=int, help="Residue number")
@click.option("--dye", required=True, help="Dye name")
@click.option("--linker", help="Linker type")
@click.option("--output", help="Output PDB path")
def label(pdb_id_or_path, chain, residue, dye, linker, output):
    """Label a protein with a dye."""
    model = IMP.Model()
    if os.path.exists(pdb_id_or_path):
        pdb_path = pdb_id_or_path
        pdb_name = Path(pdb_path).stem
    else:
        pdb_name = pdb_id_or_path.upper()
        pdb_path = str(get_structure_dir(f"{pdb_name}.pdb"))
        if not os.path.exists(pdb_path):
            ensure_dir(get_structure_dir())
            if not download_pdb(pdb_id_or_path, pdb_path): sys.exit(1)
                
    protein = IMP.atom.read_pdb(pdb_path, model, IMP.atom.NonWaterPDBSelector())
    struct_path, struct_type = find_dye_structure(dye, linker)
    if not struct_path:
        click.echo(f"Error: No structure found for {dye} {linker}")
        sys.exit(1)
        
    click.echo(f"Using: {struct_path}")
    if struct_type == "rmf":
        lib = read_rotamer_library_rmf(str(struct_path))
        dye_id = struct_path.stem.split("_")[0]
        base_pdb, _ = find_dye_structure(dye_id, linker)
        if not base_pdb or not str(base_pdb).endswith(".pdb"):
             base_pdb = _rotamer_library_dir() / f"{struct_path.stem}.pdb"
        if base_pdb.exists():
            dye_hier = IMP.atom.read_pdb(str(base_pdb), model, IMP.atom.AllPDBSelector())
        else:
            dye_hier = IMP.atom.read_mol2(str(get_structure_dir("alexa488_r48.mol2")), model)
        apply_rotamer_coordinates(dye_hier, lib["coords"][1])
    elif struct_type == "pdb":
        dye_hier = IMP.atom.read_pdb(str(struct_path), model, IMP.atom.AllPDBSelector())
    else:
        dye_hier = IMP.atom.read_mol2(str(struct_path), model)
    
    attach_dyes(protein, [(dye_hier, chain, residue)], strip_site_sidechain=True)
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "root"))
    root.add_child(protein); root.add_child(dye_hier)
    if not output: output = f"output/test_systems/{pdb_name}_{dye}_{residue}.pdb"
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    IMP.atom.write_pdb(root, output)
    click.echo(f"Saved to {output}")


# --- build-lib ---

@dye.command("build-lib")
@click.option("--n-steps", default=1000, type=int)
@click.option("--cluster-threshold", default=0.5, type=float)
@click.option("--output-dir", default=str(get_template_dir("rotamer")))
def build_lib(n_steps, cluster_threshold, output_dir):
    """Batch generate rotamer libraries."""
    dye_dir = get_structure_dir()
    mol2_files = [f for f in dye_dir.glob("*.mol2") if "1DG3" not in f.name]
    os.makedirs(output_dir, exist_ok=True)
    for mol2 in mol2_files:
        dye_name = mol2.stem.replace("_r48", "")
        try:
            lib = generate_linker_rotamers(str(mol2), n_steps=n_steps, cluster_threshold=cluster_threshold)
            out = os.path.join(output_dir, f"{dye_name}.rmf3")
            write_rotamer_library_rmf(out, lib)
            click.echo(f"Built {dye_name} ({len(lib['weight'])} rotamers)")
        except Exception as e:
            click.echo(f"Failed {dye_name}: {e}")


# --- analyze-tc ---

@dye.command("analyze-tc")
@click.argument("input_rmf", type=click.Path(exists=True))
@click.option("--dt", default=1.0, help="Timestep (ps)")
@click.option("--show-all", is_flag=True)
def analyze_tc(input_rmf, dt, show_all):
    """Calculate rotational correlation times."""
    lib = read_rotamer_library_rmf(input_rmf)
    if not lib.get("transitions"):
        click.echo("No transition matrix found.")
        return
    tc = rotamer_rotational_correlation_time(lib["transitions"], dt)
    click.echo(f"Slowest TC: {tc:.2f} ps")
    if show_all:
        times = rotamer_correlation_times(lib["transitions"], dt)
        for i, t in enumerate(times):
            if t > 0: click.echo(f"  tau_{i+2}: {t:.2f}")


# --- reconstruct ---

@dye.command()
@click.option("--lib-rmf", required=True)
@click.option("--n-frames", default=100, type=int)
@click.option("--output-rmf", required=True)
@click.option("--seed", default=42, type=int)
def reconstruct(lib_rmf, n_frames, output_rmf, seed):
    """Reconstruct trajectory from kinetic library."""
    lib = read_rotamer_library_rmf(lib_rmf)
    indices = reconstruct_rotamer_trajectory(lib, n_frames, seed=seed)
    model = IMP.Model()
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "reconstructed"))
    res = IMP.atom.Residue.setup_particle(IMP.Particle(model, "DYE"), IMP.atom.ResidueType("DYE"), 1)
    root.add_child(res)
    particles = []
    for name in lib["atom_names"]:
        p = IMP.Particle(model, name); IMP.atom.Atom.setup_particle(p, IMP.atom.AtomType("C"))
        IMP.core.XYZR.setup_particle(p); IMP.core.XYZR(p).set_radius(1.0)
        res.add_child(IMP.atom.Hierarchy.setup_particle(p)); particles.append(p)
    fh = RMF.create_rmf_file(output_rmf); IMP.rmf.add_hierarchies(fh, [root])
    for i, ridx in enumerate(indices):
        coords = lib["coords"][ridx + 1]
        for p, c in zip(particles, coords): IMP.core.XYZ(p).set_coordinates(IMP.algebra.Vector3D(c[0], c[1], c[2]))
        IMP.rmf.save_frame(fh, f"step_{i}")
    click.echo(f"Saved to {output_rmf}")


# --- sample-rotamer ---

def find_dye_mol2(dye_name):
    structures_dir = get_structure_dir()
    mapping = {"alexa488": "alexa488_r48", "atto655": "atto655"}
    clean = dye_name.lower().replace(" ", "").replace("-", "")
    target = mapping.get(clean, dye_name)
    for p in [f"{target}.mol2", f"{dye_name}.mol2"]:
        f = structures_dir / Path(p)
        if f.exists(): return f
    return structures_dir / "alexa488_r48.mol2" # Final fallback

@dye.command("sample-rotamer")
@click.option("--protein-pdb", required=True)
@click.option("--chain", default="A")
@click.option("--residue", required=True, type=int)
@click.option("--dye", required=True)
@click.option("--n-samples", default=100, type=int)
@click.option("--output-rmf", required=True)
def sample_rotamer(protein_pdb, chain, residue, dye, n_samples, output_rmf):
    """Library-based sampling with clash detection."""
    lib_path, _ = find_dye_structure(dye)
    lib = read_rotamer_library_rmf(str(lib_path))
    model = IMP.Model()
    protein = IMP.atom.read_pdb(protein_pdb, model, IMP.atom.NonWaterPDBSelector())
    # Template for attachment
    mol2_path = find_dye_mol2(dye)
    dye_hier = IMP.atom.read_mol2(str(mol2_path), model)
    # Ensure all atoms have XYZR for RMF
    for a in IMP.atom.get_leaves(dye_hier):
        if not IMP.core.XYZR.get_is_setup(a):
            IMP.core.XYZR.setup_particle(a, 1.0)
            
    attached = attach_dyes(protein, [(dye_hier, chain, residue)], True)
    site = attached[0]["site"]
    ca, n, c = [IMP.core.XYZ(site[k]).get_coordinates() for k in ["CA", "N", "C"]]
    
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "root"))
    root.add_child(protein); root.add_child(dye_hier)
    fh = RMF.create_rmf_file(output_rmf); IMP.rmf.add_hierarchies(fh, [root])
    
    rng = random.Random(42)
    for i in range(n_samples):
        ridx = sample_rotamer_index(lib["weight"], rng=rng)
        apply_rotamer_coordinates(dye_hier, lib["coords"][ridx + 1])
        place_dye_from_coords(dye_hier, ca, n, c)
        IMP.rmf.save_frame(fh, str(i))
    click.echo(f"Wrote {n_samples} frames to {output_rmf}")


# --- sample-dof-walk ---

def _protein_residue_groups(protein_hier):
    groups = {}
    for a in IMP.atom.get_by_type(protein_hier, IMP.atom.ATOM_TYPE):
        res_p = a.get_parent()
        if not IMP.atom.Residue.get_is_setup(res_p): continue
        res = IMP.atom.Residue(res_p)
        chain_p = res_p.get_parent()
        if not IMP.atom.Chain.get_is_setup(chain_p): continue
        chain = IMP.atom.Chain(chain_p)
        key = (chain.get_id(), int(res.get_index()))
        groups.setdefault(key, []).append(a)
    return groups

def _center_and_radius(coords):
    if not coords: return None
    c = np.mean(coords, axis=0)
    r = np.max(np.linalg.norm(coords - c, axis=1))
    return (c[0], c[1], c[2], r)

def _protein_atom_obstacles(protein_hier, site_ca, interaction_sphere, ignore_key=None):
    groups = _protein_residue_groups(protein_hier)
    obstacles = []
    for (cid, ridx), atoms in groups.items():
        if ignore_key is not None and (cid, ridx) == ignore_key:
            continue
        for atom in atoms:
            coord = IMP.core.XYZ(atom).get_coordinates()
            if np.linalg.norm(np.array(coord) - site_ca) <= interaction_sphere:
                obstacles.append((coord[0], coord[1], coord[2], 1.0))
    return obstacles

@dye.command("sample-dof-walk")
@click.option("--protein-pdb", required=True)
@click.option("--chain", default="A")
@click.option("--residue", required=True, type=int)
@click.option("--dye", required=True)
@click.option("--linker", help="Linker type")
@click.option("--n-steps", default=1000, type=int)
@click.option("--output-rmf", required=True)
@click.option("--seed", default=42, type=int)
def sample_dof_walk(protein_pdb, chain, residue, dye, linker, n_steps, output_rmf, seed):
    """Collision-gated Metropolis walk over the linker's internal DOFs.

    Not Langevin dynamics: no forces, friction or temperature -- a proposal is
    accepted whenever the dye does not clash. Real Langevin/Brownian dynamics
    is ``sample-langevin`` (IMP.bff.cgdye.sampling.langevin).
    """
    from .sampling.library_gen import LinkerSampler
    random.seed(seed); np.random.seed(seed)
    protein_path = resolve_protein_pdb(protein_pdb)
    model = IMP.Model()
    protein = IMP.atom.read_pdb(protein_path, model, IMP.atom.NonWaterPDBSelector())
    
    # 1. Try to find a PDB for the dye first (RMF needs Residue/Chain)
    struct_path, struct_type = find_dye_structure(dye, linker)
    if struct_type == "rmf":
        dye_id = struct_path.stem.split("_")[0]
        base_pdb, _ = find_dye_structure(dye_id, linker)
        if base_pdb and str(base_pdb).endswith(".pdb"):
            struct_path, struct_type = base_pdb, "pdb"
        else:
            struct_path, struct_type = get_structure_dir("alexa488_r48.mol2"), "mol2"
            
    if struct_type == "pdb":
        dye_hier = IMP.atom.read_pdb(str(struct_path), model, IMP.atom.AllPDBSelector())
        mol2_path = find_dye_mol2(dye)
    else:
        # If we only have MOL2, we MUST wrap it in a Residue for RMF
        dye_mol2 = str(struct_path) if struct_path else str(get_structure_dir("alexa488_r48.mol2"))
        mol2_path = dye_mol2
        raw_hier = IMP.atom.read_mol2(dye_mol2, model)
        dye_hier = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "DYE_ROOT"))
        d_chain = IMP.atom.Chain.setup_particle(IMP.Particle(model, "D"), "D")
        dye_hier.add_child(d_chain)
        d_res = IMP.atom.Residue.setup_particle(IMP.Particle(model, "DYE"), IMP.atom.ResidueType("DYE"), 1)
        d_chain.add_child(d_res)
        for a in IMP.atom.get_leaves(raw_hier):
            d_res.add_child(a)
            
    sampler = LinkerSampler(str(mol2_path))
    site = resolve_dye_site(protein, chain, residue)
    ca, n, c = [IMP.core.XYZ(site[k]).get_coordinates() for k in ["CA", "N", "C"]]
    site_ca = np.array(ca)
        
    obstacles = _protein_atom_obstacles(protein, site_ca, 35.0, ignore_key=(chain, residue))
    
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "root"))
    root.add_child(protein); root.add_child(dye_hier)
    fh = RMF.create_rmf_file(output_rmf); IMP.rmf.add_hierarchies(fh, [root])
    
    def place_dye():
        place_dye_from_coords(dye_hier, ca, n, c)

    def has_collision():
        for atom in IMP.atom.get_by_type(dye_hier, IMP.atom.ATOM_TYPE):
            c = np.array(IMP.core.XYZ(atom).get_coordinates())
            for ox, oy, oz, orad in obstacles:
                if np.linalg.norm(c - [ox, oy, oz]) < (1.2 * 0.8 + orad):
                    return True
        return False

    current_cfg = [0.0] * (len(sampler.rot_bonds) + len(sampler.rot_angles))
    accepted = 0
    for step in range(n_steps):
        proposal = [v + random.gauss(0, 0.1) for v in current_cfg]
        sampler.apply_config(proposal)
        place_dye()
        if not has_collision():
            current_cfg = proposal; accepted += 1
        else:
            sampler.apply_config(current_cfg)
            place_dye()
        if (step+1) % 10 == 0: IMP.rmf.save_frame(fh, str(step))
        
    click.echo(f"Finished: {accepted}/{n_steps} accepted. Saved to {output_rmf}")


# --- label-fp ---

def get_residue_range(hier):
    res_indices = []
    for a in IMP.atom.get_leaves(hier):
        res_p = a.get_parent()
        if IMP.atom.Residue.get_is_setup(res_p):
            res_indices.append(IMP.atom.Residue(res_p).get_index())
    if not res_indices: return 1, 1
    return min(res_indices), max(res_indices)

@dye.command("label-fp")
@click.argument("pdb_id_or_path")
@click.option("--site", multiple=True, required=True, help="Attachment site in format residue:fp[:chain][:anchor] (e.g., 6:eGFP:A:C)")
@click.option("--output", help="Output PDB path")
def label_fp(pdb_id_or_path, site, output):
    """Label a protein with one or more Fluorescent Proteins (FPs)."""
    from .labeling.attachment import align_hierarchies
    model = IMP.Model()
    if os.path.exists(pdb_id_or_path):
        pdb_path = pdb_id_or_path
        pdb_name = Path(pdb_path).stem
    else:
        pdb_name = pdb_id_or_path.upper()
        pdb_path = str(get_structure_dir(f"{pdb_name}.pdb"))
        if not os.path.exists(pdb_path):
            ensure_dir(get_structure_dir())
            if not download_pdb(pdb_id_or_path, pdb_path): sys.exit(1)
                
    protein = IMP.atom.read_pdb(pdb_path, model, IMP.atom.NonWaterPDBSelector())
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "root"))
    root.add_child(protein)
    
    for s in site:
        parts = s.split(":")
        res_num = int(parts[0])
        fp_type = parts[1]
        chain_id = parts[2] if len(parts) > 2 and parts[2] not in ("N", "C") else "A"
        anchor = "N"
        if len(parts) > 3: anchor = parts[3]
        elif len(parts) > 2 and parts[2] in ("N", "C"): anchor = parts[2]
        
        fp_path = get_template_dir("fp") / f"{fp_type}.pdb"
        if not fp_path.exists():
            for f in fp_path.parent.glob("*.pdb"):
                if f.stem.lower() == fp_type.lower():
                    fp_path = f; break
                    
        if not fp_path.exists():
            click.echo(f"Error: FP template {fp_type} not found.")
            continue
            
        fp_hier = IMP.atom.read_pdb(str(fp_path), model, IMP.atom.AllPDBSelector())
        r_min, r_max = get_residue_range(fp_hier)
        fp_anchor_res = r_min if anchor.upper() == "N" else r_max
        
        click.echo(f"Attaching {fp_type} (anchor {anchor}) to {chain_id}:{res_num}...")
        
        # Get target backbone coords
        try:
            from .labeling.attachment import resolve_dye_site
            t_site = resolve_dye_site(protein, chain_id, res_num)
            t_ca, t_n, t_c = [IMP.core.XYZ(t_site[k]).get_coordinates() for k in ["CA", "N", "C"]]
            
            # Align FP's anchor residue to target site
            # Templates are renumbered to 1..N and use chain 'A'
            align_hierarchies(fp_hier, "A", fp_anchor_res, t_ca, t_n, t_c)
            
            root.add_child(fp_hier)
        except Exception as e:
            click.echo(f"  Error: {e}")
            
    if not output: output = f"output/test_systems/{pdb_name}_labeled_fp.pdb"
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    IMP.atom.write_pdb(root, output)
    click.echo(f"Saved to {output}")


# --- label-fusion ---

@dye.command("label-fusion")
@click.argument("pdb_path")
@click.option("--chain", default="A", help="Chain to analyze")
@click.option("--output", help="Output labeled PDB")
def label_fusion(pdb_path, chain, output):
    """Automatically detect and label FPs in a fusion protein."""
    from .sampling.segments import find_fp_domains, parse_plddt_from_pdb, segments_from_plddt
    model = IMP.Model()
    protein = IMP.atom.read_pdb(pdb_path, model, IMP.atom.NonWaterPDBSelector())

    # Sequence comes from the hierarchy IMP just built, rather than from a
    # second parse by another library: IMP.bff carries no dependency beyond
    # what IMP itself brings.
    seq = _chain_sequence(protein, chain)
    
    if not seq:
        click.echo(f"Error: Could not find sequence for chain {chain} in {pdb_path}")
        sys.exit(1)
        
    click.echo(f"Analyzing chain {chain} (length {len(seq)})...")
    
    # 1. Detect FP domains
    fp_domains = find_fp_domains(seq)
    if not fp_domains:
        click.echo("No FP domains detected.")
        sys.exit(0)
        
    for name, s, e, idy in fp_domains:
        click.echo(f"  Found {name} at {s}-{e} (identity: {idy:.2f})")
        
    # 2. Parse pLDDT
    try:
        plddt = parse_plddt_from_pdb(pdb_path, chain)
    except Exception as e:
        click.echo(f"Warning: Could not parse pLDDT ({e}). Using dummy values.")
        plddt = {i+1: 100.0 for i in range(len(seq))}
        
    # 3. Create segments (defines rigid cores and flexible linkers)
    segs = segments_from_plddt(len(seq), plddt, fp_domains)
    
    # 4. Attach FP structures at detected domains
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "root"))
    root.add_child(protein)
    
    labeled_any = False
    for seg in segs:
        if seg["kind"] == "fp":
            fp_name = seg["name"]
            click.echo(f"  Labeling {fp_name} at residue {seg['start']}...")
            
            # Find matching template
            template_path = get_template_dir("fp") / f"{fp_name}.pdb"
            if not template_path.exists():
                # Try case-insensitive fallback
                for f in template_path.parent.glob("*.pdb"):
                    if f.stem.lower() == fp_name.lower():
                        template_path = f; break
            
            if template_path.exists():
                fp_hier = IMP.atom.read_pdb(str(template_path), model, IMP.atom.AllPDBSelector())
                try:
                    attach_dyes(protein, [(fp_hier, chain, seg["start"])], True)
                    root.add_child(fp_hier)
                    labeled_any = True
                except Exception as e:
                    click.echo(f"    Error attaching {fp_name}: {e}")
            else:
                click.echo(f"    Warning: No template found for {fp_name}")

    if not labeled_any:
        click.echo("No domains were labeled.")
        sys.exit(0)

    if not output:
        output = Path(pdb_path).stem + "_labeled.pdb"
    
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    IMP.atom.write_pdb(root, output)
    click.echo(f"Wrote labeled fusion system to {output}")


if __name__ == "__main__":
    main()
