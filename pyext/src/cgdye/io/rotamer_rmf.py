"""Read/write rotamer libraries in RMF format.

Each rotamer is stored as a frame in the RMF file.
Weights are stored as frame-level attributes (Score category).
"""

import os
import numpy as np
import IMP
import IMP.atom
import IMP.rmf
import RMF


def write_rotamer_library_rmf(path, library):
    """Write a rotamer library to an RMF file.
    
    Args:
        path: Path to the .rmf3 file
        library: Dict with 'weight' (list), 'atom_names' (list), 
                 'coords' (dict {id: array})
    """
    if not path.endswith(".rmf3"):
        path += ".rmf3"
        
    model = IMP.Model()
    root = IMP.atom.Hierarchy.setup_particle(IMP.Particle(model, "rotamers"))
    
    # Create a dummy chain and residue to hold the atoms
    chain = IMP.atom.Chain.setup_particle(IMP.Particle(model, "A"), "A")
    root.add_child(chain)
    res = IMP.atom.Residue.setup_particle(IMP.Particle(model, "DYE"), IMP.atom.ResidueType("DYE"), 1)
    chain.add_child(res)

    # Create atoms once
    atom_names = library["atom_names"]
    particles = []
    for name in atom_names:
        p = IMP.Particle(model, name)
        # Use Mass instead of Atom to preserve name in RMF
        IMP.atom.Mass.setup_particle(p, 1.0)
        IMP.core.XYZR.setup_particle(p)
        IMP.core.XYZR(p).set_radius(1.0)
        root_hier = IMP.atom.Hierarchy.setup_particle(p)
        res.add_child(root_hier)
        particles.append(p)
        
    fh = RMF.create_rmf_file(path)
    IMP.rmf.add_hierarchies(fh, [root])
    
    # Setup score category for weights
    score_cat = fh.get_category("score")
    weight_key = fh.get_key(score_cat, "weight", RMF.FloatTag())
    
    # Setup kinetic category for transitions
    kinetic_cat = fh.get_category("kinetic")
    trans_key = fh.get_key(kinetic_cat, "transition_matrix", RMF.IntsTag())
    
    root_node = fh.get_root_node()
    
    # Save transition matrix as a static value on the root node
    if "transitions" in library:
        # Flatten 2D matrix to 1D for RMF IntsTag
        flat_trans = np.array(library["transitions"]).flatten().tolist()
        root_node.set_static_value(trans_key, flat_trans)
    
    ids = sorted(library["coords"].keys())
    for i, rid in enumerate(ids):
        coords = library["coords"][rid]
        for p, c in zip(particles, coords):
            IMP.core.XYZ(p).set_coordinates(IMP.algebra.Vector3D(c[0], c[1], c[2]))
        
        # Save weight to frame
        IMP.rmf.save_frame(fh, str(rid))
        root_node.set_frame_value(weight_key, float(library["weight"][i]))
        
    del fh # Close file


def read_rotamer_library_rmf(path):
    """Read a rotamer library from an RMF file.
    
    Returns dict with 'weight', 'atom_names', 'coords', 'transitions'.
    """
    if not os.path.exists(path):
        if os.path.exists(path + ".rmf3"):
            path += ".rmf3"
        else:
            raise FileNotFoundError(f"RMF library not found: {path}")
            
    model = IMP.Model()
    fh = RMF.open_rmf_file_read_only(path)
    hs = IMP.rmf.create_hierarchies(fh, model)
    root = hs[0]
    
    # Use get_leaves to get all atoms in the hierarchy
    atoms = IMP.atom.get_leaves(root)
    # Filter to only keep particles that are likely atoms (have a name and XYZ)
    atom_names = [a.get_name() for a in atoms]
    
    score_cat = fh.get_category("score")
    weight_key = fh.get_key(score_cat, "weight", RMF.FloatTag())
    
    kinetic_cat = fh.get_category("kinetic")
    trans_key = fh.get_key(kinetic_cat, "transition_matrix", RMF.IntsTag())
    
    root_node = fh.get_root_node()
    
    n_frames = fh.get_number_of_frames()
    weights = []
    coords_dict = {}
    
    for i in range(n_frames):
        IMP.rmf.load_frame(fh, RMF.FrameID(i))
        weights.append(root_node.get_frame_value(weight_key))
        
        frame_coords = np.zeros((len(atoms), 3))
        for j, a in enumerate(atoms):
            c = IMP.core.XYZ(a).get_coordinates()
            frame_coords[j] = [c[0], c[1], c[2]]
        coords_dict[i + 1] = frame_coords
        
    # Read transition matrix if present
    transitions = None
    if root_node.get_has_value(trans_key):
        flat_trans = root_node.get_static_value(trans_key)
        n_clusters = int(np.sqrt(len(flat_trans)))
        transitions = np.array(flat_trans).reshape((n_clusters, n_clusters)).tolist()
        
    return {
        "id": list(range(1, n_frames + 1)),
        "weight": weights,
        "atom_names": atom_names,
        "coords": coords_dict,
        "transitions": transitions
    }
