import os
import sys

# Add local path to be able to import cgdye
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from IMP.bff.cgdye.io.nmr_cif import write_nmr_restraints

# Host sCX4 Ha atoms (NMR naming: Ha)
ha_atoms = {
    # Subunit 1
    "H61": "Ha(1)", "H57": "Ha(1)", 
    # Subunit 2
    "H58": "Ha(2)", "H62": "Ha(2)", 
    # Subunit 3
    "H60": "Ha(3)", "H64": "Ha(3)", 
    # Subunit 4
    "H63": "Ha(4)", "H59": "Ha(4)"
}

# Guest Atto655 NOESY atoms
noesy_groups = {
    12: { "atoms": ["H33"], "limit": 3.5, "error_pos": 0.5, "intensity": "Strong" },
    6:  { "atoms": ["H7"], "limit": 3.55, "error_pos": 0.5, "intensity": "Weak" },
    9:  { "atoms": ["H6"], "limit": 3.5, "error_pos": 0.5, "intensity": "Medium" },
    27: { "atoms": ["H3", "H1", "H2"], "limit": 3.5, "error_pos": 0.5, "intensity": "Strong" },
    24: { "atoms": ["H37", "H36"], "limit": 5.0, "error_pos": 0.5, "intensity": "Weak" },
    23: { "atoms": ["H34", "H35"], "limit": 5.0, "error_pos": 0.5, "intensity": "Medium" },
    25: { "atoms": ["H39", "H38"], "limit": 5.0, "error_pos": 0.5, "intensity": "Medium" },
    26: { "atoms": ["H5", "H4"], "limit": 5.0, "error_pos": 0.5, "intensity": "Medium" }
}

# Guest Atto655 CSP (repulsive) atoms
csp_groups = {
    3:  { "atoms": ["H14", "H15"], "limit": 4.5, "error_neg": 0.5 },
    4:  { "atoms": ["H12", "H13"], "limit": 4.5, "error_neg": 0.5 },
    6:  { "atoms": ["H7"], "limit": 4.5, "error_neg": 0.5 },
    18: { "atoms": ["H10", "H11"], "limit": 4.5, "error_neg": 0.5 },
    20: { "atoms": ["H25", "H27", "H26"], "limit": 4.5, "error_neg": 0.5 },
    21: { "atoms": ["H30", "H28", "H29"], "limit": 4.5, "error_neg": 0.5 }
}

data = {
    "default_restraint_set": "",
    "restraint_sets": {
        "common_prefix": "CSP_",
        "single_s1_prefix": "S1ONLY_",
        "single_s1_common_prefix": ""
    },
    "Positions": {},
    "Distances": {}
}

# Generate unique restraint set prefixes for each NOESY peak
for peak_id in noesy_groups.keys():
    data["restraint_sets"][f"subunit_peak{peak_id}_prefix"] = f"SUB_PEAK{peak_id}_"


# Add Host Atoms to Positions
for atom in ha_atoms.keys():
    pos_name = f"SCX4_{atom}"
    data["Positions"][pos_name] = {
        "simulation_type": "Atom",
        "atom_name": atom,
        "residue_name": "CX4",
        "component_id": "fixed",
        "chain_identifier": "A"
    }

added_guest_atoms = set()
def add_guest_atom(atom):
    if atom not in added_guest_atoms:
        pos_name = f"DYE_{atom}"
        data["Positions"][pos_name] = {
            "simulation_type": "Atom",
            "atom_name": atom,
            "residue_name": "DDS",
            "component_id": "mobile",
            "chain_identifier": "B"
        }
        added_guest_atoms.add(atom)

comments = {}
rid_counter = 1

# Generate NOESY distances (Upper Bound)
for nmr_guest_id, group_info in noesy_groups.items():
    for dye_atom in group_info["atoms"]:
        add_guest_atom(dye_atom)
        for ha_atom, ha_nmr_name in ha_atoms.items():
            dist_name = f"SUB_PEAK{nmr_guest_id}_{rid_counter:04d}"
            data["Distances"][dist_name] = {
                "distance": group_info["limit"],
                "error_pos": group_info["error_pos"],
                "error_neg": 1.0, 
                "position1_name": f"DYE_{dye_atom}",
                "position2_name": f"SCX4_{ha_atom}",
                "distance_type": "AtomUpperBound"
            }
            # Comment exactly as requested: NMR naming scheme
            comments[dist_name] = f"NOESY: Atto655 '{nmr_guest_id}' <-> sCx4 '{ha_nmr_name}'"
            rid_counter += 1

# Generate CSP distances (Lower Bound)
for nmr_guest_id, group_info in csp_groups.items():
    for dye_atom in group_info["atoms"]:
        add_guest_atom(dye_atom)
        for ha_atom, ha_nmr_name in ha_atoms.items():
            dist_name = f"CSP_PEAK{nmr_guest_id}_{rid_counter:04d}"
            data["Distances"][dist_name] = {
                "distance": group_info["limit"], 
                "error_pos": 1.0,
                "error_neg": group_info["error_neg"],
                "position1_name": f"DYE_{dye_atom}",
                "position2_name": f"SCX4_{ha_atom}",
                "distance_type": "AtomLowerBound"
            }
            # Comment exactly as requested: NMR naming scheme
            comments[dist_name] = f"CSP Repulsion: Atto655 '{nmr_guest_id}' <-/-> sCx4 '{ha_nmr_name}'"
            rid_counter += 1

out_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "inputs", "restraints", "nmr.restraints.cif"))
write_nmr_restraints(out_path, data)

# Map dynamic ids from cgdye's sorted outputs
sorted_dnames = sorted(data["Distances"].keys())
file_row_id_to_comment = {i + 1: comments[dname] for i, dname in enumerate(sorted_dnames)}

# Post-processing to add comments into the CIF file
with open(out_path, "r") as f:
    lines = f.readlines()

with open(out_path, "w") as f:
    columns = 0
    in_distance_loop = False
    
    for line in lines:
        if line.startswith("data_nmr_restraints"):
            f.write(line)
            f.write("#\n")
            f.write("# ------------------------------------------------------------\n")
            f.write("# NMR Restraints using original NMR naming scheme\n")
            f.write("# sCX4 atoms mapping   : Ha(1-4)\n")
            f.write("# Atto655 atoms mapping: 1-27\n")
            f.write("# ------------------------------------------------------------\n")
            continue
            
        if line.startswith("_ihm_derived_distance_restraint."):
            columns += 1
            in_distance_loop = True
            f.write(line)
            continue
            
        if line.startswith("loop_"):
            in_distance_loop = False
            f.write(line)
            continue
            
        if in_distance_loop and columns > 0 and line.strip() and not line.startswith("#"):
            # We are inside the restraint loop values
            parts = line.split()
            if len(parts) >= columns:
                try:
                    row_id = int(parts[0])
                    if row_id in file_row_id_to_comment:
                        line = line.rstrip() + f" # {file_row_id_to_comment[row_id]}\n"
                except ValueError:
                    pass
            f.write(line)
            continue

        if columns > 0 and line.strip() == "#" and in_distance_loop:
            columns = 0
            in_distance_loop = False
            
        f.write(line)

print(f"Successfully generated new nmr.restraints.cif at {out_path} with {rid_counter-1} restraints and NMR-naming comments.")
