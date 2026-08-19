\brief Bayesian Fluorescence Framework

# imp_bff {#imp_bff}

Bayesian Fluorescence Framework (BFF) processes and analyzes 
fluorescence data. BFF contains functions to computed label 
distributions and forward models. Among other data the program 
imp_bff samples can sample biomolecular conformations restrained 
by experimental data by  simulating fluorophore distributions 
around attachment sites and  by comparing simulated observables 
to experimental data.

Labels are modelled two ways: as **accessible volumes** (`IMP.bff.AV`,
`AVNetworkRestraint`, the `fret` docking layer) or as **explicit dyes** — an
atomistic dye + linker placed on a residue from a rotamer library (the
FRETpredict libraries ship as module data) or sampled over its linker degrees
of freedom, screened against the protein and turned into R0/κ²/FRET
efficiencies (`IMP.bff.RotamerFRET`, `IMP.bff.attach_dyes`, ...). Every public
name is reachable flat as `IMP.bff.<Name>`; the manual page
`doc/manual/structure/structure_cgdye.ipynb` walks the explicit route, and
`okf/cgdye.md` records how the code is organised.


## Inter-label distance score usage:

First, import the module:

```python
import IMP.bff
import IMP.bff.restraints
```

Then, select the "score set", i.e., a set of distances that are used for 
score calculation from a FPS.JSON file:

```python
fps_json_fn = str(root_dir / "screening.fps.json")
score_set = "inter"
```

Finally, create the restraint and add it to the model.

```python
fret_restraint = IMP.bff.restraints.AVNetworkRestraintWrapper(
    hier, fps_json_fn,
    mean_position_restraint=True,
    score_set=score_set
)
fret_restraint.add_to_model()
output_objects.append(fret_restraint)
```


The command tree is one tree. `imp_bff --help` lists every command:
`flexfit` and `rmsd` fit against distance restraints, `decays` runs the
automated decay analysis, `dye` is explicit-dye labelling and sampling,
`rotamer` is rotamer-library FRET, and `av-vs-rotamer` regenerates the
comparison note. Two of those groups used to live *inside* the package and
were reachable only as `python -m IMP.bff.cgdye.cli` and
`python -m IMP.bff.cli`; a click command is a decorated function, so a library
module carrying one cannot be imported without click, which is why command
trees belong in `bin/`.

# imp_bff_traj2bcif: convert a trajectory to BinaryCIF {#imp_bff_traj2bcif}

Converts a DCD or XTC trajectory to a BinaryCIF `_atom_site` coordinate
category, which is the format the shipped rotamer libraries use. Lossless
float32 by default: the FRETpredict pins are sensitive to dipole *directions*
between atoms about 1.7 A apart, so a quantisation grid that looks harmless as
a displacement is not one as an angle. `--grid` opts into quantisation for
corpora where that does not hold. Reading a DCD needs only `IMP.bff`; reading
an XTC needs `mdtraj`, which `imp_bff` already imports.

# imp_bff_dye_pdb2cif: convert a dye PDB to mmCIF {#imp_bff_dye_pdb2cif}

Writes the `_atom_site` records for a dye structure, deriving the element from
the atom name where the PDB does not carry one. The conversion itself is
`IMP.bff.io.structure.convert_pdb_to_cif`; this is its command-line driver.

# Info

_Author(s)_: Thomas-Otavio Peulen

_Maintainer_: `tpeulen`

_License_: [LGPL](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
- None
