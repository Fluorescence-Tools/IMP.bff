# The connection layer to IMP

Everything built from this directory is an IMP *concept* -- a decorator, a
restraint, a scoring function, an `IMP::atom` simulator, an assembly of rigid
bodies -- and exists only when IMP is present. The rest of the module is the
core: it links against IMP's kernel and algebra today and is being made
independent of them (PRD-137), but its own API never names a `Model`, a
`Particle` or a `Hierarchy`.

Sources live here and reach the IMP build through one flat file,
`src/ImpLayer.cpp`, which includes them: IMP's tooling compiles `src/*.cpp`
and `src/internal/*.cpp` into the module and leaves other subdirectories
alone (`src/standalone/` relies on that to stay *out* of the IMP build;
this directory needs the opposite). `ImpLayer.cpp` is their only entry in
`src/Files.cmake`. Their headers stay **flat** in `include/` -- IMP's
tooling links only `include/*.h` and `include/internal/*.h` into the build
tree, so a public subdirectory is not an option -- and are named in
`Headers.cmake` here, the list the standalone build reads to leave them out.

| header | source | is |
|---|---|---|
| `AV.h` | `AV.cpp` | the accessible-volume **decorator** on a labelled particle, the sampler over `IMP::Model`, and the doors that read a PDB with `IMP::atom` (`resample_av`, `get_av_from_structure`, `get_avs_for_structure`); the array door `get_av` is core |
| `AVOccupancyMap.h` | `AVOccupancyMap.cpp` | the particle view of the core's `OccupancyGrid` raster, and the registry that shares one raster between volumes over the same particles |
| `ProbeNetworkRestraint.h` | `ProbeNetworkRestraint.cpp` | an `IMP::Restraint` over a network of AVs |
| `AVMeanDistanceRestraint.h` | `AVMeanDistanceRestraint.cpp` | an `IMP::Restraint` on the mean AV distance |
| `Potentials.h` | `Potentials.cpp` | `IMP::Restraint` factories for coarse-grained potentials, and the restraint factories over a typed dye system (formerly the IMP half of `Scoring.h`) |
| `Docking.h` | `Docking.cpp` | rigid-body docking assemblies over `IMP::atom::Hierarchy` |
| `ProbeDynamics.h` | `ProbeDynamics.cpp` | `IMP::atom::Simulator` (Langevin / Brownian) for an attached probe |
| `FPSProject.h`, `FPSExport.h` | `FPSProject.cpp`, `FPSExport.cpp` | the fps.json project and its exports, which drive `Docking.h` |
| `ProbeAttachment.h` | `ProbeAttachment.cpp` | attaching a probe hierarchy to a protein hierarchy: site resolution, backbone frame, placement, alignment (all over `IMP::atom::Hierarchy`) |
| `HierarchyBridge.h` | `HierarchyBridge.cpp` | the `IMP::atom::Hierarchy` / `IMP::Particle` overloads of core functions -- a `ProteinFrame` from a hierarchy, a PDB into a Model, a selection expression on a hierarchy, the strip mask, coordinates written back, Olga's radii per particle, a `PathMap`'s spheres from particles |

The SWIG topic files that wrap these -- `IMP_bff.av.i`, `avmeandistance.i`,
`potentials.i`, `scoring.i`, `docking.i`, `sampling.i` (its `ProbeDynamics`
part), `fpsexport.i`, `fpsproject.i` -- are the layer's Python surface and
stay in `pyext/` at the positions SWIG's ordering needs.

No core *header* names an IMP particle, hierarchy or model any more, and no
core file includes a layer header (PRD-137 step 5); a core file must not gain
one. `include/PathMap.h` names the decorator as a friend by forward
declaration only.

What is left (step 5c) is implementation, not interface: four core sources
still *read a PDB through `IMP::atom`* behind an IMP-free declaration --
`load_protein_frames` (HierarchyFrame.cpp), `load_structure` (StructureIO.cpp),
the PDB branch of the trajectory loader (ProbeSampling.cpp), and
`selection_from_expression`, which is declared here but defined in
SelectionExpression.cpp because it compiles the parser's private AST against
a hierarchy. Each is marked `PRD-137 step 5c residue` at its definition. The
replacement for the readers is the core's own `read_pdb_records`
(AVBuilder.h); the selection needs its AST in an internal header.
