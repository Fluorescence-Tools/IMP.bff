/*
 * The connection layer's Python surface: the IMP concepts -- the AV decorator
 * and its restraints, the probe attachment and dynamics over IMP particles,
 * docking, the fps.json project and exports that drive it, the coarse-grained
 * potentials as IMP restraints, and the Hierarchy/em/algebra bridges. Wrapped
 * after core.i, in the order these used to sit among the core's files; every
 * type they take from the core is known by then, and nothing in the core
 * takes a type from here. The standalone build does not include this file.
 */

%include "IMP_bff.av.i"

/* Where a probe sits on a structure and how it gets there, over IMP::atom
   hierarchies (formerly part of label.i). */
%include "IMP_bff.probeattachment.i"

/* The dynamics of an explicit probe: IMP::atom::Simulator over IMP particles
   (formerly part of sampling.i). */
%include "IMP_bff.dyedynamics.i"
%include "IMP_bff.probedynamics.i"

/* The FRET-restrained docking engine. */
%include "IMP_bff.docking.i"

/*
 * FPS's exports: the results table and the five files `SaveForm` writes.
 * `FPSExport.h` owns the byte layout, the number formatting and the three
 * defects of FPS's own writers that are deliberately not reproduced (the stale
 * `BestFitRotation`, the relative `_tmp.pdb`, the camera-space transforms).
 */
%include "IMP_bff.fpsexport.i"

/*
 * The project document (PRD-121 G1): structures in body order, the labelling
 * source, the selected distances, the five per-mode parameter blocks, the AV
 * globals and the current poses. After `docking.i`, whose `DockingParameters`
 * a project hands back.
 */
%include "IMP_bff.fpsproject.i"

/* The mean-position FRET restraint, and what repeated docking says about a
   model's precision. */
%include "IMP_bff.avmeandistance.i"

/* The coarse-grained protein potentials: contacts, sterics, solvation. After
   the accessible-surface kernels it composes and the residue types it reads. */
%include "IMP_bff.potentials.i"

/* The connection layer's Hierarchy and Particle overloads of core functions,
   last: they return the core's values and take the core's lattice. */
%include "IMP_bff.hierarchybridge.i"
