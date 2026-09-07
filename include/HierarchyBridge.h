/**
 *  \file IMP/bff/HierarchyBridge.h
 *  \brief The core's values read from, and written to, an IMP::atom::Hierarchy.
 *
 * Every function here is the `IMP::atom::Hierarchy` / `IMP::Particle`
 * overload of something the core does over arrays and its own values: a
 * #IMP::bff::ProteinFrame from a hierarchy, a PDB read into an `IMP::Model`,
 * a selection expression evaluated on a hierarchy, the FPS strip mask applied
 * to one, coordinates written back into one, Olga's radii for particles, and
 * the obstacle spheres of a #IMP::bff::PathMap taken from particles. They
 * were declared in the core headers beside their array twins and moved here
 * in PRD-137 step 5, so that the core's headers name no IMP particle and the
 * connection layer is the one place that does.
 *
 * Each section is marked with the header it came from; the Python names are
 * unchanged.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_HIERARCHYBRIDGE_H
#define IMPBFF_HIERARCHYBRIDGE_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/Base.h>
#include <IMP/bff/HierarchyFrame.h>
#include <IMP/bff/StructureIO.h>
#include <IMP/bff/SelectionExpression.h>
#include <IMP/bff/StripMask.h>
#include <IMP/bff/ProbeSampling.h>
#include <IMP/bff/VdwRadii.h>
#include <IMP/bff/PathMap.h>

#include <IMP/Model.h>
#include <IMP/Particle.h>
#include <IMP/atom/Atom.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Selection.h>
#include <IMP/atom/bond_decorators.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

// -------- from HierarchyFrame.h --------
//! Coordinates and residue index of every XYZ leaf, in hierarchy order.
/*!
    \param[in] hierarchy the frame to read
    \return four values per atom: x, y, z, and the residue index as a double
            (-1 where the leaf has no residue parent)
*/
IMPBFFEXPORT void hierarchy_atom_coordinates(
        IMP::atom::Hierarchy hierarchy,
        double** out_view, int* n_out_view);

IMPBFFEXPORT std::string atom_name(IMP::atom::Atom atom);

/*!
    The atom name is the last whitespace-separated field of IMP's particle
    name: IMP names an atom `"Atom CB"`, and the caller wants `CB`.
*/
IMPBFFEXPORT std::vector<std::string> hierarchy_atom_metadata(
        IMP::atom::Hierarchy hierarchy);

IMPBFFEXPORT ProteinFrame protein_frame_from_hierarchy(
        IMP::atom::Hierarchy hierarchy);


// -------- from StructureIO.h --------
//! A PDB read into a model: the hierarchy, its leaves and their coordinates.
class IMPBFFEXPORT LoadedStructure {
    IMP::atom::Hierarchy hierarchy_;
    IMP::ParticlesTemp leaves_;
    std::vector<double> coords_;

 public:
    LoadedStructure() {}
    LoadedStructure(IMP::atom::Hierarchy h, const IMP::ParticlesTemp& leaves,
                    const std::vector<double>& coords)
        : hierarchy_(h), leaves_(leaves), coords_(coords) {}

    //! By value: IMP's value types may not be handed out by non-const ref.
    IMP::atom::Hierarchy get_hierarchy() const { return hierarchy_; }
    //! The leaf particles, in hierarchy order -- three coordinates each.
    IMP::ParticlesTemp get_leaves() const { return leaves_; }
    int get_number_of_leaves() const {
        return static_cast<int>(leaves_.size());
    }
    //! `(n_leaves, 3)` in Angstrom.
    void get_coords(double** output, int* n_output1, int* n_output2) const;
};

//! What a FlexFit block names: flexible residues and the bonds to add.
class IMPBFFEXPORT FlexFitSelection {
    IMP::ParticlesTemp residues_;
    IMP::atom::Bonds bonds_;

 public:
    FlexFitSelection() {}
    FlexFitSelection(const IMP::ParticlesTemp& residues,
                     const IMP::atom::Bonds& bonds)
        : residues_(residues), bonds_(bonds) {}

    //! One particle per flexible residue, in the order the block lists them.
    IMP::ParticlesTemp get_residues() const { return residues_; }
    //! The bonds the block declares between named atoms, created in the model.
    IMP::atom::Bonds get_bonds() const { return bonds_; }
};

//! Read a PDB into \p m and return its hierarchy.
/*! `NonWaterPDBSelector`, which is what every caller of this wanted. */
IMPBFFEXPORT IMP::atom::Hierarchy read_pdb_hierarchy(const std::string& path,
                                                     IMP::Model* m);

//! `(N, 3)` coordinates of a hierarchy's XYZ leaves, in hierarchy order.
IMPBFFEXPORT void structure_coordinates(IMP::atom::Hierarchy hierarchy,
                                        double** out_view, int* n_out_view);

IMPBFFEXPORT LoadedStructure load_structure_with_particles(
        const std::string& path, IMP::Model* model);

IMPBFFEXPORT FlexFitSelection read_angle_file(IMP::atom::Hierarchy hier,
                                              const std::string& flexfit_json);

IMPBFFEXPORT IMP::ParticleIndexes select_flexible_residues(
        IMP::atom::Hierarchy hierarchy,
        const std::vector<AtomReference>& residues);

//! Create the single bonds a flexfit file names between two named atoms.
/*! \param[in] atom_pairs two entries per bond, in order */
IMPBFFEXPORT IMP::ParticleIndexes create_named_bonds(
        IMP::atom::Hierarchy hierarchy,
        const std::vector<AtomReference>& atom_pairs);


// -------- from SelectionExpression.h --------
IMPBFFEXPORT IMP::atom::Selection selection_from_expression(
        IMP::atom::Hierarchy hierarchy, const std::string& expression);

//! The atoms of \p hierarchy, in hierarchy order, as a selection sees them.
IMPBFFEXPORT std::vector<SelectionAtom> selection_atoms(
        IMP::atom::Hierarchy hierarchy);

IMPBFFEXPORT std::vector<int> select_atom_indices(
        IMP::atom::Hierarchy hierarchy, const std::string& expression);


// -------- from StripMask.h --------
IMPBFFEXPORT std::vector<double> strip_obstacles(
        IMP::atom::Hierarchy hierarchy, const std::string& mask,
        const std::string& keep = "");

IMPBFFEXPORT StripReport strip_report(IMP::atom::Hierarchy hierarchy,
                                      const std::string& mask);


// -------- from ProbeSampling.h --------
/*!
    The hierarchy's leaves must be the dye's atoms and \p coords three values
    per leaf. \throw ValueException on an atom-count mismatch.
*/
IMPBFFEXPORT void apply_coordinates(const IMP::atom::Hierarchy hierarchy,
                                            const std::vector<double>& coords);


// -------- from VdwRadii.h --------
IMPBFFEXPORT double olga_vdw_particle_radius(IMP::Particle *p);

IMPBFFEXPORT std::vector<std::string> olga_vdw_unknown_atom_names(
        const IMP::ParticlesTemp &ps);

//! Olga's radius for every particle of `ps`, in order (Angstrom).
IMPBFFEXPORT std::vector<double> olga_vdw_radii(const IMP::ParticlesTemp &ps);

// -------- from PathMap.h --------
//! Take a path map's obstacle spheres from particles, and keep taking them.
/*! Installs a sphere source on \p map that reads the XYZR of \p ps at every
    sample_obstacles(), so the map samples the particles where they are now.
    This is what PathMap::set_particles used to be; the lattice itself no
    longer knows a particle. */
IMPBFFEXPORT void set_path_map_particles(PathMap* map, const IMP::ParticlesTemp& ps);

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_HIERARCHYBRIDGE_H
