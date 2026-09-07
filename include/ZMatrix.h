/**
 *  \file IMP/bff/ZMatrix.h
 *  \brief A molecule's internal coordinates: measure them, and rebuild from
 *         them.
 *
 *  The Z-matrix every chemistry toolkit knows: each atom below a base triple
 *  is placed from three already-placed atoms by a bond length, a bond angle
 *  and one dihedral. One dihedral per atom is the atom's only degree of
 *  freedom, which is what makes the representation rotamer-shaped -- a
 *  conformer *is* its dihedral vector, and any dihedral vector is a conformer
 *  (continuous sampling, not just stored bins).
 *
 *  This is general over the molecule: an amino-acid side chain, a dye+linker
 *  (the `.drot` store of PRD-118), any connected graph of atoms. Bond
 *  connectivity can be perceived from coordinates and elements (validated
 *  against a MOL2 bond table on Alexa488+linker: identical 87 bonds) or
 *  supplied by the caller.
 *
 *  The placement primitive is FASPR's ``Internal2Cartesian`` construction
 *  (see \c src/faspr/Utility.cpp) in double precision and the standard
 *  (IUPAC) dihedral convention -- the exact inverse of \c dihedral_deg, to
 *  ~1e-13, so encode/decode round-trips are lossless. The convention itself
 *  is pinned by test: a dihedral that round-trips against its own placement
 *  proves nothing about the world's convention: standard+180 round-trips
 *  just as cleanly, and only an independent A/B catches it (okf/log.md
 *  2026-08-22).
 *
 *  \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 *  \acknowledgment The construction follows FASPR: Huang X, Pearce R, Zhang Y.
 *  "FASPR: an open-source tool for fast and accurate protein side-chain
 *  packing." Bioinformatics 2020;36:3758-3765 (MIT;
 *  https://github.com/tommyhuangthu/FASPR).
 */

#ifndef IMPBFF_ZMATRIX_H
#define IMPBFF_ZMATRIX_H

#include <IMP/bff/bff_config.h>

#include <IMP/algebra/Vector3D.h>

#include <string>
#include <utility>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! The signed dihedral of (a-b-c-d) in degrees, standard (IUPAC) convention.
/** The praxeolitic form: project the outer bonds onto the plane perpendicular
 *  to the central axis b-c and take the signed angle between the projections.
 *  Matches IMP.atom.get_dihedral's convention and FASPR's Dihedral. */
IMPBFFEXPORT double dihedral_deg(const algebra::Vector3D& a,
                                 const algebra::Vector3D& b,
                                 const algebra::Vector3D& c,
                                 const algebra::Vector3D& d);

//! The angle at `b` in the triple (a, b, c), radians.
/** When an arm has zero length the angle is undefined and this returns the
    **tetrahedral** angle, 1.9106 rad. That is a choice, made because a
    harmonic whose minimum sits at a collapsed angle is not a default anyone
    wants; the alternative spelling returned 0.0. Neither fires on real input
    -- 6,364 angles across the shipped structures have no zero-length arm. */
IMPBFFEXPORT double bond_angle_rad(const algebra::Vector3D& a,
                                   const algebra::Vector3D& b,
                                   const algebra::Vector3D& c);

//! The same angle in degrees. One formula, two units, each in its name.
IMPBFFEXPORT double bond_angle_deg(const algebra::Vector3D& a,
                                   const algebra::Vector3D& b,
                                   const algebra::Vector3D& c);

//! Place one atom from internal coordinates (FASPR's construction).
/** The new atom bonds to `parent` with length `bond_length`; `bond_angle_deg`
 *  is the angle at `parent` between `gp` and the new atom; `phi_deg` is the
 *  signed dihedral (ggp, gp, parent, new) as \c dihedral_deg measures it.
 *  Exact inverse of measuring those three values back, to ~1e-13. */
IMPBFFEXPORT algebra::Vector3D internal2cartesian(
        const algebra::Vector3D& ggp, const algebra::Vector3D& gp,
        const algebra::Vector3D& parent,
        double bond_length, double bond_angle_deg_, double phi_deg);

//! The covalent radius of an element, Angstrom; carbon's for an unknown one.
/** One table for the whole module. `StructureIO`'s bond inference had a
    second one (Cordero radii under a multiplicative cutoff) that agreed with
    this one on every molecule in the tree but not on the cutoffs themselves
    -- C-H at 1.49 A against 1.31 A -- so a stretched bond fell on different
    sides depending on which entry point was asked. */
IMPBFFEXPORT double covalent_radius(const std::string& element);

//! The tolerance #perceive_bonds adds to the sum of two covalent radii.
IMPBFFEXPORT double bond_tolerance();

//! Covalent bonds perceived from geometry: pairs (i, j), i < j.
/** Bonded when |xi - xj| < r_i + r_j + 0.35 A over the covalent radii of
 *  H, C, N, O, P, S. Validated against the MOL2 bond table of the shipped
 *  Alexa488 C1R template: the identical 87 bonds. `elements` holds element
 *  symbols parallel to `xyz`. */
IMPBFFEXPORT std::vector<std::pair<int, int> > perceive_bonds(
        const std::vector<algebra::Vector3D>& xyz,
        const std::vector<std::string>& elements);

//! A molecule's spanning-tree internal coordinates, with encode and decode.
/** Built once from a template (coordinates, and either elements -- bonds are
 *  perceived -- or an explicit bond list); used many times:
 *
 *  - \c encode(coords) -> one dihedral per row, the conformer's chi vector;
 *  - \c decode(chi) -> full coordinates, template bond lengths/angles and
 *    the given dihedrals (encode(decode(chi)) == chi exactly);
 *  - \c decode(chi, base) seeds the base atoms from a caller's frame.
 *
 *  Rows are (atom, parent, grandparent, great-grandparent) in breadth-first
 *  order from `root`, so every atom's three references are placed before it.
 *  The base is the root plus the atoms too close to it to have a complete
 *  reference triple. Row count is n_atoms - n_base; the dihedral a row
 *  controls is the IUPAC dihedral of its four atoms.
 *
 *  \note The representation is exact for any connected molecule, but it is
 *  torsion-parametric: decode rebuilds with the *template's* bond lengths
 *  and angles. A conformer whose bonds stretch (real MD does, by ~0.004 A)
 *  encodes and decodes to the torsion-equivalent rigid-geometry conformer.
 *  For a lossless store of such ensembles keep per-conformer angles too
 *  (the `.drot` v2 layout, PRD-118). */
class IMPBFFEXPORT ZMatrix {
public:
    ZMatrix() {}

    //! Build from a template; bonds perceived from `elements`.
    /** `root` picks the spanning root (default 0). The template is stored
     *  and defines the reference geometry of every decode. */
    void set_template(const std::vector<algebra::Vector3D>& xyz,
                      const std::vector<std::string>& elements,
                      int root = -1);

    //! Build from a template with an explicit bond list (pairs i < j).
    void set_template_bonds(const std::vector<algebra::Vector3D>& xyz,
                            const std::vector<std::pair<int, int> >& bonds,
                            int root = -1);

    //! Number of atoms of the template.
    unsigned int get_number_of_atoms() const { return (unsigned int) n_; }

    //! The base atoms (root + atoms without a complete reference triple).
    std::vector<int> get_base() const { return base_; }

    //! All rows, flattened: (atom, parent, gp, ggp) per row.
    std::vector<int> get_rows() const { return rows_flat_; }

    //! The template's own dihedral vector (what decode() with no argument
    //! rebuilds).
    std::vector<double> get_template_chi() const { return chi_template_; }

    //! The dihedral vector of a conformer, degrees, one entry per row.
    std::vector<double> encode(const std::vector<algebra::Vector3D>& xyz)
            const;

    //! Full coordinates of the conformer with dihedral vector `chi`.
    /** Bond lengths and angles are the template's; `chi` must hold one value
     *  per row. Base atoms keep their template positions unless `base`
     *  (full length, parallel to the template) is given, whose base-atom
     *  rows seed the frame instead. */
    std::vector<algebra::Vector3D> decode(
            const std::vector<double>& chi,
            const std::vector<algebra::Vector3D>* base = nullptr) const;

    //! Per-row bond lengths and angles of the template.
    std::vector<double> get_template_lengths() const { return lens_; }
    std::vector<double> get_template_angles() const { return angs_; }

    //! The (length, angle, dihedral) a conformer actually has, per row.
    /** Where \c encode measures only the dihedral and lets \c decode supply
     *  the template's rigid geometry, this measures all three degrees of
     *  freedom -- so the conformer can be rebuilt exactly, bond stretching
     *  and angle bending included. Real MD flexes both (~0.004 A, ~1 deg),
     *  and dropping them is the single largest error in a torsion-only
     *  store: 0.016 A RMSD on the shipped dye libraries, which is enough to
     *  move a FRET efficiency in the third decimal. This is what the `.drot`
     *  v9 grids hold (PRD-118).
     *
     *  All three vectors are resized to the row count. */
    void frame_internals(const std::vector<algebra::Vector3D>& xyz,
                         std::vector<double>& lengths,
                         std::vector<double>& angles,
                         std::vector<double>& chi) const;

    //! Rebuild a conformer from per-row internals, not the template's.
    /** The exact inverse of \c frame_internals: with that function's output
     *  and the frame's own base atoms, the deviation from the original
     *  coordinates is ~1e-13 A. `base` seeds the base atoms (full length,
     *  parallel to the template); without it they keep template positions. */
    std::vector<algebra::Vector3D> decode_internals(
            const std::vector<double>& lengths,
            const std::vector<double>& angles,
            const std::vector<double>& chi,
            const std::vector<algebra::Vector3D>* base = nullptr) const;

    //! The perceived/explicit bonds the tree was built from.
    std::vector<std::pair<int, int> > get_bonds() const { return bonds_; }

private:
    void build_(int root);

    unsigned int n_ = 0;
    std::vector<std::pair<int, int> > bonds_;       //!< i < j
    std::vector<std::vector<int> > adj_;            //!< ascending neighbors
    std::vector<int> base_;
    std::vector<int> rows_flat_;                    //!< 4 per row
    std::vector<double> lens_, angs_;               //!< per row, template
    std::vector<double> chi_template_;              //!< per row
    std::vector<algebra::Vector3D> template_;
};

IMPBFF_END_NAMESPACE

#endif /* IMPBFF_ZMATRIX_H */