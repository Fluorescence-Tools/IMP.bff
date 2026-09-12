#ifndef IMPBFF_LABELIZER_FEATURES_H
#define IMPBFF_LABELIZER_FEATURES_H

/**
 *  \file IMP/bff/LabelizerFeatures.h
 *  \brief The structural quantities a label-site score is computed from:
 *         secondary structure, half-sphere exposure, relative solvent
 *         accessibility, residue depth, and an imported conservation grade.
 *
 * These were the parts of the Labelizer package (Gebhardt *et al.*,
 * *Nat. Commun.* **16**, 3305, 2025) that did not run: `ss` shelled out to
 * DSSP and refused to start on macOS at all
 * (`secondary_structure.py:126`, `RuntimeError("Unknown platform")`), while
 * `se` and `me` shelled out to MSMS, whose bundled binaries are 32-bit
 * ppc/i386 Mach-O. A score model whose default solvent-exposure term cannot be
 * evaluated is not a model anyone can check, so the two external programs are
 * replaced here by native kernels and the difference is *measured* rather than
 * asserted — see `okf/validation/labelizer_ab.md`.
 *
 * Nothing here re-implements a kernel the module already has: the accessible
 * surface comes from #IMP::bff::solvent_accessible_surface_area, the structure
 * from #IMP::bff::read_pdb_records.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/bff_config.h>

#include <IMP/bff/Base.h>
#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

// ---------------------------------------------------------------------------
// The residue view
// ---------------------------------------------------------------------------

//! One residue, and where its atoms are in the flat arrays of #LabelizerStructure.
/*! The backbone indices are resolved once, at read time, because every kernel
    below wants them and each would otherwise re-scan the atom list. A missing
    atom is `-1` and stays `-1`: a residue without a CA cannot be scored, and
    substituting a neighbouring atom yields a plausible, wrong answer. */
struct IMPBFFEXPORT LabelizerResidue {
    //! Author chain identifier, as the PDB writes it.
    std::string chain;
    //! Author residue number.
    int seq_id;
    //! Three-letter residue name, e.g. `"TRP"`.
    std::string comp_id;
    //! Indices into #LabelizerStructure::xyz (as triples) and the parallel arrays.
    std::vector<int> atoms;
    //! Backbone atom indices, or -1 when the atom is absent.
    int n, ca, c, o, cb;

    LabelizerResidue() : seq_id(0), n(-1), ca(-1), c(-1), o(-1), cb(-1) {}
    IMP_SHOWABLE_INLINE(LabelizerResidue,
                        out << "LabelizerResidue(" << chain << seq_id << ":"
                            << comp_id << ")");
};
IMP_VALUES(LabelizerResidue, LabelizerResidues);

//! A structure as the scoring layer sees it: flat coordinates, grouped residues.
struct IMPBFFEXPORT LabelizerStructure {
    //! `3N` coordinates, xyzxyz…
    std::vector<double> xyz;
    //! `N` van der Waals radii.
    std::vector<double> vdw;
    //! `N` atom names, stripped of padding.
    std::vector<std::string> atom_name;
    //! Residues in file order.
    std::vector<LabelizerResidue> residues;

    //! Number of atoms.
    std::size_t size() const { return vdw.size(); }
};

//! The twenty standard amino acids, three-letter, uppercase.
IMPBFFEXPORT const std::vector<std::string>& labelizer_standard_residues();

//! Read a structure and group it into residues.
/*!
    Reads through #IMP::bff::read_pdb_records, so the parse and its cache are
    shared with the AV builder.

    \param[in] pdb_path the structure
    \param[in] protein_only keep only the residues named by
               #labelizer_standard_residues. This is how the reference's
               `remove_hetatoms` (`pdbhelper.py:16`, which detaches every
               residue whose hetflag is not blank) is reproduced: the record
               type is not carried on a #IMP::bff::PDBAtomRecord, and for a
               protein the two tests differ only on modified residues such as
               MSE, which the reference also drops.
    \param[in] model for a multi-model file, which model to take. The reference
               takes `structure[0]` unconditionally (`pdbhelper.py:35`).
    \return the grouped structure
    \throw IOException when the file has no coordinates
*/
IMPBFFEXPORT LabelizerStructure labelizer_read_structure(const std::string& pdb_path,
                                           bool protein_only = true,
                                           int model = 0);

// ---------------------------------------------------------------------------
// Geometry the reference needed and biopython supplied
// ---------------------------------------------------------------------------

//! The Cbeta position of a residue, real or reconstructed for glycine.
/*!
    A glycine has no Cbeta, and the reference builds a virtual one rather than
    skipping the residue: the N vector is centred on CA, rotated by -120
    degrees about the CA->C axis, and added back to CA
    (`fret_score.py:293`, itself lifted from `Bio/PDB/HSExposure.py`).

    \param[in] s the structure
    \param[in] residue index into `s.residues`
    \param[out] out_view,n_out_view three coordinates, or zero-length when the
                residue lacks the backbone atoms to place one
*/
IMPBFFEXPORT void labelizer_cbeta_position(const LabelizerStructure& s, int residue,
                                    double** out_view, int* n_out_view);

//! Half-sphere exposure about Cbeta, the `(up, down)` neighbour counts.
/*!
    For each residue the CA-CB direction splits a sphere of the given radius in
    two; `up` counts the CA atoms of other residues on the side the side chain
    points to, `down` the rest. It is a cheap burial measure that needs no
    surface, which is why the reference offers it as an alternative to MSMS
    (`solvent_exposure.py:141`, `Bio.PDB.HSExposure.HSExposureCB`).

    \param[in] s the structure
    \param[in] radius sphere radius in Angstrom; the reference uses 13.0 for
               the tryptophan and charge terms and 10/13/16 for the `I_SE*`
               solvent-exposure tables
    \param[out] out_view,n_out_view `2 * n_residues` values, up then down for
                each residue in `s.residues` order
*/
IMPBFFEXPORT void labelizer_half_sphere_exposure(const LabelizerStructure& s, double radius,
                                          double** out_view, int* n_out_view);

// ---------------------------------------------------------------------------
// Secondary structure
// ---------------------------------------------------------------------------

//! DSSP secondary structure, one 8-state code per residue.
/*!
    Kabsch and Sander (1983), implemented natively. The electrostatic hydrogen
    bond energy between the C=O of residue `i` and the N-H of residue `j` is

        E = 332 * 0.42 * 0.20 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN)   kcal/mol

    with a bond declared below #LABELIZER_DSSP_HBOND_ENERGY. The amide hydrogen is not
    in the file and is placed the way DSSP places it: on N, along the direction
    opposite the preceding residue's C=O.

    The eight states are `H` (alpha helix), `B` (isolated bridge), `E` (strand),
    `G` (3-10 helix), `I` (pi helix), `T` (turn), `S` (bend) and `-` (none),
    assigned in that order of precedence.

    \param[in] s the structure
    \return one character per residue, in `s.residues` order
*/
IMPBFFEXPORT std::string labelizer_dssp(const LabelizerStructure& s);

//! The hydrogen-bond energy below which a bond is declared, kcal/mol.
IMPBFFEXPORT extern const double LABELIZER_DSSP_HBOND_ENERGY;

// ---------------------------------------------------------------------------
// Solvent exposure
// ---------------------------------------------------------------------------

//! Which maximum-accessibility scale a relative accessibility is taken against.
enum LabelizerMaxAsa {
    //! Tien *et al.* (2013) theoretical maxima. The reference's `N_SE1`.
    LABELIZER_MAXASA_WILKE = 0,
    //! Kabsch and Sander (1983). The reference's `N_SE2`.
    LABELIZER_MAXASA_SANDER = 1,
    //! Miller *et al.* (1987). The reference's `N_SE3`.
    LABELIZER_MAXASA_MILLER = 2
};

//! The maximum solvent-accessible area per residue type, Angstrom squared.
IMPBFFEXPORT std::map<std::string, double> labelizer_max_asa(LabelizerMaxAsa scale);

//! Solvent-accessible surface area per residue, Angstrom squared.
/*!
    A Shrake-Rupley sum over the residue's atoms, through the module's existing
    #IMP::bff::solvent_accessible_surface_area.

    \param[in] s the structure
    \param[in] probe_radius rolling-probe radius; 1.4 is the DSSP convention
    \param[in] n_sphere_points sampling density per atom
    \param[out] out_view,n_out_view one area per residue
*/
IMPBFFEXPORT void labelizer_residue_sasa(const LabelizerStructure& s,
                                  double probe_radius, int n_sphere_points,
                                  double** out_view, int* n_out_view);

//! Relative solvent accessibility per residue, a fraction.
/*!
    #labelizer_residue_sasa divided by the residue's maximum. A residue type absent
    from the scale yields a negative value, which the caller must treat as
    "not computed" rather than as zero exposure.

    \param[in] s the structure
    \param[in] scale which maximum-accessibility table
    \param[in] probe_radius,n_sphere_points passed to #labelizer_residue_sasa
    \param[out] out_view,n_out_view one fraction per residue
*/
IMPBFFEXPORT void labelizer_relative_solvent_accessibility(
        const LabelizerStructure& s, LabelizerMaxAsa scale, double probe_radius,
        int n_sphere_points, double** out_view, int* n_out_view);

//! Mean distance from a residue's atoms to the solvent-excluded surface.
/*!
    The reference's default solvent-exposure term is MSMS residue depth
    (`N_SE11_MEAN_SURFACE_DIST`, `solvent_exposure.py:157`), and MSMS does not
    run on any machine this package is developed on. This computes the same
    quantity from a numerically-built solvent-excluded surface: the probe
    centre is rolled over the atoms on the golden-spiral point set
    #IMP::bff::sphere_points, the points that no atom occludes are kept, and
    each is pulled back onto the surface by one probe radius.

    It is *not* MSMS and is not claimed to be. What it is measured against is
    the reference's own published output for 1DDB; the agreement is recorded in
    `okf/validation/labelizer_ab.md`.

    \param[in] s the structure
    \param[in] probe_radius rolling-probe radius, Angstrom
    \param[in] n_sphere_points points per atom used to build the surface
    \param[out] out_view,n_out_view one depth per residue, Angstrom
*/
IMPBFFEXPORT void labelizer_residue_depth(const LabelizerStructure& s, double probe_radius,
                                   int n_sphere_points, double** out_view,
                                   int* n_out_view);

//! Distance from each residue's Cbeta to the solvent-excluded surface.
/*!
    The observable of the `N_SE10_CB_SURFACE_DIST` table, and a different
    quantity from #labelizer_residue_depth: that averages over every atom of the
    residue, this measures the one atom a label is attached at. The bins differ
    accordingly -- the Cbeta table starts at 0.69 Angstrom where the mean-atom
    table starts at 1.46.

    Glycine has no Cbeta and gets the reconstructed one from
    #labelizer_cbeta_position, so every residue has a value.

    \param[in] s the structure
    \param[in] probe_radius,n_sphere_points as #labelizer_residue_depth
    \param[out] out_view,n_out_view one distance per residue, Angstrom;
                negative where no Cbeta could be placed
*/
IMPBFFEXPORT void labelizer_cbeta_depth(const LabelizerStructure& s, double probe_radius,
                                 int n_sphere_points, double** out_view,
                                 int* n_out_view);

//! The whole protein's formal charge: net, positive and negative.
/*!
    A charged dye is not indifferent to the charge of what it is attached to,
    and the reference computes this (`charge_environment.py:160`) — though
    nothing there consumes it, so it is offered rather than used.

    Formal charges only, from the same table the `ce` term counts with:
    ARG, LYS and HIS `+1`, ASP and GLU `-1`, everything else zero. Histidine
    counted as fully charged is the reference's choice and is wrong at
    physiological pH, where it is mostly neutral; it is reproduced because the
    charge environment term depends on it and the two must agree.

    The reference also offers a variant reading a per-atom partial charge out
    of the PDB **occupancy** column, for a PQR file loaded as a PDB. That is
    not ported: `IMP::bff::PDBAtomRecord` does not carry occupancy, and reading
    a charge out of a field that means "fraction of the site this atom
    occupies" is the same category of abuse as the score-in-the-B-factor this
    package's container exists to end. A caller wanting partial charges should
    read a PQR as a PQR.

    \param[in] s the structure
    \param[out] out_view,n_out_view three values: net, positive, negative.
                The negative total is reported as a negative number, so the
                three sum as `net = positive + negative`.
*/
IMPBFFEXPORT void labelizer_global_charge(const LabelizerStructure& s, double** out_view,
                                   int* n_out_view);

// ---------------------------------------------------------------------------
// Conservation, imported
// ---------------------------------------------------------------------------

//! A conservation grade per residue, keyed `"<chain><seq_id>"`.
/*!
    Nothing here computes conservation: an alignment and a rate estimate are a
    different program, and the reference imports them too
    (`conservation_score.py`, whose ConSurf call lives in the web backend, not
    in the package). Both shapes the reference reads are supported.

    \param[in] path a ConSurf `.grades` table, or a PDB carrying the normalised
               grade in the B-factor column
    \return grade by residue key; an empty map when the file holds none
    \throw IOException when the file cannot be read
*/
IMPBFFEXPORT std::map<std::string, double> labelizer_read_consurf(
        const std::string& path);

//! The residue key the score tables and the reference CSVs are indexed by.
/*! `"<chain><seq_id>"`, e.g. `"A123"` — the reference's own convention
    (`labeling_parameter.py:129`), kept so its output files can be compared
    row for row. */
IMPBFFEXPORT std::string labelizer_residue_key(const std::string& chain, int seq_id);

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_LABELIZER_FEATURES_H
