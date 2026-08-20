/**
 *  \file IMP/bff/ModelPrecision.h
 *  \brief How far a docked body wanders between independent runs.
 *
 * The FPS notion of a model's *precision*: dock the same two bodies from several
 * independent starts, superpose the results on the fixed body, and measure how
 * much each atom of the mobile body moves. That per-atom RMSF is the positional
 * uncertainty, and it is the number a docking result should be quoted with.
 *
 * It is a **precision, not an accuracy**: it says the restraints did not pin the
 * body down further, not that the body is there.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_MODELPRECISION_H
#define IMPBFF_MODELPRECISION_H

#include <IMP/bff/bff_config.h>

#include <IMP/showable_macros.h>
#include <IMP/value_macros.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! What repeated docking says about a model's precision.
struct IMPBFFEXPORT PositionUncertainty {
    //! Models that took part -- those whose atom count matched the first.
    int n_models;
    //! Per-atom RMSF about the mean structure, A.
    std::vector<double> rmsf;
    //! Chain of each atom, in file order.
    std::vector<std::string> chains;
    //! The mean structure, flat, three per atom.
    std::vector<double> mean_coords;
    //! Which atoms were superposed on.
    std::vector<int> fixed_mask;

    PositionUncertainty() : n_models(0) {}

    double get_rmsf_mean() const;
    double get_rmsf_max() const;
    //! The mean RMSF over the atoms *not* superposed on -- the moving body.
    double get_mobile_rmsf_mean() const;

    void get_rmsf(double** out_view, int* n_out_view) const;
    void get_mean_coords(double** out_view, int* n_out_view) const;

    IMP_SHOWABLE_INLINE(PositionUncertainty,
                        out << "PositionUncertainty(" << n_models
                            << " models, <rmsf> = " << get_rmsf_mean() << " A)");
};
IMP_VALUES(PositionUncertainty, PositionUncertainties);

//! Superpose docked models on the fixed body and report per-atom RMSF.
/*!
    \param[in] pdb_paths the best-scoring PDB of each run. They are the same
               structure docked from different starts, so the atom ordering is
               shared; a file whose atom count differs is **skipped** rather
               than aligned by name, because a differing count means it is not
               the same structure and pairing atoms by index would be wrong.
    \param[in] fixed_chains chain identifiers of the reference body. When none of
               them appears, every atom is superposed on -- there is no fixed
               body, and aligning on all of it is the only frame available.
    \return a record whose `n_models` is 0 or 1 when fewer than two files could
            be read; an RMSF over one structure is not a precision, so the rest
            comes back empty rather than as zeros.
*/
IMPBFFEXPORT PositionUncertainty estimate_position_uncertainty(
        const std::vector<std::string>& pdb_paths,
        const std::vector<std::string>& fixed_chains);

//! The distinct chain identifiers of a PDB's ATOM/HETATM records, sorted.
/*! What names the fixed body when the caller has a reference structure rather
    than a list of chains. \throw IOException when the file cannot be read. */
IMPBFFEXPORT std::vector<std::string> pdb_chain_ids(const std::string& path);

//! Write the mean structure with the RMSF in the B-factor column.
/*! Columns outside the coordinates and the B-factor are copied from the first
    input file, so the output is that file with two fields replaced. */
IMPBFFEXPORT void write_position_uncertainty_pdb(
        const PositionUncertainty& uncertainty,
        const std::string& template_pdb, const std::string& out_pdb);

//! Write `atom_index,chain,rmsf`.
IMPBFFEXPORT void write_position_uncertainty_csv(
        const PositionUncertainty& uncertainty, const std::string& out_csv);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_MODELPRECISION_H
