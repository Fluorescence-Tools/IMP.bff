/**
 *  \file IMP/bff/RotamerFret.h
 *  \brief FRET over a trajectory from two screened rotamer libraries.
 *
 * The driver FRETpredict is: place a rotamer library at each of two labelled
 * residues in every frame, screen both against the protein, and report the
 * pair's \f$\langle\kappa^2\rangle\f$ and its efficiency in the three
 * averaging limits, frame by frame.
 *
 * **The parameter names here are FRETpredict's**, deliberately -- `fixed_R0`,
 * `ign_H`, `libname_1`, `r0lib` -- because the parity harness
 * (`test/cgprobe/rotamer/expensive_test_fretpredict_parity.py`) hands the *same*
 * keyword dictionary to this class and to FRETpredict and compares the files
 * they write. A rename here would be a rename of the experiment.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_ROTAMERFRET_H
#define IMPBFF_ROTAMERFRET_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/RotamerEnsemble.h>
#include <IMP/bff/RotamerLibrary.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! What one frame of a trajectory says about a labelled pair.
struct IMPBFFEXPORT FRETFrameResult {
    //! The two partition functions, donor first: how much of each library the
    //! site leaves accessible.
    double z_donor, z_acceptor;
    double kappa2_avg;
    //! The three averaging limits (see #IMP::bff::FRETPairEfficiencies).
    double static_efficiency, dynamic1, dynamic2;

    FRETFrameResult()
        : z_donor(0), z_acceptor(0), kappa2_avg(2.0 / 3.0),
          static_efficiency(0), dynamic1(0), dynamic2(0) {}

    IMP_SHOWABLE_INLINE(FRETFrameResult,
                        out << "FRETFrameResult(E=" << static_efficiency
                            << ", k2=" << kappa2_avg << ")");
};
IMP_VALUES(FRETFrameResult, FRETFrameResults);

//! FRET over a trajectory from two screened rotamer libraries.
class IMPBFFEXPORT RotamerFRET {
    std::vector<ProteinFrame> frames_;
    std::vector<int> residues_;
    std::vector<std::string> chains_;
    std::string donor_, acceptor_, libname_1_, libname_2_, r0lib_;
    std::string output_prefix_;
    RotamerSiteOptions site_;
    RotamerLibrary lib_1_, lib_2_;
    bool fixed_r0_;
    double r0_, z_cutoff_;
    std::vector<double> user_weights_, weights_;

    std::vector<double> z_values_, k2_values_, estatic_, edynamic1_,
            edynamic2_;

    void load_libraries();
    void write_summary(const std::string& prefix,
                       const std::vector<double>& k2,
                       const std::vector<double>& es,
                       const std::vector<double>& ed1,
                       const std::vector<double>& ed2,
                       const std::vector<double>& weights) const;

    //! The frames-in constructor #from_frames uses.
    RotamerFRET(const std::vector<ProteinFrame>& frames,
                const std::vector<int>& residues,
                const std::vector<std::string>& chains,
                const std::string& donor, const std::string& acceptor,
                const std::string& libname_1, const std::string& libname_2,
                const RotamerSiteOptions& site, bool fixed_R0, double r0,
                const std::string& r0lib, double z_cutoff,
                const std::string& output_prefix,
                const std::vector<double>& user_weights);

public:
    //! \param[in] protein a PDB (or multi-MODEL PDB) path
    /*! \param[in] residues the two labelled residue numbers
        \param[in] chains their chains; empty matches any
        \param[in] donor,acceptor the dye names the spectra are looked up by
        \param[in] libname_1,libname_2 the two rotamer libraries
        \param[in] temperature K
        \param[in] electrostatic add the Debye-Hueckel term
        \param[in] potential `lj` or `gauss`
        \param[in] ign_H mask hydrogens out of the steric term
        \param[in] sigma_scaling,epsilon_scaling dye-probe LJ scaling
        \param[in] fixed_R0 keep \p r0 rather than computing one per frame
        \param[in] r0 the Förster radius, **Angström** (it was nanometres
                   until 2026-08-25, when R0 became one unit across the
                   package: a fixed 5.5 nm is 55 here)
        \param[in] r0lib a dye-library CIF to take the spectra from
        \param[in] z_cutoff frames where either partition function falls below
                   this report NaN rather than a number the screening does not
                   support
        \param[in] output_prefix what #save writes its files under
        \param[in] max_frames stop after this many; negative reads all
        \param[in] user_weights one per frame, for the summary table
        \throw ValueException unless exactly two residues are given */
    RotamerFRET(const std::string& protein,
                const std::vector<int>& residues = std::vector<int>(),
                const std::vector<std::string>& chains =
                        std::vector<std::string>(),
                const std::string& donor = "AlexaFluor 488",
                const std::string& acceptor = "AlexaFluor 594",
                const std::string& libname_1 = "AlexaFluor 488 C1R cutoff30",
                const std::string& libname_2 = "AlexaFluor 594 C1R cutoff30",
                double temperature = 300.0, bool electrostatic = false,
                const std::string& potential = "lj", bool ign_H = true,
                double sigma_scaling = 0.5, double epsilon_scaling = 1.0,
                bool fixed_R0 = false, double r0 = 54.0,
                const std::string& r0lib = "", double z_cutoff = 0.05,
                const std::string& output_prefix = "res",
                int max_frames = -1,
                const std::vector<double>& user_weights =
                        std::vector<double>());

    //! The same, over frames a caller already has.
    /*! The door for the formats this module does not read itself: an RMF
        trajectory becomes `ProteinFrame` values (`protein_frames_from_rmf`)
        and comes in here.

        A named factory and not a second constructor, because SWIG turns
        keyword arguments off for anything overloaded and the constructor's
        twenty parameters are keyword arguments in every caller. */
    static RotamerFRET from_frames(
            const std::vector<ProteinFrame>& frames,
                const std::vector<int>& residues,
                const std::vector<std::string>& chains =
                        std::vector<std::string>(),
                const std::string& donor = "AlexaFluor 488",
                const std::string& acceptor = "AlexaFluor 594",
                const std::string& libname_1 = "AlexaFluor 488 C1R cutoff30",
                const std::string& libname_2 = "AlexaFluor 594 C1R cutoff30",
                double temperature = 300.0, bool electrostatic = false,
                const std::string& potential = "lj", bool ign_H = true,
                double sigma_scaling = 0.5, double epsilon_scaling = 1.0,
                bool fixed_R0 = false, double r0 = 54.0,
                const std::string& r0lib = "", double z_cutoff = 0.05,
            const std::string& output_prefix = "res",
            const std::vector<double>& user_weights = std::vector<double>());

    //! What one frame says: both partition functions, ⟨κ²⟩ and the three
    //! efficiencies.
    /*! \note With `fixed_R0` false this **sets** the object's `r0` to the one
        this frame's \f$\langle\kappa^2\rangle\f$ gives, which is what
        FRETpredict does and why `r0` is readable afterwards. */
    FRETFrameResult frame_fret(const ProteinFrame& frame);

    //! Every frame, into the value arrays.
    /*! A frame where either partition function is at or below `z_cutoff`
        keeps its partition functions and reports NaN for the rest: the
        screening left the dye nowhere to be, and a number computed from that
        would be a number about the fallback, not about the structure. */
    void trajectory_analysis();

    //! Write the FRETpredict-style output files.
    /*!
        `<prefix>-Z-<r1>-<r2>.dat` and one each for `w_s`, `k2`, `Es`, `Ed1`,
        `Ed2`, plus a labelled `-data-` table of (average, SD, SE) per
        quantity. A text table and not a pickled DataFrame: a `.pkl` cannot be
        read without the library that wrote it.

        \param[in] output_prefix overrides the constructor's
        \throw ValueException when `user_weights` is not one per frame
    */
    void save(const std::string& output_prefix = "");

    //! Recompute the summary table from the saved files with other weights.
    /*!
        \param[in] boltzmann_weights weight each frame by the product of its
                   two partition functions (read back from the `-Z-` file)
        \param[in] user_weights one per frame, multiplied into the above
        \param[in] output_prefix where the table goes; the inputs are always
                   read from the constructor's prefix
    */
    void reweight(bool boltzmann_weights = false,
                  const std::vector<double>& user_weights =
                          std::vector<double>(),
                  const std::string& output_prefix = "");

    //! #trajectory_analysis followed by #save.
    void run();

    //! `(n_frames, 2)` partition functions, donor first.
    void get_z_values(double** out_view, int* n_out_view) const;
    //! One per frame; NaN where the screening left nothing.
    void get_k2_values(double** out_view, int* n_out_view) const;
    void get_estatic_values(double** out_view, int* n_out_view) const;
    void get_edynamic1_values(double** out_view, int* n_out_view) const;
    void get_edynamic2_values(double** out_view, int* n_out_view) const;

    //! The Förster radius in Å -- the fixed one, or the last frame's.
    double get_r0() const { return r0_; }
    //! The two labelled residues, and the chains they are on.
    std::vector<int> get_residues() const { return residues_; }
    std::vector<std::string> get_chains() const { return chains_; }
    //! The two dyes, and the two libraries screened at their sites.
    std::string get_donor() const { return donor_; }
    std::string get_acceptor() const { return acceptor_; }
    std::string get_libname_1() const { return libname_1_; }
    std::string get_libname_2() const { return libname_2_; }
    int get_n_frames() const { return static_cast<int>(frames_.size()); }
    std::string get_output_prefix() const { return output_prefix_; }

    IMP_SHOWABLE_INLINE(RotamerFRET,
                        out << "RotamerFRET(" << frames_.size() << " frames, "
                            << donor_ << " -> " << acceptor_ << ")");
};

//! Build a #IMP::bff::RotamerFRET from an fps.json distance.
/*! The file names the two positions, their dyes and their libraries; anything
    it does not name is left at the constructor's default.

    \param[in] fps_path the fps.json file
    \param[in] protein the structure its positions refer to
    \param[in] distance_name which distance; empty takes the first
    \param[in] output_prefix what #IMP::bff::RotamerFRET::save writes under
    \param[in] fixed_R0,r0 as in the constructor (\p r0 in Å)
    \param[in] temperature,electrostatic as in the constructor
*/
IMPBFFEXPORT RotamerFRET rotamer_fret_from_fps(
        const std::string& fps_path, const std::string& protein,
        const std::string& distance_name = "",
        const std::string& output_prefix = "res", bool fixed_R0 = false,
        double r0 = 54.0, double temperature = 300.0,
        bool electrostatic = false);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_ROTAMERFRET_H
