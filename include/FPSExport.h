/**
 *  \file IMP/bff/FPSExport.h
 *  \brief What a docking or screening run is written out as: FPS's six exports.
 *
 * FPS's `SaveForm` writes a PyMOL script per result, an Overlay, an
 * OverlayStates, an R table, a \f$\chi^2\f$ table and a .NET `BinaryFormatter`
 * dump. The first five are here. The sixth is not and will not be: it is
 * MS-NRBF carrying the field layout of a long-dead assembly, and it is
 * meaningless without the project it was written beside
 * (`okf/references/fps-export-formats.md` §8). #IMP::bff::FPSResultTable::get_json
 * is what replaces it.
 *
 * Three of FPS's defects are **not** reproduced, and each is a named decision
 * rather than a silent improvement:
 *
 * 1. **The best fit is computed where it is used.** FPS reads
 *    `SimulationResult.BestFitRotation` in all three script writers and none of
 *    them computes it -- it is a side effect of the GUI's RMSD column, so an
 *    exported "Overlay" is a chain of *pairwise* fits under the default
 *    "RMSD vs previous", and an all-zero matrix (which `AngleAndAxis` turns
 *    into a spurious 180 degree rotation) if the best-fit box happened to be
 *    unticked. See #IMP::bff::add_best_fit and §7 of the reference.
 * 2. **`_tmp.pdb` is written where the export is**, not into PyMOL's current
 *    working directory, and it is deleted at the end of the script.
 * 3. **`rotate` and `translate` carry `camera=0`.** PyMOL interprets both in
 *    *camera* space by default; FPS writes model-frame vectors and gets away
 *    with it only because the camera is identity in a fresh session.
 *
 * And two of its column conventions are kept but named:
 *
 * * the \f$\chi^2\f$ column is **not reduced** (FPS's file is not; its GUI
 *   shows one that is), so the file carries **both**, as `chi2` and `chi2_r`;
 * * `chi2_bond` is a **subset** of `chi2`, never an addition to it.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_FPSEXPORT_H
#define IMPBFF_FPSEXPORT_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/Docking.h>

#include <IMP/algebra/Transformation3D.h>
#include <IMP/algebra/Vector3D.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! One input structure, as FPS's `Molecule` is used by the export.
struct IMPBFFEXPORT FPSMolecule {
    //! `Molecule.FullFileName`: the path a `load` line names.
    std::string path;
    //! `Molecule.Name`: the basename without extension, i.e. the PyMOL object.
    /*! FPS uses it unquoted and unsanitised, so two inputs with the same
        basename from different directories silently become one object and a
        `name_1`. #fps_molecules leaves the name alone but a caller can rename
        a molecule before exporting. */
    std::string name;
    //! `Molecule.CM`: the origin every `rotate` line is taken about.
    /*! Computed the way the assembly computes it -- the rigid body's own
        reference-frame origin after `IMP::atom::create_rigid_body` -- rather
        than by a second formula, so the `origin=` of the script and the frame
        the pose was captured in cannot disagree. */
    IMP::algebra::Vector3D center;
    //! The body's reference frame **as the untransformed PDB gives it**.
    /*!
        Needed, and easy to leave out -- which is a bug, so it is written down.
        `IMP::atom::create_rigid_body` does **not** start a body at the
        identity: it sets the frame from the members' principal axes, so the
        pose #IMP::bff::capture_poses records is not the transform to apply to
        the coordinates a `load` line reads. The transform a PyMOL script needs
        is \f$T_{pose} \circ T_{input}^{-1}\f$, and this is
        \f$T_{input}\f$. Measured on HIV-RT: the body-0 quaternion at the
        input pose is (0.8870, -0.1975, -0.4164, 0.0301), an 55 degree
        rotation, not (1, 0, 0, 0).
    */
    IMP::algebra::Transformation3D frame;
    int n_atoms;

    FPSMolecule()
        : center(0, 0, 0),
          frame(IMP::algebra::get_identity_transformation_3d()), n_atoms(0) {}

    IMP_SHOWABLE_INLINE(FPSMolecule,
                        out << "FPSMolecule(" << name << ", " << n_atoms
                            << " atoms)");
};
IMP_VALUES(FPSMolecule, FPSMolecules);

//! Read the structures an export will `load`, and where their frames sit.
/*! \param[in] pdb_paths the same list `build_docking_assembly` was given
    \throw IOException when a file is missing */
IMPBFFEXPORT FPSMolecules fps_molecules(
        const std::vector<std::string>& pdb_paths);

//! One labelling position, ready to become a PyMOL pseudoatom.
struct IMPBFFEXPORT FPSLabelPosition {
    std::string name;
    //! `Donor`, `Acceptor` or empty for FPS's `DyeType.Unknown`.
    /*! Only the first two get a `color` line, as in FPS. */
    std::string dye;
    //! Which molecule the position rides on; its transform is applied.
    int body_id;
    //! The mean position **in the input frame**, i.e. before any pose.
    IMP::algebra::Vector3D position;

    FPSLabelPosition() : body_id(0), position(0, 0, 0) {}

    IMP_SHOWABLE_INLINE(FPSLabelPosition,
                        out << "FPSLabelPosition(" << name << ")");
};
IMP_VALUES(FPSLabelPosition, FPSLabelPositions);

//! The labelling positions of an assembly, at the pose it currently holds.
/*!
    Call it on a **freshly built** assembly: the coordinates are recorded as
    they stand, and the export re-applies each result's pose to them, so a
    position read at an already-docked pose would be transformed twice.

    Volumes are resampled first, so the coordinate is the volume's mean
    position and not wherever the particle happened to be left. `XYZ` and
    `ATOM` positions are taken as they are -- they have no volume to sample.
*/
IMPBFFEXPORT FPSLabelPositions fps_label_positions(
        const DockingAssembly& assembly);

//! One row of FPS's results grid -- one repetition of a run.
struct IMPBFFEXPORT FPSResultRow {
    //! FPS's `InternalNumber`: 1-based, stable, and what every file name uses.
    int number;
    //! FPS's `E`: the restraint sum, **raw, not reduced, clash excluded**.
    /*! `SimulationResult.cs:19` says so, and `SaveEnergyTable` writes it
        undivided while the GUI divides by `dof` -- so the file and the screen
        disagree. Both are written here, as `chi2` and `chi2_r`; the divisor is
        #FPSResultTable::get_dof. */
    double chi2;
    //! FPS's `Ebond`: the part of #chi2 the atom-to-atom restraints contribute.
    /*! **A subset of #chi2, never an addition to it.** */
    double chi2_bond;
    //! FPS's `Eclash`: the excluded-volume term, which is **not** in #chi2.
    double chi2_clash;
    bool converged;
    //! FPS's `ParentStructure`: which row this one branched from; 0 = none.
    /*! The tree shape of a run lives only here -- which replica came from
        which docked solution -- and it is the one thing FPS's exports drop
        entirely (§9 of the reference). It is a column in this port's table. */
    int parent;
    //! FPS's `SimulationMethod`, e.g. `Docking`, `ErrorEstimation`.
    std::string method;
    //! RMSD against the previous exported row, and against the reference, Å.
    /*! Filled by #add_rmsd_columns; `rmsd_previous` of the first row is
        written as `---`, which is what FPS writes. */
    double rmsd_previous, rmsd_reference;
    //! The pose, as #IMP::bff::capture_poses writes it.
    std::string poses;
    //! The experiment-against-model table this row was scored on.
    PairDistances pairs;

    //! The superposition onto the table's reference, as PyMOL applies it.
    /*! **Translate first, then rotate about the origin** -- FPS's ordering,
        which `SimulationResult.RMSD` also assumes. The default is the
        identity (angle 0), not FPS's default-constructed all-zero matrix,
        whose axis-angle is a 180 degree rotation about z. */
    IMP::algebra::Vector3D best_fit_translation;
    IMP::algebra::Vector3D best_fit_axis;
    //! Degrees.
    double best_fit_angle;

    FPSResultRow()
        : number(0), chi2(0), chi2_bond(0), chi2_clash(0), converged(false),
          parent(0), method("Unknown"), rmsd_previous(-1), rmsd_reference(0),
          best_fit_translation(0, 0, 0), best_fit_axis(0, 0, 1),
          best_fit_angle(0) {}

    IMP_SHOWABLE_INLINE(FPSResultRow,
                        out << "FPSResultRow(" << number << ", chi2 " << chi2
                            << ")");
};
IMP_VALUES(FPSResultRow, FPSResultRows);

//! Every repetition of a run, and what the whole table needs to be written.
struct IMPBFFEXPORT FPSResultTable {
    FPSResultRows rows;
    //! How many rigid bodies, for the reduced \f$\chi^2\f$'s degrees of freedom.
    int n_molecules;
    //! The `number` of the row RMSDs and best fits are taken against.
    /*! Non-positive means the **first exported row**, which is what FPS falls
        back to when the user has picked no reference
        (`SaveForm.cs:359`) -- after its in-place sort, the lowest
        `InternalNumber` in the exported set and *not* the best-scoring one. */
    int reference_number;
    //! The `DistanceDataType` the R table's file name carries.
    std::string distance_type;

    FPSResultTable()
        : n_molecules(1), reference_number(0), distance_type("Rmp") {}

    //! FPS's own \f$dof = \max(N_{dist} - 6(N_{mol}-1), 1)\f$.
    /*! `MainForm.cs:683`. Six rigid-body degrees of freedom per body beyond
        the first, and a floor of 1 so a small project does not divide by
        zero or by a negative. */
    int get_dof() const;
    //! The whole table as JSON -- what replaces `SimulationResults.bin`.
    /*! Everything the binary carried except the molecules, which it did not
        carry either, plus the molecule manifest the binary's "re-attach by
        position" hazard needed and never had. */
    std::string get_json() const;

    IMP_SHOWABLE_INLINE(FPSResultTable,
                        out << "FPSResultTable(" << rows.size() << " rows)");
};
IMP_VALUES(FPSResultTable, FPSResultTables);

//! One docking result as one row of the table.
/*! \param[in] result what a run returned
    \param[in] number the 1-based `InternalNumber`
    \param[in] method the provenance, e.g. `Docking`
    \param[in] parent the row this one branched from; 0 = none

    #FPSResultRow::chi2 is `score - e_clash`, i.e. FPS's `E`: this module's
    score is the objective the minimiser descended and has the clash term in
    it, FPS's has not. */
IMPBFFEXPORT FPSResultRow fps_result_row(const DockingResult& result,
                                         int number,
                                         const std::string& method = "Docking",
                                         int parent = 0);

//! An error-estimation run as a results table: the parent, then its replicas.
/*! Row 1 is the parent (`method` `Docking`, `parent` 0); rows 2..N+1 are the
    replicas (`method` `ErrorEstimation`, `parent` 1). The RMSD-vs-reference
    column is already the replicas' RMSD against the parent -- the statistic
    the run exists to report -- because #FPSResultTable::reference_number is
    set to the parent's. */
IMPBFFEXPORT FPSResultTable fps_bootstrap_table(const BootstrapResult& result,
                                                int n_molecules = 1);

//! Fill both RMSD columns of a table from its poses.
/*!
    \param[in] table the rows, whose `poses` must be set
    \param[in] assembly an assembly built from the same inputs
    \param[in] fps_sign_convention accumulate FPS's \f$|Ur - t|^2\f$ instead of
               \f$|Ur + t|^2\f$ (#IMP::bff::pose_rmsd)

    "Previous" is the previous **row of this table**, not the previous in
    simulation time -- FPS's meaning, and worth knowing before reading the
    column.
*/
/*! \return the table with both columns filled

    Returned rather than filled in place because a value type crossing SWIG
    cannot be a non-const reference parameter -- IMP's `IMP_SWIG_VALUE` refuses
    it by name, and the refusal is right: a Python caller would be mutating a
    copy. */
IMPBFFEXPORT FPSResultTable add_rmsd_columns(const FPSResultTable& table,
                                             const DockingAssembly& assembly,
                                             bool fps_sign_convention = false);

//! Compute every row's superposition onto the table's reference.
/*!
    This is the fix for FPS's biggest export trap: there, the best fit baked
    into an exported script is whatever the GUI last computed -- pairwise
    fits, or an all-zero matrix -- because no `Save*` method computes it. Here
    it is computed against #FPSResultTable::reference_number at the moment it
    is needed, so an Overlay is an overlay.

    \param[in] table the rows, whose `poses` must be set
    \param[in] assembly an assembly built from the same inputs
    \return the table with every row's best fit filled in
*/
IMPBFFEXPORT FPSResultTable add_best_fit(const FPSResultTable& table,
                                        const DockingAssembly& assembly);

//! What the writers are told beyond the table itself.
struct IMPBFFEXPORT FPSExportOptions {
    //! The file-name prefix. FPS defaults to `structure` in Dock mode and
    //! `screening_` in Filter mode, and then ignores it for four of the six
    //! outputs; here it is used for **all** of them.
    std::string prefix;
    //! Emit the best-fit superposition into the PyMOL scripts.
    bool best_fit;
    //! Add `select all` / `save <name>.pdb, sele` to each per-result script.
    bool save_pdb;
    //! Write `camera=0` on every `rotate`/`translate`. See the file comment.
    bool camera_zero;
    //! Number the Filter R table's rows FPS's way: the 0-based array index.
    /*! **FPS's off-by-one.** Everything else in FPS -- the Dock R table, the
        chi2 table, every file name -- prints the 1-based `InternalNumber`;
        `SaveFilterDistances` prints `j`, the position in the exported array
        (`SaveForm.cs:463`). It is not even stable, because exporting a
        selection renumbers it. The default here is 1-based; **join on `File`
        either way.** */
    bool fps_filter_number;
    //! Where OverlayStates writes its scratch PDB; empty is `<dir>/_tmp.pdb`.
    /*! FPS writes the bare relative name `_tmp.pdb`, so it lands in PyMOL's
        current working directory, fails silently if that is not writable, and
        is left behind. */
    std::string tmp_pdb;
    //! The pseudoatoms to add; empty adds none.
    FPSLabelPositions labels;

    FPSExportOptions()
        : prefix("structure"), best_fit(false), save_pdb(false),
          camera_zero(true), fps_filter_number(false) {}

    IMP_SHOWABLE_INLINE(FPSExportOptions,
                        out << "FPSExportOptions(prefix=" << prefix << ")");
};
IMP_VALUES(FPSExportOptions, FPSExportOptionsList);

//! `<dir>/<prefix><number>.pml`: one docking solution, reconstructed in PyMOL.
/*! Loads the untransformed inputs and applies the pose as `rotate`/`translate`
    commands. **All \f$N\f$ transforms are written**, including molecule 0's:
    FPS skips it, which is valid only while its own engine normalises the pose
    so that molecule 0 is the identity, and is wrong for any pose that does not
    come from that engine.
    \return the path written
    \throw IOException when it cannot be written */
IMPBFFEXPORT std::string write_fps_pymol_script(
        const std::string& directory, const FPSResultRow& row,
        const FPSMolecules& molecules,
        const FPSExportOptions& options = FPSExportOptions());

//! One #write_fps_pymol_script per row. \return the paths written
IMPBFFEXPORT std::vector<std::string> write_fps_pymol_scripts(
        const std::string& directory, const FPSResultTable& table,
        const FPSMolecules& molecules,
        const FPSExportOptions& options = FPSExportOptions());

//! `<dir>/<prefix>Overlay.pml`: every result as its own single-state object.
/*! \return the path written */
IMPBFFEXPORT std::string write_fps_overlay(
        const std::string& directory, const FPSResultTable& table,
        const FPSMolecules& molecules,
        const FPSExportOptions& options = FPSExportOptions());

//! `<dir>/<prefix>OverlayStates.pml`: one object, one state per result.
/*! Each state also gets a `set_title` naming its `InternalNumber`; FPS loses
    that mapping entirely, leaving only the header comment and the assumption
    that the reader knows the order.
    \return the path written */
IMPBFFEXPORT std::string write_fps_overlay_states(
        const std::string& directory, const FPSResultTable& table,
        const FPSMolecules& molecules,
        const FPSExportOptions& options = FPSExportOptions());

//! `<dir>/<prefix>Rtable_<DataType>.txt`: model distances, one row per result.
/*!
    Tab separated. Header `Structure`, `Number`, then one column per distance
    named `<position1>_<position2>` in the table's own order.

    **One file, not FPS's two.** FPS writes `Rtable_Rmp.txt` plus, when the
    project's data type is not `Rmp`, a second file whose values are the *same*
    \f$R_{mp}\f$ pushed through a global polynomial -- an estimate, not a
    simulation, while the identically-named Filter-mode file *is* simulated.
    Here the model distance is simulated end to end, so there is one file and
    its name says which observable is in it.

    \return the path written */
IMPBFFEXPORT std::string write_fps_r_table(
        const std::string& directory, const FPSResultTable& table,
        const FPSExportOptions& options = FPSExportOptions());

//! `<dir>/<prefix>chi2table.txt`: the results grid, one row per repetition.
/*!
    Tab separated, ten columns:

        File  chi2  chi2_r  chi2_bond  chi2_clash  Converged  Method  Parent
              RMSD vs previous  RMSD vs <M>

    against FPS's nine. The differences are all deliberate: `chi2_r` is added
    because FPS's file writes the raw sum under a header that reads as a
    reduced one and its GUI shows the reduced one; `Parent` is added because
    the branch structure of a run is otherwise unrecoverable; and FPS's
    `RMSD vs <M> (selected)` is dropped because it is `sqrt(0/0)` = NaN
    whenever no molecule is flagged, which is the common case (its D17).

    \return the path written */
IMPBFFEXPORT std::string write_fps_chi2_table(
        const std::string& directory, const FPSResultTable& table,
        const FPSExportOptions& options = FPSExportOptions());

//! Every Dock-mode export in one call. \return the paths written, in order
IMPBFFEXPORT std::vector<std::string> write_fps_exports(
        const std::string& directory, const FPSResultTable& table,
        const FPSMolecules& molecules,
        const FPSExportOptions& options = FPSExportOptions());

//! `<dir>/<prefix>Rtable_<DataType>.txt` for a screen. \return the path
/*! Columns `File`, `Number`, then one per distance. `File` is the structure's
    basename and is the **only** safe join key -- see
    #FPSExportOptions::fps_filter_number. */
IMPBFFEXPORT std::string write_fps_screening_r_table(
        const std::string& directory,
        const std::vector<ScreenedStructure>& structures,
        const FPSExportOptions& options = FPSExportOptions());

//! `<dir>/<prefix>chi2table.txt` for a screen. \return the path
/*!
    FPS's seven Filter-mode columns exactly:
    `File`, `Chi2r`, `NaNs`, `RefRMSD`, `>1sigma`, `>2sigma`, `>3sigma`.

    Two things to read together. `Chi2r` **is** reduced here -- by the number
    of *scored* distances, not by a degrees-of-freedom count -- so it is not
    comparable with the Dock-mode `chi2`. And distances with no model value are
    excluded from it and counted in `NaNs` instead, so a structure with many
    empty volumes gets a flatteringly low `Chi2r` from very few restraints.
    The sigma counts nest: a 3.5-sigma outlier is in all three.

    \return the path written */
IMPBFFEXPORT std::string write_fps_screening_chi2_table(
        const std::string& directory,
        const std::vector<ScreenedStructure>& structures,
        const FPSExportOptions& options = FPSExportOptions());

IMPBFF_END_NAMESPACE

#endif //IMPBFF_FPSEXPORT_H
