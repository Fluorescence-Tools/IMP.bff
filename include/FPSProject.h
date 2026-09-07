/**
 *  \file IMP/bff/FPSProject.h
 *  \brief The document that binds a FRET docking or screening run together.
 *
 * A run needs more than a command line: the structures in body order, where
 * the labelling and the distances come from, which distances are selected, the
 * per-mode search parameters, the AV globals, and the pose the bodies are
 * currently in. FPS calls that a *project* and stores it in one file; this
 * module had no such object, so every run was a twelve-flag invocation and
 * nothing about it was repeatable (PRD-121 G1).
 *
 * **Two doors, and only one of them writes** (PRD-121 D3):
 *
 * 1. #IMP::bff::read_fps_project_txt reads the **text twin** FPS writes beside
 *    every project it saves (`MainForm.cs:559-560` calls `om.Save(f)` and then
 *    `om.Export(f + ".txt")`, so the twin always exists). The `.bin` itself is
 *    an explicit non-goal: it is a gzipped .NET `BinaryFormatter` dump of a
 *    `Dictionary<string, Dictionary<string, Option>>` carrying
 *    assembly-qualified type names, and a C++ reader for MS-NRBF is not worth
 *    the readable file sitting next to it.
 * 2. #IMP::bff::read_fps_project and #IMP::bff::write_fps_project read and
 *    write an **extended fps.json** carrying the project as its own `Project`
 *    section, beside `Positions`, `Distances` and `χ²`. That file can be
 *    self-contained -- one path for a whole run.
 *
 * Writing the legacy format is a non-goal, exactly as `FPSIO.h` argues for the
 * other two: the formats stay readable so a decade of measurements is not
 * lost, and nothing should produce another one.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_FPSPROJECT_H
#define IMPBFF_FPSPROJECT_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/Docking.h>

#include <IMP/bff/Base.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! FPS's five run modes, spelled as its project file spells them.
/*! `ProjectData.cs:10-11`: the first four are `DockModes`, the fifth is the
    lone `FilterModes` entry. The strings are the keys of the `parameters`
    object in the fps.json form and the names #fps_mode_parameters takes, so
    they are the one spelling and there is no second table to drift. */
IMPBFFEXPORT std::vector<std::string> fps_mode_names();

//! One mode's search parameters -- FPS's `FPSParameters` (`SimulationBase.cs:11`).
/*!
    The ten FPS fields are carried **verbatim**, including the ones this port
    does not use, so a legacy project survives a read and a caller can see what
    the run that produced a published number was told. What each one does:

    * `viscosity_factor`, `time_step_factor`, `e_tolerance`, `k_tolerance`,
      `f_tolerance`, `t_tolerance` drive FPS's damped rigid-body integrator.
      **This port does not reproduce that integrator** (PRD-121 D1), so they
      are stored and not read. They are not dropped because dropping them makes
      the file lossy and the provenance unrecoverable.
    * `max_iterations` is FPS's integrator step budget, and is **not** IMP's
      iteration count -- 200 000 damped Verlet steps and 500 conjugate-gradient
      iterations are not the same quantity and converting between them would be
      a fabricated number. #iterations is this port's own budget and sits
      beside it.
    * `max_force` and `clash_tolerance` **do** change the answer here and are
      honoured: #IMP::bff::DockingParameters::max_force and
      #IMP::bff::DockingParameters::clash_tolerance (PRD-121 G7).
    * `rkt` is \f$1/kT\f$ and is read by Sample mode alone.
    * `optimize_selected` is `"Selected"`, `"All"` or `"SelectedThenAll"`;
      #IMP::bff::DockingParameters normalises the spelling, so FPS's own
      casing passes straight through.
*/
struct IMPBFFEXPORT FPSModeParameters {
    double viscosity_factor, time_step_factor;
    int max_iterations;
    double max_force, clash_tolerance;
    //! FPS's `rkT`, i.e. \f$1/kT\f$. Sample mode is the only reader.
    double rkt;
    //! `ETolerance` is dead in every FPS mode; it is stored, not obeyed.
    double e_tolerance, k_tolerance, f_tolerance, t_tolerance;
    std::string optimize_selected;
    //! This port's iteration budget -- **not** FPS's #max_iterations.
    int iterations;

    FPSModeParameters()
        : viscosity_factor(1.0), time_step_factor(1.0), max_iterations(200000),
          max_force(400.0), clash_tolerance(1.0), rkt(10.0),
          e_tolerance(100.0), k_tolerance(0.001), f_tolerance(0.001),
          t_tolerance(0.02), optimize_selected("Selected"), iterations(500) {}

    //! The block as a JSON object, as the `Project.parameters` section holds it.
    std::string get_json() const;

    IMP_SHOWABLE_INLINE(FPSModeParameters,
                        out << "FPSModeParameters(max_force=" << max_force
                            << ", clash_tolerance=" << clash_tolerance << ", "
                            << optimize_selected << ")");
};
IMP_VALUES(FPSModeParameters, FPSModeParametersList);

//! The parameters FPS ships for one mode (`ProjectData.cs:61-149`).
/*! The table is in `okf/references/fps-sampling-and-protocol.md` §6. An
    unknown mode name gives the `Dock` block, which is FPS's own base case.
    \param[in] mode one of #fps_mode_names */
IMPBFFEXPORT FPSModeParameters fps_mode_parameters(const std::string& mode);

//! FPS's `ConversionParameters` (`ConversionForm.cs:12`).
struct IMPBFFEXPORT FPSConversionParameters {
    //! `R0`, the Förster radius the conversion polynomial is fitted at, Å.
    double forster_radius;
    //! `PolynomOrder`.
    int polynomial_order;

    FPSConversionParameters() : forster_radius(52.0), polynomial_order(3) {}

    IMP_SHOWABLE_INLINE(FPSConversionParameters,
                        out << "FPSConversionParameters(R0="
                            << forster_radius << ")");
};
IMP_VALUES(FPSConversionParameters, FPSConversionParametersList);

//! FPS's `AVGlobalParameters` (`AVEngine.cs:32`).
/*!
    `grid_size` is **relative**: the grid spacing follows from the AV
    parameters as \f$dg = \max(\min(0.2L, 0.2W, 0.4R_1, 0.4R_2, 0.4R_3),
    0.4)\f$ (`AVEngine.cs:568`), where 0.2 is #grid_size and 0.4 is
    #min_grid_size. Both shipped reference clouds come out at dg = 0.6 Å and
    this module builds the same grid (PRD-121, phase 0).

    #linker_initial_sphere is **not** this module's clearance and must not be
    copied into one: FPS's 0.5 and `IMP.bff`'s derived
    `allowed_sphere_radius` are different quantities, which is the trap that
    made a linker width of 4.5 Å return a silently empty volume. It is stored
    for provenance.
*/
struct IMPBFFEXPORT FPSAVGlobalParameters {
    double grid_size, min_grid_size, linker_initial_sphere;
    int link_search_nodes, e_samples;

    FPSAVGlobalParameters()
        : grid_size(0.2), min_grid_size(0.4), linker_initial_sphere(0.5),
          link_search_nodes(3), e_samples(200000) {}

    IMP_SHOWABLE_INLINE(FPSAVGlobalParameters,
                        out << "FPSAVGlobalParameters(grid_size=" << grid_size
                            << ")");
};
IMP_VALUES(FPSAVGlobalParameters, FPSAVGlobalParametersList);

//! Everything a run is told, in one document.
/*!
    Default-constructed it is FPS's shipped project: the five parameter blocks
    of `ProjectData()` and its AV and conversion globals, with no structures
    and no labelling source yet.
*/
struct IMPBFFEXPORT FPSProject {
    //! Where this project was read from; empty for one built in memory.
    /*! It is what a relative #structures or #labelling_json entry is resolved
        against, and it is what #get_labelling_path falls back to for a
        self-contained file. */
    std::string path;
    //! FPS's `ProjectFPSMode`: `"None"`, `"Dock"` or `"Filter"`.
    /*! FPS's own two-way split -- everything that moves bodies is `Dock`,
        screening is `Filter`. Kept as FPS spells it. */
    std::string mode;
    //! Provenance: the first `#` line of a legacy text export, or empty.
    std::string source;

    //! One structure per rigid body, **in body order**.
    /*! FPS's `MoleculesPaths`. The order is the body order everywhere:
        `--fixed-body 0` is this list's first entry, and a captured pose is
        indexed by it. */
    std::vector<std::string> structures;

    //! The labelling source, door one: an fps.json path.
    /*! Empty and #positions_path empty means the project file **is** the
        labelling file -- the self-contained form. */
    std::string labelling_json;
    //! The labelling source, door two: FPS's legacy positions `.txt`.
    std::string positions_path;
    //! The distances `.txt` beside it. FPS's reader wants it named
    //! `Distances.txt` and living in the positions file's directory.
    std::string distances_path;

    //! The fps.json score set this run scores, or empty for every distance.
    std::string score_set;

    //! Which distances are selected, **by name**.
    /*! Empty means "whatever #score_set says", which is the ordinary case for
        an fps.json project. This is the in-tree spelling because an fps.json
        distance has a name; FPS's is positional (see #selected_flags) and a
        position in a list is not a durable identifier. */
    std::vector<std::string> selected_distances;
    //! FPS's `SelectedDistances`, a `Boolean[]` parallel to the distances file.
    /*! Filled by the legacy text reader and by nothing else. 1 selected, 0
        not. It survives a round trip so a legacy project is not silently
        re-interpreted, and #resolve_selected_distances turns it into names
        once the distance list is known. */
    std::vector<int> selected_flags;

    //! The `body_id` held fixed while the others move.
    int fixed_body;
    //! Weight of the excluded-volume term.
    double ev_weight;
    //! Width of the mean-position transfer function, Å.
    double sigma_da;
    //! Initial random displacement of the mobile bodies, Å; 0 does not shuffle.
    double shuffle;
    //! FPS-style refinement cycles after a dock.
    int refine_av_cycles;

    //! Clash detection on one coarse bead per residue rather than all atoms.
    /*! #IMP::bff::DockingParameters::coarse_clash. Stored because it decides
        whether #clash_radii_source and #clash_radii_scale can be honoured at
        all -- a bead is not an atom -- so a project carrying one without the
        other would resume into an exception. */
    bool coarse_clash;
    //! Which van der Waals radii the excluded-volume term measures overlap by.
    /*! `"imp"` (the default, the particles' own united-atom radii) or
        `"olga"` (the Bondi-scale table FPS's `ClashTolerance` was calibrated
        against). #IMP::bff::DockingParameters::clash_radii_source. */
    std::string clash_radii_source;
    //! Multiplies whatever #clash_radii_source gives; 1.0 changes nothing.
    double clash_radii_scale;

    //! The current rigid-body poses, as #IMP::bff::capture_poses writes them.
    /*! A JSON array of `{"body_id", "t", "q"}`. `"[]"` means "the input
        structures as they stand"; anything else is a state a run reached, so
        storing the project stores the run and #IMP::bff::apply_poses continues
        it. */
    std::string poses;

    FPSModeParameters dock, refine, error_estimation, sample, screening;
    FPSConversionParameters conversion;
    FPSAVGlobalParameters av;

    FPSProject();

    //! The parameter block of one mode. \param[in] mode one of #fps_mode_names
    /*! An unknown name gives #dock, matching #fps_mode_parameters. */
    FPSModeParameters get_parameters(const std::string& mode) const;
    //! Replace one mode's block. An unknown name is ignored.
    void set_parameters(const std::string& mode,
                        const FPSModeParameters& parameters);

    //! The file the readers should be handed for labelling and distances.
    /*! #labelling_json if set, else #positions_path (the legacy pair, which
        `read_fps_json` dispatches on by extension), else #path -- the
        self-contained case, where the project file carries its own
        `Positions` and `Distances`.

        **Paths are stored as written and resolved only here.** A project
        saying `hiv_rt.fps.json` still says so after a read and a write, and
        the directory can be moved or shipped. Resolving at read time would
        write absolute paths back out and pin the file to one machine. */
    std::string get_labelling_path() const;

    //! #distances_path, resolved against #path.
    std::string get_distances_path() const;

    //! #structures, resolved against #path where they were written relative.
    /*! A Windows absolute path (`C:\...`) counts as absolute and is left
        alone: gluing a POSIX directory in front of one produces a path that
        is wrong in a way nobody can read. */
    std::vector<std::string> get_structure_paths() const;

    //! The mode's block, translated into what a run here is actually told.
    /*!
        Only the parameters this port honours cross over: #max_force,
        #clash_tolerance and #optimize_selected from the block, and the
        project's #fixed_body, #ev_weight, #sigma_da, #shuffle,
        #refine_av_cycles and #score_set. FPS's integrator settings do not,
        because there is no integrator here to give them to (PRD-121 D1), and
        #max_iterations does not, because #iterations is the budget of a
        different optimiser.
    */
    DockingParameters get_docking_parameters(
            const std::string& mode = "Dock") const;

    //! Turn #selected_flags into names, given the distance list they index.
    /*!
        \param[in] distance_names the distances in file order -- the order
               FPS's `Boolean[]` was written in
        \return the selected names; empty when there are no flags. A flag
                array shorter than the list selects only what it covers, and a
                longer one is truncated: a legacy file may have been written
                against a distance file that has since changed, and guessing
                past the end would invent a selection.
    */
    std::vector<std::string> resolve_selected_distances(
            const std::vector<std::string>& distance_names) const;

    //! What is wrong with this project, in the order a user would fix it.
    /*! Empty means it can be run. A missing file is reported by name; a
        parameter out of range is reported with the value. This is a *report*,
        not an exception: a project being edited is allowed to be incomplete
        and a caller wants the whole list, not the first problem. */
    std::vector<std::string> get_problems() const;

    //! The `Project` section, as JSON text.
    std::string get_json() const;

    IMP_SHOWABLE_INLINE(FPSProject,
                        out << "FPSProject(" << structures.size()
                            << " structures, mode " << mode << ")");
};
IMP_VALUES(FPSProject, FPSProjects);

//! Build a project from the `Project` section of an fps.json, as JSON text.
/*! \param[in] project_json the section; `{}` gives the defaults
    \param[in] path what relative entries are resolved against, and what
           #FPSProject::path is set to
    \throw ValueException when the text is not JSON */
IMPBFFEXPORT FPSProject fps_project_from_json(const std::string& project_json,
                                              const std::string& path = "");

//! Read the text twin of a legacy FPS `.bin` project (`OptionsManager::Export`).
/*!
    The format, in full, because it has no specification anywhere else
    (`GPFileTools/OptionsManager.cs:448-489`):

    ~~~
    # FPS v. 1.1.0.0                      <- banner, then the save date
    # Monday, May 19, 2014 11:22:33 AM

    # Data                                <- a category
    ProjectFPSMode = Dock                  <- Key = value.ToString()
    MoleculesPaths = System.String[]       <- an array prints its *type*...
     protein.pdb dna.pdb                   <- ...and its elements on the NEXT line
    LabelingPositionsPath = LPs.txt

    # Search parameters
    DockParameters = Fps.FPSParameters     <- a struct prints its type...
     ViscosityFactor = 1                   <- ...then one indented line per field
     MaxIterations = 200000
    ~~~

    Four things the writer does that a reader has to know, all of them found by
    reading `Export` rather than by guessing:

    1. **An array's elements are space-separated on the following line**, and
       nothing is quoted or escaped (`sr.Write(" " + ov.ToString())`). A
       structure path containing a space is therefore **not recoverable** from
       this format -- FPS loses it, and so does this. The `.bin` twin has it;
       reading that is the non-goal above. A path that splits is reported by
       #FPSProject::get_problems as a missing file, which is the honest
       symptom.
    2. **A zero-length array still writes its line**, an empty one. The
       element line is consumed unconditionally after an array key for exactly
       that reason.
    3. **A struct and an array are told apart by the value on the key line**
       -- `System.X[]` versus a type name -- not by looking at what follows.
       Both are indented by one space and ` Field = value` is a legal array
       element.
    4. **Numbers are written in the writer's culture.** `Double.ToString()`
       with a German locale writes `0,0005`. A value that will not parse with
       a `.` is retried with `,` as the decimal point, because the alternative
       is silently reading a tolerance as zero.

    Unknown categories and unknown keys are skipped rather than refused: this
    is an options dump of a GUI, other builds of it carried other keys, and a
    project that mostly parses is worth more than an exception.

    \param[in] path the `.txt` twin
    \throw IOException when the file cannot be read
*/
IMPBFFEXPORT FPSProject read_fps_project_txt(const std::string& path);

//! Read a project from either door, dispatching on the extension.
/*! `.json` (so `.fps.json` too) is read as an extended fps.json and its
    `Project` section becomes the project; anything else is read as a legacy
    text export. A `.json` with no `Project` section is not an error -- it is a
    plain labelling file, and the project that comes back is the defaults
    pointed at it, which is what "run this fps.json" should mean.
    \throw IOException when the file cannot be read
    \throw ValueException when a `.json` file is not JSON */
IMPBFFEXPORT FPSProject read_fps_project(const std::string& path);

//! Write an extended fps.json carrying the project in its own section.
/*!
    The three labelling sections are written exactly as
    #IMP::bff::write_fps_json writes them, so a project file **is** an
    fps.json: every existing reader takes it, and a self-contained project is
    one path for a whole run.

    \param[in] path where to write
    \param[in] project the project; its `Project` section is written under the
           key `Project`
    \param[in] positions_json,distances_json,score_sets_json the labelling
           sections. All three `{}` writes a project that points at its
           labelling source instead of carrying it.
    \param[in] extra_json any other top-level keys to keep
    \param[in] validate check the assembled payload and refuse to write a
           non-conforming file
    \throw ValueException when \p validate is set and the payload does not
           conform
    \throw IOException when the file cannot be written
*/
IMPBFFEXPORT void write_fps_project(const std::string& path,
                                    const FPSProject& project,
                                    const std::string& positions_json = "{}",
                                    const std::string& distances_json = "{}",
                                    const std::string& score_sets_json = "{}",
                                    const std::string& extra_json = "{}",
                                    bool validate = true);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_FPSPROJECT_H
