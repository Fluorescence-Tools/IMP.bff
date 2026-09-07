/**
 * \file FPSExport.cpp
 * \brief What a docking or screening run is written out as: FPS's six exports.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/FPSExport.h>

#include <IMP/bff/AV.h>
#include <IMP/bff/FPSIO.h>
#include <IMP/bff/internal/Text.h>
#include <IMP/bff/internal/json.h>

#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/pdb.h>
#include <IMP/core/XYZ.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/algebra/Rotation3D.h>
#include <IMP/bff/Base.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! Fixed-point text in the C locale, whatever the process locale is.
/*! FPS forces `NumberDecimalSeparator = "."` at startup (`Program.cs:19-22`)
    and the `F` specifier never groups thousands, so a C-locale fixed-point
    write reproduces every number it prints. `ostringstream` would otherwise
    follow the global locale, which a *caller* -- a GUI, a notebook -- may have
    changed, and a comma decimal separator inside a comma-separated PyMOL
    vector is not a formatting nit, it is a wrong script. */
std::string fixed_text(double value, int digits) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    if (!std::isfinite(value)) {
        // `Double.NaN.ToString("F3")` is the literal `NaN`, and a reader that
        // sees it learns something a blank cell would have hidden.
        return value != value ? std::string("NaN")
                              : std::string(value > 0 ? "Inf" : "-Inf");
    }
    out << std::fixed << std::setprecision(digits) << value;
    std::string text = out.str();
    // `-0.000` is a sign on a zero. It reaches a file whenever a component
    // rounds to zero from below -- an identity transform written as
    // `translate [-0.000, -0.000, -0.000]` -- and it makes two identical
    // exports compare unequal byte for byte. Dropped, on the display side
    // only: no arithmetic sees it.
    if (text.size() > 1 && text[0] == '-' &&
        text.find_first_not_of("-0.") == std::string::npos) {
        text.erase(0, 1);
    }
    return text;
}

//! `Vector3.ToString()`: three decimals, comma and space (`MatrixVector3.cs:62`).
std::string vector_text(const IMP::algebra::Vector3D& v) {
    return fixed_text(v[0], 3) + ", " + fixed_text(v[1], 3) + ", " +
           fixed_text(v[2], 3);
}

std::string integer_text(int v) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << v;
    return out.str();
}

//! The basename of a path, extension removed -- FPS's `Molecule.Name`.
std::string export_stem_of(const std::string& path) {
    const std::size_t slash = path.find_last_of("/\\");
    const std::string base =
            slash == std::string::npos ? path : path.substr(slash + 1);
    const std::size_t dot = base.find_last_of('.');
    return dot == std::string::npos || dot == 0 ? base : base.substr(0, dot);
}

//! The basename of a path, extension kept -- FPS's `ShortFileName`.
std::string export_basename_of(const std::string& path) {
    const std::size_t slash = path.find_last_of("/\\");
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

//! `<prefix><number>`, the name every FPS file column and file name uses.
std::string result_name(const FPSExportOptions& options, int number) {
    return options.prefix + integer_text(number);
}

//! `, camera=0` or nothing. See the header: PyMOL's default is camera space.
std::string camera(const FPSExportOptions& options) {
    return options.camera_zero ? std::string(", camera=0") : std::string();
}

//! One molecule's pose as PyMOL applies it: rotate about CM, then translate.
struct PymolTransform {
    IMP::algebra::Vector3D axis;
    double angle_degrees;
    IMP::algebra::Vector3D origin;
    IMP::algebra::Vector3D translation;
    //! The same rotation as a rotation, for transporting a labelling position.
    IMP::algebra::Rotation3D rotation;
    PymolTransform()
        : axis(0, 0, 1), angle_degrees(0), origin(0, 0, 0), translation(0, 0, 0),
          rotation(IMP::algebra::get_identity_rotation_3d()) {}
};

//! The transforms of one result, one per molecule, in molecule order.
/*!
    `capture_poses` writes one entry per rigid body in body order and
    `build_docking_assembly` creates one body per input PDB in the order given,
    so entry *i* is molecule *i*. The `body_id` is carried in the JSON and is
    checked against nothing here on purpose: matching by position is what the
    assembly itself does, and a second matching rule would be a second place
    for the two to disagree.

    **The pose is not the transform.** `IMP::atom::create_rigid_body` starts a
    body in its principal-axis frame, not at the identity, so the recorded pose
    \f$T_{pose}\f$ has to be composed with the *input* frame's inverse before
    it means anything to a script that loads the untransformed PDB:
    \f$M = T_{pose} \circ T_{input}^{-1}\f$. Getting this wrong is silent --
    every structure is rotated by the body's principal-axis rotation and the
    overlay still looks like an overlay, because every structure is rotated by
    the *same* amount.

    PyMOL is then given the pair FPS writes, `rotate` about `CM` followed by
    `translate`: from \f$x \mapsto R x + v\f$ and
    \f$x \mapsto R(x - CM) + CM + T\f$, the translation is
    \f$T = v + R\,CM - CM\f$.
*/
std::vector<PymolTransform> transforms_of(const std::string& poses_json,
                                          const FPSMolecules& molecules) {
    std::vector<PymolTransform> out(molecules.size());
    for (std::size_t i = 0; i < molecules.size(); ++i) {
        out[i].origin = molecules[i].center;
    }
    const nlohmann::json poses =
            nlohmann::json::parse(poses_json.empty() ? "[]" : poses_json, NULL,
                                  false);
    if (poses.is_discarded() || !poses.is_array()) return out;
    for (std::size_t i = 0; i < poses.size() && i < out.size(); ++i) {
        const nlohmann::json& pose = poses[i];
        if (!pose.is_object() || !pose.contains("t") || !pose.contains("q")) {
            continue;
        }
        const std::vector<double> t = pose["t"].get<std::vector<double> >();
        const std::vector<double> q = pose["q"].get<std::vector<double> >();
        if (t.size() != 3 || q.size() != 4) continue;
        const IMP::algebra::Transformation3D pose_frame(
                IMP::algebra::Rotation3D(q[0], q[1], q[2], q[3]),
                IMP::algebra::Vector3D(t[0], t[1], t[2]));
        const IMP::algebra::Transformation3D move =
                pose_frame * molecules[i].frame.get_inverse();
        const IMP::algebra::Rotation3D rotation = move.get_rotation();
        const std::pair<IMP::algebra::Vector3D, double> axis_angle =
                IMP::algebra::get_axis_and_angle(rotation);
        out[i].axis = axis_angle.first;
        out[i].angle_degrees = axis_angle.second * 180.0 / IMP::algebra::PI;
        out[i].rotation = rotation;
        out[i].translation = move.get_translation() +
                             rotation.get_rotated(molecules[i].center) -
                             molecules[i].center;
    }
    return out;
}

//! The `load` block every script starts with.
void write_loads(std::ostream& out, const FPSMolecules& molecules) {
    for (std::size_t i = 0; i < molecules.size(); ++i) {
        // Unquoted, as FPS writes it -- and as PyMOL's `load` wants it. A path
        // with a space breaks the script either way; the export directory is
        // the caller's choice.
        out << "load " << molecules[i].path << "\n";
    }
}

//! The per-molecule `rotate`/`translate` pair for one result.
void write_pose(std::ostream& out, const std::vector<PymolTransform>& t,
                const FPSMolecules& molecules, const std::string& target,
                const FPSExportOptions& options) {
    for (std::size_t i = 0; i < molecules.size() && i < t.size(); ++i) {
        const std::string name = target.empty() ? molecules[i].name : target;
        out << "rotate [" << vector_text(t[i].axis) << "], "
            << fixed_text(t[i].angle_degrees, 3) << ", " << name
            << ", origin=[" << vector_text(t[i].origin) << "]"
            << camera(options) << "\n";
        out << "translate [" << vector_text(t[i].translation) << "], " << name
            << camera(options) << "\n";
    }
}

//! The best-fit block: translate to the reference frame, then rotate about 0.
/*! FPS's ordering, kept because `SimulationResult.RMSD` assumes it and because
    a reader of an FPS script should recognise the shape. What is *not* FPS's
    is where the numbers come from -- see #IMP::bff::add_best_fit. */
void write_best_fit(std::ostream& out, const FPSResultRow& row,
                    const std::string& selection,
                    const FPSExportOptions& options) {
    out << "translate [" << vector_text(row.best_fit_translation) << "], "
        << selection << camera(options) << "\n";
    out << "rotate [" << vector_text(row.best_fit_axis) << "], "
        << fixed_text(row.best_fit_angle, 3) << ", " << selection
        << ", origin=[0, 0, 0]" << camera(options) << "\n";
}

//! The row a table's RMSD and best fit are taken against.
const FPSResultRow* reference_row(const FPSResultTable& table) {
    if (table.rows.empty()) return NULL;
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        if (table.rows[i].number == table.reference_number &&
            table.reference_number > 0) {
            return &table.rows[i];
        }
    }
    return &table.rows[0];
}

//! The distance columns of a table: first-seen order over every row.
/*! Rows come from one network and so agree, but a table assembled by hand may
    not; taking the union rather than row 0's list means a distance that only
    some rows carry still gets a column, with `NaN` where it is missing. */
std::vector<std::string> distance_columns(const FPSResultTable& table,
                                          std::vector<std::string>* headers) {
    std::vector<std::string> names;
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        const PairDistances& pairs = table.rows[i].pairs;
        for (std::size_t k = 0; k < pairs.size(); ++k) {
            if (std::find(names.begin(), names.end(), pairs[k].name) !=
                names.end()) {
                continue;
            }
            names.push_back(pairs[k].name);
            // FPS's column name, `l1.Name + '_' + l2.Name` -- unescaped, so a
            // position name containing `_` makes the header ambiguous. Kept as
            // FPS spells it because a decade of readers split on it.
            headers->push_back(pairs[k].position1 + "_" + pairs[k].position2);
        }
    }
    return names;
}

double model_distance_of(const PairDistances& pairs, const std::string& name) {
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        if (pairs[i].name == name) return pairs[i].distance_model;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

//! An output stream in the C locale, or an exception.
/*! Named for this file: the module is built as one translation unit, so an
    anonymous namespace does not isolate a helper from a same-named one in
    another `.cpp` -- `StructureIO.cpp` has an `open_for_write` of its own. */
std::ofstream export_stream(const std::string& path) {
    std::ofstream out(path.c_str());
    if (!out) IMP_THROW("cannot write " << path, IOException);
    out.imbue(std::locale::classic());
    return out;
}

}  // namespace

/* ------------------------------------------------------------------------
 * The inputs an export names
 * ------------------------------------------------------------------------ */

FPSMolecules fps_molecules(const std::vector<std::string>& pdb_paths) {
    FPSMolecules out;
    for (std::size_t i = 0; i < pdb_paths.size(); ++i) {
        if (!internal::file_exists(pdb_paths[i])) {
            IMP_THROW("PDB file not found: " << pdb_paths[i], IOException);
        }
        // The same read and the same rigid-body construction
        // `build_docking_assembly` performs, so the centre recorded here is
        // the frame the pose was captured in -- by construction, not by a
        // second formula that could drift from it.
        IMP_NEW(IMP::Model, model, ());
        IMP::atom::Hierarchy h = IMP::atom::read_pdb(
                pdb_paths[i], model,
                new IMP::atom::NonWaterNonHydrogenPDBSelector());
        FPSMolecule molecule;
        molecule.path = pdb_paths[i];
        molecule.name = export_stem_of(pdb_paths[i]);
        molecule.n_atoms = static_cast<int>(IMP::atom::get_leaves(h).size());
        if (molecule.n_atoms > 0) {
            IMP::core::RigidBody rb = IMP::atom::create_rigid_body(h);
            molecule.frame =
                    rb.get_reference_frame().get_transformation_to();
            molecule.center = molecule.frame.get_translation();
        }
        out.push_back(molecule);
    }
    return out;
}

FPSLabelPositions fps_label_positions(const DockingAssembly& assembly) {
    FPSLabelPositions out;
    ProbeNetworkRestraint* network = assembly.get_network();
    if (network == NULL) return out;
    IMP::Model* model = assembly.get_model();

    std::map<std::string, IMP::algebra::Vector3D> coordinates;
    IMP::bff::AVs avs = network->get_used_avs();
    for (std::size_t i = 0; i < avs.size(); ++i) {
        IMP::bff::AV av = avs[i];
        av.resample();
        coordinates.insert(std::make_pair(av.get_particle()->get_name(),
                                          IMP::core::XYZ(av).get_coordinates()));
    }
    const std::map<std::string, IMP::ParticleIndex> points =
            network->get_point_positions();
    for (std::map<std::string, IMP::ParticleIndex>::const_iterator it =
                 points.begin();
         it != points.end(); ++it) {
        coordinates[it->first] =
                IMP::core::XYZ(model, it->second).get_coordinates();
    }

    nlohmann::json positions = nlohmann::json::object();
    if (!assembly.get_fps_json_path().empty()) {
        const FPSDocument document = read_fps_json(assembly.get_fps_json_path());
        positions = nlohmann::json::parse(document.positions, NULL, false);
        if (positions.is_discarded()) positions = nlohmann::json::object();
    }
    // File order first, so the pseudoatoms come out in the order a user wrote
    // them; anything the file does not mention follows, so nothing is dropped.
    std::vector<std::string> order;
    for (nlohmann::json::const_iterator it = positions.begin();
         it != positions.end(); ++it) {
        if (coordinates.count(it.key()) > 0) order.push_back(it.key());
    }
    for (std::map<std::string, IMP::algebra::Vector3D>::const_iterator it =
                 coordinates.begin();
         it != coordinates.end(); ++it) {
        if (std::find(order.begin(), order.end(), it->first) == order.end()) {
            order.push_back(it->first);
        }
    }
    for (std::size_t i = 0; i < order.size(); ++i) {
        FPSLabelPosition label;
        label.name = order[i];
        label.position = coordinates[order[i]];
        if (positions.contains(order[i]) && positions[order[i]].is_object()) {
            label.body_id = positions[order[i]].value("body_id", 0);
            // fps.json has no dye column of its own; a converted legacy file
            // may carry one, and FPS colours a pseudoatom only when it does.
            label.dye = positions[order[i]].value("dye", std::string());
        }
        out.push_back(label);
    }
    return out;
}

/* ------------------------------------------------------------------------
 * The table
 * ------------------------------------------------------------------------ */

int FPSResultTable::get_dof() const {
    int n_distances = 0;
    if (!rows.empty()) n_distances = static_cast<int>(rows[0].pairs.size());
    const int used = n_distances - 6 * (std::max(1, n_molecules) - 1);
    return used > 1 ? used : 1;
}

std::string FPSResultTable::get_json() const {
    nlohmann::json out;
    out["n_molecules"] = n_molecules;
    out["reference_number"] = reference_number;
    out["distance_type"] = distance_type;
    out["dof"] = get_dof();
    nlohmann::json list = nlohmann::json::array();
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const FPSResultRow& row = rows[i];
        nlohmann::json entry;
        entry["number"] = row.number;
        entry["chi2"] = row.chi2;
        entry["chi2_r"] = row.chi2 / static_cast<double>(get_dof());
        entry["chi2_bond"] = row.chi2_bond;
        entry["chi2_clash"] = row.chi2_clash;
        entry["converged"] = row.converged;
        entry["parent"] = row.parent;
        entry["method"] = row.method;
        entry["rmsd_previous"] = row.rmsd_previous;
        entry["rmsd_reference"] = row.rmsd_reference;
        entry["poses"] =
                nlohmann::json::parse(row.poses.empty() ? "[]" : row.poses);
        nlohmann::json pairs = nlohmann::json::array();
        for (std::size_t k = 0; k < row.pairs.size(); ++k) {
            pairs.push_back(nlohmann::json::parse(row.pairs[k].get_json()));
        }
        entry["pairs"] = pairs;
        list.push_back(entry);
    }
    out["rows"] = list;
    return out.dump();
}

FPSResultRow fps_result_row(const DockingResult& result, int number,
                            const std::string& method, int parent) {
    FPSResultRow row;
    row.number = number;
    // FPS's `E` excludes the clash term; this module's score includes it,
    // because the clash term is part of the objective the minimiser descended.
    row.chi2 = result.score - result.e_clash;
    row.chi2_bond = result.e_bond;
    row.chi2_clash = result.e_clash;
    row.converged = result.converged;
    row.parent = parent;
    row.method = method;
    row.poses = result.poses;
    row.pairs = result.pairs;
    return row;
}

FPSResultTable fps_bootstrap_table(const BootstrapResult& result,
                                   int n_molecules) {
    FPSResultTable table;
    table.n_molecules = n_molecules;
    table.reference_number = 1;  // the parent
    if (!result.truth.empty()) {
        table.distance_type = result.truth[0].distance_type;
    }

    FPSResultRow parent;
    parent.number = 1;
    parent.chi2 = result.parent_score;
    parent.method = "Docking";
    parent.parent = 0;
    parent.converged = true;
    parent.poses = result.parent_poses;
    parent.pairs = result.truth;
    table.rows.push_back(parent);

    for (std::size_t i = 0; i < result.replicas.size(); ++i) {
        const BootstrapReplica& replica = result.replicas[i];
        FPSResultRow row;
        row.number = static_cast<int>(i) + 2;
        row.chi2 = replica.score - replica.e_clash;
        row.chi2_bond = replica.e_bond;
        row.chi2_clash = replica.e_clash;
        row.converged = replica.converged;
        row.parent = 1;
        row.method = "ErrorEstimation";
        row.poses = replica.poses;
        row.pairs = replica.pairs;
        // Already measured against the parent by the bootstrap itself; filling
        // it here means the table is complete without a second assembly.
        row.rmsd_reference = replica.rmsd_to_parent;
        table.rows.push_back(row);
    }
    return table;
}

FPSResultTable add_rmsd_columns(const FPSResultTable& source,
                                const DockingAssembly& assembly,
                                bool fps_sign_convention) {
    FPSResultTable table = source;
    const FPSResultRow* reference = reference_row(table);
    const std::string reference_poses =
            reference != NULL ? reference->poses : std::string();
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        FPSResultRow& row = table.rows[i];
        row.rmsd_reference = pose_rmsd(assembly, row.poses, reference_poses,
                                       fps_sign_convention);
        // "Previous" is the previous row of *this table*, which is what FPS
        // means and is not the previous in simulation time.
        row.rmsd_previous = i == 0 ? -1.0
                                   : pose_rmsd(assembly, row.poses,
                                               table.rows[i - 1].poses,
                                               fps_sign_convention);
    }
    return table;
}

FPSResultTable add_best_fit(const FPSResultTable& source,
                            const DockingAssembly& assembly) {
    FPSResultTable table = source;
    const FPSResultRow* reference = reference_row(table);
    const std::string reference_poses =
            reference != NULL ? reference->poses : std::string();
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        FPSResultRow& row = table.rows[i];
        const IMP::algebra::Transformation3D fit =
                pose_superposition(assembly, row.poses, reference_poses);
        const std::pair<IMP::algebra::Vector3D, double> axis_angle =
                IMP::algebra::get_axis_and_angle(fit.get_rotation());
        row.best_fit_axis = axis_angle.first;
        row.best_fit_angle = axis_angle.second * 180.0 / IMP::algebra::PI;
        // PyMOL translates first and then rotates about the origin, so the
        // vector it must be given is R^-1 v, not v: R(x + R^-1 v) = Rx + v.
        row.best_fit_translation =
                fit.get_rotation().get_inverse().get_rotated(
                        fit.get_translation());
    }
    return table;
}

/* ------------------------------------------------------------------------
 * The writers
 * ------------------------------------------------------------------------ */

std::string write_fps_pymol_script(const std::string& directory,
                                   const FPSResultRow& row,
                                   const FPSMolecules& molecules,
                                   const FPSExportOptions& options) {
    internal::make_directory(directory);
    const std::string name = result_name(options, row.number);
    const std::string path = directory + "/" + name + ".pml";
    std::ofstream out = export_stream(path);

    out << "# Energy = " << fixed_text(row.chi2, 8)
        << (row.converged ? " (converged)" : " (not converged)") << "\n";
    write_loads(out, molecules);
    // Every molecule, molecule 0 included. FPS starts at 1, which is valid
    // only while its own engine normalises molecule 0 to the identity.
    const std::vector<PymolTransform> transforms =
            transforms_of(row.poses, molecules);
    write_pose(out, transforms, molecules, std::string(), options);

    if (options.best_fit) {
        out << "select all\n";
        write_best_fit(out, row, "sele", options);
    }
    if (options.save_pdb) {
        out << "select all\n";
        out << "save " << directory << "/" << name << ".pdb, sele\n";
    }

    std::vector<std::string> label_names;
    for (std::size_t i = 0; i < options.labels.size(); ++i) {
        const FPSLabelPosition& label = options.labels[i];
        std::size_t body = 0;
        for (std::size_t k = 0; k < transforms.size(); ++k) {
            if (static_cast<int>(k) == label.body_id) body = k;
        }
        IMP::algebra::Vector3D position = label.position;
        if (body < transforms.size()) {
            // `R*(lp - CM) + CM + T` -- the same arithmetic the `rotate`
            // /`translate` pair above performs on the molecule, so a
            // pseudoatom lands where its atoms do.
            const PymolTransform& t = transforms[body];
            position = t.rotation.get_rotated(label.position - t.origin) +
                       t.origin + t.translation;
        }
        out << "pseudoatom " << label.name << ", pos=[" << vector_text(position)
            << "]\n";
        out << "label " << label.name << ", \"" << label.name << "\"\n";
        out << "show spheres, " << label.name << "\n";
        if (label.dye == "Donor") out << "color green, " << label.name << "\n";
        if (label.dye == "Acceptor") out << "color red, " << label.name << "\n";
        label_names.push_back(label.name);
    }
    if (options.best_fit && !label_names.empty()) {
        out << "deselect\n";
        out << "select " << label_names[0];
        for (std::size_t i = 1; i < label_names.size(); ++i) {
            out << " + " << label_names[i];
        }
        out << "\n";
        write_best_fit(out, row, "sele", options);
    }
    out << "deselect\n";
    return path;
}

std::vector<std::string> write_fps_pymol_scripts(
        const std::string& directory, const FPSResultTable& table,
        const FPSMolecules& molecules, const FPSExportOptions& options) {
    std::vector<std::string> out;
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        out.push_back(write_fps_pymol_script(directory, table.rows[i],
                                             molecules, options));
    }
    return out;
}

namespace {

//! The header and `load` block Overlay and OverlayStates share verbatim.
void write_overlay_head(std::ostream& out, const FPSResultTable& table,
                        const FPSMolecules& molecules) {
    out << "# Overlay of structures";
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        out << (i == 0 ? " " : ", ") << table.rows[i].number;
    }
    out << "\n";
    write_loads(out, molecules);
}

//! `copy _tmp<i>` per molecule, then the pose -- the body of both overlays.
void write_overlay_copies(std::ostream& out, const FPSResultRow& row,
                          const FPSMolecules& molecules,
                          const FPSExportOptions& options) {
    const std::vector<PymolTransform> transforms =
            transforms_of(row.poses, molecules);
    for (std::size_t i = 0; i < molecules.size(); ++i) {
        const std::string tmp = "_tmp" + integer_text(static_cast<int>(i));
        out << "copy " << tmp << ", " << molecules[i].name << "\n";
        if (i >= transforms.size()) continue;
        out << "rotate [" << vector_text(transforms[i].axis) << "], "
            << fixed_text(transforms[i].angle_degrees, 3) << ", " << tmp
            << ", origin=[" << vector_text(transforms[i].origin) << "]"
            << camera(options) << "\n";
        out << "translate [" << vector_text(transforms[i].translation) << "], "
            << tmp << camera(options) << "\n";
    }
}

}  // namespace

std::string write_fps_overlay(const std::string& directory,
                              const FPSResultTable& table,
                              const FPSMolecules& molecules,
                              const FPSExportOptions& options) {
    internal::make_directory(directory);
    const std::string path = directory + "/" + options.prefix + "Overlay.pml";
    std::ofstream out = export_stream(path);
    write_overlay_head(out, table, molecules);
    for (std::size_t j = 0; j < table.rows.size(); ++j) {
        const FPSResultRow& row = table.rows[j];
        write_overlay_copies(out, row, molecules, options);
        // `select _tmp*`, not FPS's `select object _tmp*`: `object` is not a
        // selection operator, and the FPS spelling depends on the PyMOL
        // version parsing it leniently.
        out << "select _tmp*\n";
        out << "create _tmpjoin, sele\n";
        if (options.best_fit) write_best_fit(out, row, "_tmpjoin", options);
        out << "copy " << result_name(options, row.number) << ", _tmpjoin\n";
        out << "delete _tmp*\n";
    }
    for (std::size_t i = 0; i < molecules.size(); ++i) {
        out << "delete " << molecules[i].name << "\n";
    }
    return path;
}

std::string write_fps_overlay_states(const std::string& directory,
                                     const FPSResultTable& table,
                                     const FPSMolecules& molecules,
                                     const FPSExportOptions& options) {
    internal::make_directory(directory);
    const std::string path =
            directory + "/" + options.prefix + "OverlayStates.pml";
    // FPS writes the bare name `_tmp.pdb`, so the scratch file lands in
    // PyMOL's working directory and stays there. An absolute path in the
    // export directory, deleted at the end, is the same mechanism without the
    // two failure modes.
    const std::string tmp_pdb =
            options.tmp_pdb.empty() ? directory + "/_tmp.pdb" : options.tmp_pdb;
    std::ofstream out = export_stream(path);
    write_overlay_head(out, table, molecules);
    for (std::size_t j = 0; j < table.rows.size(); ++j) {
        const FPSResultRow& row = table.rows[j];
        write_overlay_copies(out, row, molecules, options);
        out << "select _tmp*\n";
        if (options.best_fit) write_best_fit(out, row, "sele", options);
        out << "save " << tmp_pdb << ", sele\n";
        out << "load " << tmp_pdb << ", " << options.prefix << "\n";
        // The state's identity, which FPS drops: states come out 1..k and the
        // only record of which result each one is is the header comment.
        out << "set_title " << options.prefix << ", "
            << integer_text(static_cast<int>(j) + 1) << ", \""
            << result_name(options, row.number) << "\"\n";
        out << "delete _tmp*\n";
    }
    for (std::size_t i = 0; i < molecules.size(); ++i) {
        out << "delete " << molecules[i].name << "\n";
    }
    out << "system rm -f " << tmp_pdb << "\n";
    return path;
}

std::string write_fps_r_table(const std::string& directory,
                              const FPSResultTable& table,
                              const FPSExportOptions& options) {
    internal::make_directory(directory);
    const std::string type =
            table.distance_type.empty() ? std::string("Rmp") : table.distance_type;
    const std::string path =
            directory + "/" + options.prefix + "Rtable_" + type + ".txt";
    std::vector<std::string> headers;
    const std::vector<std::string> names = distance_columns(table, &headers);

    std::ofstream out = export_stream(path);
    out << "Structure\tNumber";
    for (std::size_t i = 0; i < headers.size(); ++i) out << "\t" << headers[i];
    out << "\n";
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        const FPSResultRow& row = table.rows[i];
        out << result_name(options, row.number) << "\t" << row.number;
        for (std::size_t k = 0; k < names.size(); ++k) {
            out << "\t"
                << fixed_text(model_distance_of(row.pairs, names[k]), 3);
        }
        out << "\n";
    }
    return path;
}

std::string write_fps_chi2_table(const std::string& directory,
                                 const FPSResultTable& table,
                                 const FPSExportOptions& options) {
    internal::make_directory(directory);
    const std::string path =
            directory + "/" + options.prefix + "chi2table.txt";
    const FPSResultRow* reference = reference_row(table);
    const int reference_number = reference != NULL ? reference->number : 0;
    const double dof = static_cast<double>(table.get_dof());

    std::ofstream out = export_stream(path);
    out << "File\tchi2\tchi2_r\tchi2_bond\tchi2_clash\tConverged\tMethod"
           "\tParent\tRMSD vs previous\tRMSD vs " << reference_number << "\n";
    for (std::size_t i = 0; i < table.rows.size(); ++i) {
        const FPSResultRow& row = table.rows[i];
        out << result_name(options, row.number) << "\t"
            << fixed_text(row.chi2, 8) << "\t"
            << fixed_text(row.chi2 / dof, 8) << "\t"
            << fixed_text(row.chi2_bond, 4) << "\t"
            << fixed_text(row.chi2_clash, 4) << "\t"
            << (row.converged ? "True" : "False") << "\t" << row.method << "\t"
            << row.parent << "\t"
            // FPS writes the literal `---` for the first row's "previous"
            // column; a negative value here means the same thing.
            << (i == 0 || row.rmsd_previous < 0.0
                        ? std::string("---")
                        : fixed_text(row.rmsd_previous, 6))
            << "\t" << fixed_text(row.rmsd_reference, 6) << "\n";
    }
    return path;
}

std::vector<std::string> write_fps_exports(const std::string& directory,
                                           const FPSResultTable& table,
                                           const FPSMolecules& molecules,
                                           const FPSExportOptions& options) {
    std::vector<std::string> out =
            write_fps_pymol_scripts(directory, table, molecules, options);
    out.push_back(write_fps_overlay(directory, table, molecules, options));
    out.push_back(write_fps_overlay_states(directory, table, molecules, options));
    out.push_back(write_fps_r_table(directory, table, options));
    out.push_back(write_fps_chi2_table(directory, table, options));
    return out;
}

/* ------------------------------------------------------------------------
 * Filter mode
 * ------------------------------------------------------------------------ */

std::string write_fps_screening_r_table(
        const std::string& directory,
        const std::vector<ScreenedStructure>& structures,
        const FPSExportOptions& options) {
    internal::make_directory(directory);
    std::string type = "Rmp";
    std::vector<std::string> names, headers;
    for (std::size_t i = 0; i < structures.size(); ++i) {
        const PairDistances& pairs = structures[i].pairs;
        for (std::size_t k = 0; k < pairs.size(); ++k) {
            if (std::find(names.begin(), names.end(), pairs[k].name) !=
                names.end()) {
                continue;
            }
            names.push_back(pairs[k].name);
            headers.push_back(pairs[k].position1 + "_" + pairs[k].position2);
            if (names.size() == 1 && !pairs[k].distance_type.empty()) {
                type = pairs[k].distance_type;
            }
        }
    }
    const std::string path =
            directory + "/" + options.prefix + "Rtable_" + type + ".txt";
    std::ofstream out = export_stream(path);
    out << "File\tNumber";
    for (std::size_t i = 0; i < headers.size(); ++i) out << "\t" << headers[i];
    out << "\n";
    for (std::size_t i = 0; i < structures.size(); ++i) {
        // 1-based by default: FPS prints the 0-based array index here and the
        // 1-based `InternalNumber` everywhere else. Join on `File`.
        const int number = static_cast<int>(i) + (options.fps_filter_number ? 0 : 1);
        out << export_basename_of(structures[i].path) << "\t" << number;
        for (std::size_t k = 0; k < names.size(); ++k) {
            out << "\t"
                << fixed_text(model_distance_of(structures[i].pairs, names[k]),
                              3);
        }
        out << "\n";
    }
    return path;
}

std::string write_fps_screening_chi2_table(
        const std::string& directory,
        const std::vector<ScreenedStructure>& structures,
        const FPSExportOptions& options) {
    internal::make_directory(directory);
    const std::string path = directory + "/" + options.prefix + "chi2table.txt";
    std::ofstream out = export_stream(path);
    out << "File\tChi2r\tNaNs\tRefRMSD\t>1sigma\t>2sigma\t>3sigma\n";
    for (std::size_t i = 0; i < structures.size(); ++i) {
        const ScreenedStructure& s = structures[i];
        out << export_basename_of(s.path) << "\t" << fixed_text(s.chi2_r, 3) << "\t"
            << s.invalid_r << "\t" << fixed_text(s.ref_rmsd, 3) << "\t"
            << s.sigma1 << "\t" << s.sigma2 << "\t" << s.sigma3 << "\n";
    }
    return path;
}

IMPBFF_END_NAMESPACE
