/** \file Bin.cpp
 *  \brief The command line, compiled: one dispatcher, one function per sub.
 *
 *  See include/Bin.h for the why. The subs here are the ones whose subject
 *  is already the IMP-free core -- the same work the scripts in bin/ did,
 *  called straight into the library instead of through a Python shell, so
 *  they ship in the wheel and run with no IMP.
 *
 *  The grammar is CLI11 (include/internal/CLI11.h, vendored verbatim from
 *  github.com/CLIUtils/CLI11 v2.7.2, BSD-3): subcommands first-class, the
 *  help text generated, usage errors answered with exit code 2 -- the
 *  behaviour click gave the Python programs, in the library.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/Bin.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <IMP/bff/internal/CLI11.h>

#include <IMP/bff/Clustering.h>
#include <IMP/bff/RotamerLibrary.h>
#include <IMP/bff/StructureIO.h>
#include <IMP/bff/TrajectoryIO.h>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! A sub's own failure: printed as `imp_bff <sub>: <message>`, exit 1.
struct SubError : std::runtime_error {
    explicit SubError(const std::string& what) : std::runtime_error(what) {}
};

//! The calling sub, for the message prefix.
const char* g_sub = "imp_bff";

//! What the help says about the commands that stayed Python.
const char* FOOTER_TEXT =
    "The IMP- and RMF-bound modelling commands (flexfit, rmsd, select-pairs, "
    "openmm, ...) stay with the Python program bin/imp_bff of the conda "
    "package; they are not lighter for being compiled.";

//! A text file of one number per line (or whitespace-separated), as weights.
std::vector<double> read_weights(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) throw SubError(path + ": cannot read weights");
    std::vector<double> out;
    double v = 0.0;
    while (in >> v) out.push_back(v);
    return out;
}

//! View a std::vector's bytes as the (pointer, length) the writers take.
double* data_of(std::vector<double>& v) { return v.empty() ? nullptr : &v[0]; }

// ---------------------------------------------------------------------------
// pdb2cif
// ---------------------------------------------------------------------------

void run_pdb2cif(const std::string& pdb, const std::string& cif,
                 const std::string& probe_id) {
    convert_pdb_to_cif(pdb, cif, probe_id);
    std::cout << "pdb2cif: " << cif << "\n";
}

// ---------------------------------------------------------------------------
// traj2drot
// ---------------------------------------------------------------------------

const double LOSSLESS_TOL_A = 1e-4;
const double COMPACT_TOL_STEPS = 40.0;

//! Atom names, elements and residue names of a template PDB, in file order.
/*! The element column is authoritative when present -- bond perception
    depends on it, and guessing from a name calls CL chlorine or carbon
    depending on the day. */
void parse_template(const std::string& path,
                    std::vector<std::string>& names,
                    std::vector<std::string>& elements,
                    std::vector<std::string>& resnames) {
    const std::vector<PDBAtomRecord> records = read_pdb_records(path);
    for (std::size_t i = 0; i < records.size(); ++i) {
        names.push_back(records[i].atom_name);
        elements.push_back(records[i].element);
        resnames.push_back(records[i].res_name);
    }
    if (names.empty()) throw SubError(path + ": no ATOM/HETATM records");
}

//! A template PDB beside a `<stem>_cutoff<N>.<ext>` library, if there is one.
std::string template_for(const std::string& src) {
    const std::size_t slash = src.find_last_of("/\\");
    const std::string dir =
        slash == std::string::npos ? "" : src.substr(0, slash + 1);
    std::string name = slash == std::string::npos ? src : src.substr(slash + 1);
    const std::size_t dot = name.find_last_of('.');
    const std::string stem =
        dot == std::string::npos ? name : name.substr(0, dot);
    const std::size_t cut = stem.find("_cutoff");
    const std::string base = cut == std::string::npos ? stem : stem.substr(0, cut);
    const std::string candidate = dir + base + ".pdb";
    std::ifstream in(candidate.c_str());
    return in ? candidate : std::string();
}

//! Lower-cased suffix test, so `.BCIF` and `.Bcif` read as what they are.
bool ends_with_lower(const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) return false;
    std::string tail = s.substr(s.size() - suffix.size());
    for (std::size_t i = 0; i < tail.size(); ++i) {
        tail[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(tail[i])));
    }
    return tail == suffix;
}

//! `(n_frames, n_atoms, 3)` from a .bcif, .dcd or .drot trajectory.
/*! An .xtc needs mdtraj, which lives with the IMP-side program; the compiled
    sub refuses it with that message rather than half-reading one. */
std::vector<double> read_frames(const std::string& src, int n_atoms,
                                int* n_frames) {
    std::vector<double> xyz;
    if (ends_with_lower(src, ".drot.pto")) {
        const RotamerLibrary lib = read_drot(src);
        double* view = nullptr;
        int n = 0;
        lib.get_coords(&view, &n);
        xyz.assign(view, view + n);
        std::free(view);
        *n_frames = lib.n_rotamers;
        return xyz;
    }
    double* view = nullptr;
    int n_flat = 0;
    if (ends_with_lower(src, ".dcd")) {
        read_dcd(src, -1, &view, &n_flat);
    } else if (ends_with_lower(src, ".bcif")) {
        read_bcif_trajectory(src, n_atoms, "_rotamer_coord", &view, &n_flat);
    } else if (ends_with_lower(src, ".xtc")) {
        throw SubError("reading an .xtc needs mdtraj; it stays with the "
                       "IMP-side program bin/imp_bff");
    } else {
        read_trajectory(src, n_atoms, -1, &view, &n_flat);
    }
    xyz.assign(view, view + n_flat);
    std::free(view);
    if (xyz.empty() ||
        xyz.size() % (3 * static_cast<std::size_t>(n_atoms)) != 0) {
        std::ostringstream msg;
        msg << src << ": " << xyz.size()
            << " coordinates is not a whole number of frames of " << n_atoms
            << " atoms";
        throw SubError(msg.str());
    }
    *n_frames = static_cast<int>(xyz.size() / (3 * static_cast<std::size_t>(n_atoms)));
    return xyz;
}

//! Leader-cluster raw frames; replace them with the leaders, weight by population.
void leader_cluster(std::vector<double>& xyz, int* n_frames, double threshold,
                    std::vector<double>& weights) {
    const int n_atoms = *n_frames > 0
                            ? static_cast<int>(xyz.size()) / (3 * *n_frames)
                            : 0;
    int* leaders = nullptr;
    int n_leaders = 0;
    cluster_frames_leader(data_of(xyz), *n_frames, n_atoms, 3, threshold,
                          &leaders, &n_leaders);
    int* labels = nullptr;
    int n_labels = 0;
    assign_frames_to_clusters(data_of(xyz), *n_frames, n_atoms, 3, leaders,
                              n_leaders, &labels, &n_labels);
    std::vector<double> counts(static_cast<std::size_t>(n_leaders), 0.0);
    for (int i = 0; i < n_labels; ++i) {
        if (labels[i] >= 0 && labels[i] < n_leaders) counts[labels[i]] += 1.0;
    }
    std::vector<double> reduced;
    reduced.reserve(static_cast<std::size_t>(n_leaders) * 3);
    for (int i = 0; i < n_leaders; ++i) {
        const std::size_t f = static_cast<std::size_t>(leaders[i]) * 3;
        reduced.push_back(xyz[f]);
        reduced.push_back(xyz[f + 1]);
        reduced.push_back(xyz[f + 2]);
    }
    xyz.swap(reduced);
    weights.swap(counts);
    *n_frames = n_leaders;
    std::free(leaders);
    std::free(labels);
}

//! One trajectory to one `.drot.pto`, round-trip checked by default.
void traj2drot_convert(const std::string& src, const std::string& dst,
                       std::string top, const std::string& weights_path,
                       bool have_cluster, double cluster_a, bool have_grid,
                       double grid_a, double grid_deg, bool verify) {
    std::vector<std::string> names, elements, resnames;
    int n_frames = 0;
    std::vector<double> xyz;

    // A library as the source carries its own template; nothing is asked.
    if (ends_with_lower(src, ".drot.pto") && top.empty()) {
        const RotamerLibrary lib = read_drot(src);
        names = lib.atom_names;
        elements = lib.elements;
        resnames = lib.resnames;
        double* view = nullptr;
        int n = 0;
        lib.get_coords(&view, &n);
        xyz.assign(view, view + n);
        std::free(view);
        n_frames = lib.n_rotamers;
    } else {
        if (top.empty()) top = template_for(src);
        if (top.empty()) {
            throw SubError(src + ": no --top given and no template PDB "
                                 "beside it");
        }
        parse_template(top, names, elements, resnames);
        xyz = read_frames(src, static_cast<int>(names.size()), &n_frames);
    }
    std::vector<double> weights;
    if (!weights_path.empty()) weights = read_weights(weights_path);
    if (have_cluster) leader_cluster(xyz, &n_frames, cluster_a, weights);
    if (weights.empty()) weights.assign(static_cast<std::size_t>(n_frames), 1.0);
    if (weights.size() != static_cast<std::size_t>(n_frames)) {
        std::ostringstream msg;
        msg << src << ": " << weights.size() << " weights for " << n_frames
            << " frames -- pass --cluster to build weights from populations";
        throw SubError(msg.str());
    }

    DrotEncoding encoding;
    if (have_grid) {
        encoding.lossless = false;
        encoding.grid_a = grid_a;
        encoding.grid_deg = grid_deg;
    }
    write_drot(dst, data_of(xyz), static_cast<int>(xyz.size()), names, elements,
               resnames, data_of(weights), static_cast<int>(weights.size()),
               encoding);

    std::ostringstream line;
    line << dst << ": " << n_frames << " rotamers, " << names.size() << " atoms ("
         << (have_grid ? "compact grid" : "lossless") << ")";
    if (verify) {
        const RotamerLibrary lib = read_drot(dst);
        double* view = nullptr;
        int n = 0;
        lib.get_coords(&view, &n);
        std::vector<double> back(view, view + n);
        std::free(view);
        if (back.size() != xyz.size()) {
            std::ostringstream msg;
            msg << dst << ": read back " << back.size() << " coordinates, wrote "
                << xyz.size();
            throw SubError(msg.str());
        }
        double err = 0.0;
        for (std::size_t i = 0; i < back.size(); ++i) {
            err = (std::max)(err, std::fabs(back[i] - xyz[i]));
        }
        const double tol = have_grid ? COMPACT_TOL_STEPS * grid_a : LOSSLESS_TOL_A;
        if (err > tol) {
            std::ostringstream msg;
            msg << dst << ": round trip differs by " << err << " A (tolerance "
                << tol << ")";
            throw SubError(msg.str());
        }
        if (lib.atom_names != names) {
            throw SubError(dst + ": atom names did not survive the round trip");
        }
        line << ", verified to " << err << " A";
    }
    std::cout << "traj2drot: " << line.str() << "\n";
}

//! Several one-library containers into one family container.
/*! A container operation only: the payloads cross verbatim under names
    prefixed `<library>/`, and a `drot.catalog` object lists what is inside.
    Nothing is re-encoded, so bundling is a copy and the result holds exactly
    the bytes its parts did. */
void traj2drot_bundle(const std::vector<std::string>& sources,
                      const std::string& out) {
    if (sources.empty()) {
        throw SubError("--bundle needs at least one source container");
    }
    std::vector<std::string> names;
    for (std::size_t i = 0; i < sources.size(); ++i) {
        const std::size_t slash = sources[i].find_last_of("/\\");
        std::string name = slash == std::string::npos
                               ? sources[i]
                               : sources[i].substr(slash + 1);
        const std::size_t drot = name.find(".drot");
        if (drot != std::string::npos) name = name.substr(0, drot);
        names.push_back(name);
    }
    for (std::size_t i = 0; i < names.size(); ++i) {
        for (std::size_t j = i + 1; j < names.size(); ++j) {
            if (names[i] == names[j]) {
                throw SubError("two sources would take the same library name (" +
                               names[i] + ")");
            }
        }
    }
    write_drot_bundle(sources, names, out);
    const std::vector<std::string> listed = drot_catalog(out);
    std::vector<std::string> sorted_listed = listed, sorted_names = names;
    std::sort(sorted_listed.begin(), sorted_listed.end());
    std::sort(sorted_names.begin(), sorted_names.end());
    if (sorted_listed != sorted_names) {
        throw SubError(out + ": the catalog does not list what went in");
    }
    std::cout << "traj2drot: " << out << ": " << names.size()
              << " libraries\n";
}

//! The work, with the parsed values.
void run_traj2drot(const std::vector<std::string>& positional,
                   const std::string& top, const std::string& weights_path,
                   const std::string& bundle_out, double cluster_a,
                   double grid_a, double grid_deg, bool no_verify) {
    if (!bundle_out.empty()) {
        traj2drot_bundle(positional, bundle_out);
        return;
    }
    if (positional.size() < 2) {
        throw SubError("give src and dst (the bulk --all path stays with the "
                       "Python program)");
    }
    traj2drot_convert(positional[0], positional[1], top, weights_path,
                      cluster_a >= 0.0, cluster_a, grid_a >= 0.0, grid_a,
                      grid_deg, !no_verify);
}

}  // namespace

int bin_main(const std::vector<std::string>& args) {
    // The (argc, argv) overload: it skips the program name and reads the
    // words in their natural order -- the vector overload expects them
    // reversed, which nothing that arrives from Python is.
    std::vector<char*> line;
    line.push_back(const_cast<char*>("imp_bff"));
    for (std::size_t i = 0; i < args.size(); ++i) {
        line.push_back(const_cast<char*>(args[i].c_str()));
    }

    CLI::App app("imp_bff -- the command line, compiled: the bin/ scripts "
                 "whose subject is the IMP-free core");
    app.require_subcommand(1);

    int rc = 0;
    CLI::App* pdb2cif = app.add_subcommand(
        "pdb2cif", "convert a probe PDB to an mmCIF structure file");
    std::string pdb, cif, probe_id;
    pdb2cif->add_option("pdb", pdb, "input PDB")->required();
    pdb2cif->add_option("cif", cif, "output mmCIF")->required();
    pdb2cif->add_option("--probe-id", probe_id,
                        "chemical component id to record; empty takes the stem")
        ->capture_default_str();

    CLI::App* traj2drot = app.add_subcommand(
        "traj2drot",
        "a trajectory (.bcif/.dcd) to a .drot.pto library, and bundling");
    // One positional list, whatever the mode: convert reads [src, dst],
    // --bundle reads every word as a source container.
    std::vector<std::string> trajectory;
    std::string top, weights_path, bundle_out;
    double cluster_a = -1.0, grid_a = -1.0, grid_deg = 0.01;
    bool no_verify = false;
    traj2drot
        ->add_option("trajectory", trajectory,
                     "a .bcif, .dcd or .drot trajectory; with --bundle, the "
                     "source containers")
        ->take_all();
    traj2drot
        ->add_option("--top", top,
                     "template PDB: atom names, elements, residue names. "
                     "Defaults to <stem>.pdb beside the source")
        ->capture_default_str();
    traj2drot
        ->add_option("--weights", weights_path,
                     "one weight per frame, as the shipped "
                     "<stem>_cutoff<N>_weights.txt")
        ->capture_default_str();
    traj2drot
        ->add_option("--cluster", cluster_a,
                     "leader-cluster the frames at this RMSD first; the "
                     "leaders become the rotamers and the cluster populations "
                     "the weights")
        ->capture_default_str();
    traj2drot
        ->add_option("--grid", grid_a,
                     "compact rung: int16 grids at this step for base "
                     "coordinates and bond lengths. Omit for the lossless "
                     "float32 default -- the pins depend on it")
        ->capture_default_str();
    traj2drot->add_option("--grid-deg", grid_deg,
                          "angle step of the compact rung")
        ->capture_default_str();
    traj2drot
        ->add_option("--bundle", bundle_out,
                     "bundle the given .drot.pto containers into one family "
                     "container, each filed under its stem")
        ->capture_default_str();
    traj2drot->add_flag("--no-verify", no_verify, "skip the round-trip check");

    // The IMP- and RMF-bound modelling commands -- flexfit, rmsd,
    // select-pairs, openmm -- stay with the Python program bin/imp_bff of
    // the conda package; they are not lighter for being compiled.
    app.add_subcommand("help", "print help for the compiled subs")
        ->callback([&] {
            std::cout << "usage: imp_bff <sub> [options...]\n\n";
            for (const auto& sub : bin_subs()) {
                std::cout << "  " << sub.first << "\t" << sub.second << "\n";
            }
            std::cout << "\n" << FOOTER_TEXT << "\n" << std::flush;
            rc = 0;
        });

    app.footer(FOOTER_TEXT);

    g_sub = "imp_bff";
    try {
        pdb2cif->callback([&] {
            g_sub = "pdb2cif";
            run_pdb2cif(pdb, cif, probe_id);
        });
        traj2drot->callback([&] {
            g_sub = "traj2drot";
            run_traj2drot(trajectory, top, weights_path, bundle_out,
                          cluster_a, grid_a, grid_deg, no_verify);
        });
        app.parse(static_cast<int>(line.size()), line.data());
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    } catch (const SubError& e) {
        std::cerr << "imp_bff " << g_sub << ": " << e.what() << "\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "imp_bff " << g_sub << ": " << e.what() << "\n";
        return 1;
    }
    return rc;
}

std::vector<std::pair<std::string, std::string> > bin_subs() {
    std::vector<std::pair<std::string, std::string> > out;
    out.push_back(std::make_pair(
        std::string("pdb2cif"),
        std::string("convert a probe PDB to an mmCIF structure file")));
    out.push_back(
        std::make_pair(std::string("traj2drot"),
                       std::string("a trajectory (.bcif/.dcd) to a .drot.pto "
                                   "library, and bundling")));
    return out;
}

IMPBFF_END_NAMESPACE
