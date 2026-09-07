/**
 * \file StructureIO.cpp
 * \brief Structures in and out: PDB, MOL2, mmCIF, and the tables beside them.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/StructureIO.h>
#include <IMP/bff/internal/PdbFrames.h>
#include <iomanip>

#include <IMP/bff/internal/Cif.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/Text.h>

#include <IMP/bff/internal/json.h>
#include <IMP/bff/Base.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <set>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE


// Named, not anonymous: IMP compiles this module as one translation unit.
namespace structio {
using IMP::bff::internal::trimmed;

std::string field(const std::string& line, std::size_t at, std::size_t n) {
    if (at >= line.size()) return std::string();
    return trimmed(line.substr(at, std::min(n, line.size() - at)));
}

//! The Tripos atom type of an element; the element itself when unlisted.
std::string tripos_type(const std::string& element) {
    if (element == "C") return "C.3";
    if (element == "N") return "N.3";
    if (element == "O") return "O.3";
    if (element == "S") return "S.3";
    if (element == "P") return "P.3";
    if (element == "H") return "H";
    if (element == "F") return "F";
    if (element == "Cl") return "Cl";
    if (element == "Br") return "Br";
    if (element == "I") return "I";
    return element;
}

//! The element a PDB atom name implies: its first alphabetic character.
std::string element_from_name(const std::string& name) {
    for (std::size_t i = 0; i < name.size(); ++i) {
        if (std::isalpha(static_cast<unsigned char>(name[i]))) {
            return std::string(1, static_cast<char>(
                                          std::toupper(name[i])));
        }
    }
    return "C";
}

std::string stem_of(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    std::string base = slash == std::string::npos ? path : path.substr(slash + 1);
    const std::size_t dot = base.find_last_of('.');
    return dot == std::string::npos ? base : base.substr(0, dot);
}

}  // namespace structio

// --------------------------------------------------------------------------
// MOL2
// --------------------------------------------------------------------------

std::vector<AtomBond> parse_conect_bonds(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("Cannot read " << path, IOException);
    std::set<AtomBond> bonds;
    std::string line;
    while (std::getline(in, line)) {
        if (structio::field(line, 0, 6) != "CONECT") continue;
        std::istringstream fields(line);
        std::string tag;
        int a = 0;
        if (!(fields >> tag >> a)) continue;
        int b = 0;
        while (fields >> b) {
            if (a != b) bonds.insert(AtomBond(a, b));
        }
    }
    return std::vector<AtomBond>(bonds.begin(), bonds.end());
}

std::vector<AtomBond> infer_bonds(const std::vector<PDBAtomRecord>& atoms,
                                  double tolerance) {
    // Sorted by serial: the pair order the writer numbers bonds by follows
    // from it.
    std::vector<PDBAtomRecord> sorted = atoms;
    std::sort(sorted.begin(), sorted.end(),
              [](const PDBAtomRecord& l, const PDBAtomRecord& r) {
                  return l.serial < r.serial;
              });
    std::set<AtomBond> bonds;
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        const double ra = IMP::bff::covalent_radius(sorted[i].element);
        for (std::size_t j = i + 1; j < sorted.size(); ++j) {
            const double rb = IMP::bff::covalent_radius(sorted[j].element);
            const double dx = sorted[i].x - sorted[j].x;
            const double dy = sorted[i].y - sorted[j].y;
            const double dz = sorted[i].z - sorted[j].z;
            if (std::sqrt(dx * dx + dy * dy + dz * dz) <= ra + rb + tolerance) {
                bonds.insert(AtomBond(sorted[i].serial, sorted[j].serial));
            }
        }
    }
    return std::vector<AtomBond>(bonds.begin(), bonds.end());
}

void write_mol2(const std::string& path,
                const std::vector<PDBAtomRecord>& atoms,
                const std::vector<AtomBond>& bonds,
                const std::string& mol_name) {
    std::vector<PDBAtomRecord> sorted = atoms;
    std::sort(sorted.begin(), sorted.end(),
              [](const PDBAtomRecord& l, const PDBAtomRecord& r) {
                  return l.serial < r.serial;
              });
    std::map<int, int> index;  // serial -> 1-based MOL2 atom index
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        index[sorted[i].serial] = static_cast<int>(i) + 1;
    }
    std::vector<AtomBond> sorted_bonds = bonds;
    std::sort(sorted_bonds.begin(), sorted_bonds.end());

    std::ofstream out(path.c_str());
    if (!out) IMP_THROW("Cannot write " << path, IOException);
    out << "@<TRIPOS>MOLECULE\n"
        << mol_name << "\n"
        << sorted.size() << " " << sorted_bonds.size() << " 0 0 0\n"
        << "SMALL\n"
        << "NO_CHARGES\n"
        << "\n"
        << "@<TRIPOS>ATOM\n";
    char buffer[256];
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        const PDBAtomRecord& a = sorted[i];
        std::snprintf(buffer, sizeof(buffer),
                      "%6d %-8s %10.4f %10.4f %10.4f %-8s 1 %s 0.0000",
                      index[a.serial], a.atom_name.c_str(), a.x, a.y, a.z,
                      structio::tripos_type(a.element).c_str(),
                      a.res_name.c_str());
        out << buffer << "\n";
    }
    out << "@<TRIPOS>BOND\n";
    for (std::size_t i = 0; i < sorted_bonds.size(); ++i) {
        std::snprintf(buffer, sizeof(buffer), "%6d %6d %6d 1",
                      static_cast<int>(i) + 1, index[sorted_bonds[i].a],
                      index[sorted_bonds[i].b]);
        out << buffer << "\n";
    }
}

// --------------------------------------------------------------------------
// PMI stat files
// --------------------------------------------------------------------------

namespace structio {

//! Integer column keys PMI uses in the per-frame stat dicts.
const int TOTAL_SCORE_KEY = 1;
const int NFRAME_KEY = 4;

//! A flat `{key: value}` Python-literal line, as far as this needs it.
/*!
    Only integer keys and numeric values are kept: those are the two columns a
    progress bar wants, and a stat line's other values are strings and nested
    calls. Anything that does not parse yields an empty map, which is the
    contract `ast.literal_eval` inside a `try` had.
*/
bool literal_int_map(const std::string& line, std::map<int, double>& out,
                     std::map<int, std::string>& names) {
    const std::string s = trimmed(line);
    if (s.size() < 2 || s[0] != '{' || s[s.size() - 1] != '}') return false;
    std::size_t i = 1;
    const std::size_t end = s.size() - 1;
    while (i < end) {
        while (i < end && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
        if (i >= end) break;
        // key
        const std::size_t key_start = i;
        while (i < end && s[i] != ':') ++i;
        if (i >= end) return false;
        const std::string key_text = trimmed(s.substr(key_start, i - key_start));
        ++i;  // past ':'
        while (i < end && std::isspace(static_cast<unsigned char>(s[i]))) ++i;

        // value: a quoted string, or everything up to the next top-level comma
        std::string value_text;
        bool is_string = false;
        if (i < end && (s[i] == '\'' || s[i] == '"')) {
            const char quote = s[i++];
            const std::size_t value_start = i;
            while (i < end && s[i] != quote) {
                if (s[i] == '\\') ++i;
                ++i;
            }
            if (i >= end) return false;
            value_text = s.substr(value_start, i - value_start);
            is_string = true;
            ++i;  // past the closing quote
            while (i < end && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
        } else {
            int depth = 0;
            const std::size_t value_start = i;
            while (i < end) {
                const char c = s[i];
                if (c == '(' || c == '[' || c == '{') ++depth;
                else if (c == ')' || c == ']' || c == '}') --depth;
                else if (c == ',' && depth == 0) break;
                ++i;
            }
            value_text = trimmed(s.substr(value_start, i - value_start));
            // A call such as `environ(...)` is what the header line carries,
            // and it is why the header does not parse. Refuse the whole line,
            // as `literal_eval` does.
            if (value_text.find('(') != std::string::npos) return false;
        }
        if (i < end) {
            if (s[i] != ',') return false;
            ++i;
        }

        char* stop = nullptr;
        const long key = std::strtol(key_text.c_str(), &stop, 10);
        if (stop != key_text.c_str() && *stop == '\0') {
            if (is_string) {
                names[static_cast<int>(key)] = value_text;
            } else {
                char* vstop = nullptr;
                const double v = std::strtod(value_text.c_str(), &vstop);
                if (vstop != value_text.c_str() && *vstop == '\0') {
                    out[static_cast<int>(key)] = v;
                }
            }
        }
    }
    return true;
}

}  // namespace structio

int count_frames(const std::string& stat_path) {
    std::ifstream in(stat_path.c_str());
    if (!in) return 0;
    int n = 0;
    std::string line;
    while (std::getline(in, line)) ++n;
    return n > 0 ? n - 1 : 0;
}

void read_score_series(const std::string& stat_path, double** out_view,
                       int* n_out_view) {
    std::vector<double> packed;
    std::ifstream in(stat_path.c_str());
    if (in) {
        int score_key = structio::TOTAL_SCORE_KEY;
        int frame_key = structio::NFRAME_KEY;
        std::string line;
        int index = 0;
        while (std::getline(in, line)) {
            line = structio::trimmed(line);
            if (!line.empty() && line[0] == '{') {
                std::map<int, double> values;
                std::map<int, std::string> names;
                const bool parsed =
                        structio::literal_int_map(line, values, names);
                if (index == 0) {
                    // The header names the columns, in case PMI renumbers them.
                    if (parsed) {
                        for (std::map<int, std::string>::const_iterator it =
                                     names.begin();
                             it != names.end(); ++it) {
                            if (it->second == "Total_Score") score_key = it->first;
                            else if (it->second == "MonteCarlo_Nframe") {
                                frame_key = it->first;
                            }
                        }
                    }
                    ++index;
                    continue;
                }
                ++index;
                if (!parsed || values.count(score_key) == 0) continue;
                const double frame =
                        values.count(frame_key)
                                ? values.find(frame_key)->second
                                : static_cast<double>(packed.size() / 2);
                packed.push_back(frame);
                packed.push_back(values.find(score_key)->second);
            } else {
                ++index;
                // A plain `frame,score` convergence CSV, which a minimisation
                // writes; its header row simply fails to parse.
                const std::size_t comma = line.find(',');
                if (comma == std::string::npos) continue;
                const std::string a = structio::trimmed(line.substr(0, comma));
                const std::string b = structio::trimmed(line.substr(comma + 1));
                char* stop_a = nullptr;
                char* stop_b = nullptr;
                const double frame = std::strtod(a.c_str(), &stop_a);
                const double score = std::strtod(b.c_str(), &stop_b);
                if (stop_a == a.c_str() || *stop_a != '\0') continue;
                if (stop_b == b.c_str() || *stop_b != '\0') continue;
                packed.push_back(frame);
                packed.push_back(score);
            }
        }
    }
    internal::copy_to_view(packed, out_view, n_out_view);
}

// --------------------------------------------------------------------------
// coordinates
// --------------------------------------------------------------------------

namespace structio {

std::vector<double> transformed(const std::vector<double>& coords,
                                const std::vector<double>& t) {
    if (t.empty()) return coords;
    std::vector<double> out = coords;
    const std::size_t n = coords.size() / 3;
    if (t.size() == 3) {
        for (std::size_t i = 0; i < n; ++i) {
            for (int k = 0; k < 3; ++k) out[3 * i + k] += t[k];
        }
        return out;
    }
    if (t.size() == 9 || t.size() == 16) {
        const std::size_t stride = t.size() == 9 ? 3 : 4;
        for (std::size_t i = 0; i < n; ++i) {
            double v[3];
            for (int r = 0; r < 3; ++r) {
                v[r] = 0.0;
                for (int c = 0; c < 3; ++c) {
                    v[r] += t[r * stride + c] * coords[3 * i + c];
                }
                if (stride == 4) v[r] += t[r * stride + 3];
            }
            for (int k = 0; k < 3; ++k) out[3 * i + k] = v[k];
        }
        return out;
    }
    IMP_THROW("a transform is 3 values (a translation), 9 (a rotation) or 16 "
              "(a homogeneous 4x4), not " << t.size(),
              ValueException);
}

}  // namespace structio

void apply_transform(const std::vector<double>& coords,
                     const std::vector<double>& transform, double** out_view,
                     int* n_out_view) {
    internal::copy_to_view(structio::transformed(coords, transform), out_view,
                           n_out_view);
}

void write_pdb(const std::vector<double>& coords, const std::string& path,
               const std::string& chain, const std::string& res_name,
               const std::vector<double>& transform, int model_index) {
    if (coords.size() % 3 != 0) {
        IMP_THROW("coordinates come three per atom, not " << coords.size(),
                  ValueException);
    }
    const std::vector<double> moved = structio::transformed(coords, transform);
    std::ofstream out(path.c_str());
    if (!out) IMP_THROW("Cannot write " << path, IOException);
    char buffer[128];
    std::snprintf(buffer, sizeof(buffer), "MODEL     %4d", model_index);
    out << buffer << "\n";
    for (std::size_t i = 0; i < moved.size() / 3; ++i) {
        std::snprintf(buffer, sizeof(buffer),
                      "ATOM  %5d  CA %-3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
                      static_cast<int>(i) + 1, res_name.c_str(), chain.c_str(),
                      static_cast<int>(i) + 1, moved[3 * i + 0],
                      moved[3 * i + 1], moved[3 * i + 2], 1.0, 0.0);
        out << buffer << "\n";
    }
    out << "ENDMDL\nEND\n";
}

void load_structure(const std::string& path, double** out_view,
                    int* n_out_view) {
    // the first model, without water: IMP's NonWaterPDBSelector, in the core
    const std::vector<ProteinFrame> frames =
            internal::read_pdb_frames(path, internal::PDB_NON_WATER, true, 1);
    internal::copy_to_view(frames[0].coords, out_view, n_out_view);
}

double get_rmsd(const std::vector<double>& coords_a,
                    const std::vector<double>& coords_b,
                    const std::vector<int>& selection_mask, bool superpose) {
    if (coords_a.size() != coords_b.size()) {
        IMP_THROW("Shape mismatch: " << coords_a.size() << " vs "
                                     << coords_b.size(),
                  ValueException);
    }
    const std::size_t n = coords_a.size() / 3;
    std::vector<std::size_t> rows;
    if (selection_mask.empty()) {
        for (std::size_t i = 0; i < n; ++i) rows.push_back(i);
    } else {
        if (selection_mask.size() != n) {
            IMP_THROW("the selection mask has " << selection_mask.size()
                                                << " entries for " << n
                                                << " atoms",
                      ValueException);
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (selection_mask[i]) rows.push_back(i);
        }
    }
    if (rows.empty()) return 0.0;

    std::vector<double> a = coords_a;
    if (superpose) {
        Eigen::MatrixXd sel_a(rows.size(), 3), sel_b(rows.size(), 3);
        for (std::size_t i = 0; i < rows.size(); ++i) {
            for (int k = 0; k < 3; ++k) {
                sel_a(i, k) = coords_a[3 * rows[i] + k];
                sel_b(i, k) = coords_b[3 * rows[i] + k];
            }
        }
        const Eigen::Vector3d ca = sel_a.colwise().mean();
        const Eigen::Vector3d cb = sel_b.colwise().mean();
        const Eigen::Matrix3d cov =
                (sel_a.rowwise() - ca.transpose()).transpose() *
                (sel_b.rowwise() - cb.transpose());
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(
                cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3d d = Eigen::Matrix3d::Identity();
        const double det =
                (svd.matrixU() * svd.matrixV().transpose()).determinant();
        // The reflection guard: without it a near-degenerate covariance
        // superposes a structure onto its mirror image and reports a
        // plausible RMSD for a chirality flip.
        d(2, 2) = det < 0.0 ? -1.0 : 1.0;
        const Eigen::Matrix3d r = svd.matrixU() * d * svd.matrixV().transpose();
        // The rotation is fitted on the selection and applied to *all* of
        // coords_a: that is what makes a selection a reference rather than a
        // crop.
        for (std::size_t i = 0; i < n; ++i) {
            Eigen::Vector3d v(coords_a[3 * i + 0] - ca[0],
                              coords_a[3 * i + 1] - ca[1],
                              coords_a[3 * i + 2] - ca[2]);
            const Eigen::Vector3d moved = r.transpose() * v + cb;
            for (int k = 0; k < 3; ++k) a[3 * i + k] = moved[k];
        }
    }

    double total = 0.0;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        for (int k = 0; k < 3; ++k) {
            const double d = a[3 * rows[i] + k] - coords_b[3 * rows[i] + k];
            total += d * d;
        }
    }
    return std::sqrt(total / rows.size());
}

namespace {
std::ofstream open_for_write(const std::string& path) {
    std::ofstream out(path.c_str());
    if (!out) IMP_THROW("cannot write " << path, IOException);
    return out;
}
}  // namespace

void write_points_xyz(const std::string& path, const std::vector<double>& points,
                      const std::string& element, const std::string& comment) {
    std::ofstream out = open_for_write(path);
    const std::size_t n = points.size() / 4;
    out << n << "\n" << comment << "\n";
    out << std::fixed << std::setprecision(4);
    for (std::size_t i = 0; i < n; ++i) {
        out << element << " " << points[4 * i + 0] << " " << points[4 * i + 1]
            << " " << points[4 * i + 2] << " " << points[4 * i + 3] << "\n";
    }
}

void write_points_pqr(const std::string& path, const std::vector<double>& points,
                      double radius) {
    std::ofstream out = open_for_write(path);
    const std::size_t n = points.size() / 4;
    out << std::fixed;
    for (std::size_t i = 0; i < n; ++i) {
        // The PDB serial and residue-number fields are five and four digits;
        // a rastered volume overruns both, so they wrap rather than widen the
        // record and break every column after them.
        out << "ATOM  " << std::setw(5) << (int) (i % 100000 + 1) << "  AV  AV  "
            << std::setw(4) << (int) (i % 10000 + 1) << "    "
            << std::setw(8) << std::setprecision(3) << points[4 * i + 0]
            << std::setw(8) << std::setprecision(3) << points[4 * i + 1]
            << std::setw(8) << std::setprecision(3) << points[4 * i + 2]
            << std::setw(8) << std::setprecision(4) << points[4 * i + 3]
            << std::setw(7) << std::setprecision(3) << radius << "\n";
    }
    out << "END\n";
}

void write_opendx(const std::string& path, const std::vector<double>& density,
                  int nx, int ny, int nz, const std::vector<double>& origin,
                  double spacing) {
    const std::size_t want = (std::size_t) std::max(0, nx) *
                             (std::size_t) std::max(0, ny) *
                             (std::size_t) std::max(0, nz);
    if (density.size() != want) {
        IMP_THROW("density has " << density.size() << " values for a "
                                 << nx << "x" << ny << "x" << nz << " grid",
                  ValueException);
    }
    if (origin.size() != 3) {
        IMP_THROW("origin must be three coordinates", ValueException);
    }
    std::ofstream out = open_for_write(path);
    out << "# OpenDX density written by IMP.bff\n"
        << "object 1 class gridpositions counts " << nx << " " << ny << " "
        << nz << "\n";
    out << std::fixed << std::setprecision(6);
    out << "origin " << origin[0] << " " << origin[1] << " " << origin[2] << "\n"
        << "delta " << spacing << " 0 0\n"
        << "delta 0 " << spacing << " 0\n"
        << "delta 0 0 " << spacing << "\n"
        << "object 2 class gridconnections counts " << nx << " " << ny << " "
        << nz << "\n"
        << "object 3 class array type double rank 0 items " << want
        << " data follows\n";
    // OpenDX runs z fastest within y within x -- the same nesting a C-order
    // (nx, ny, nz) array has, so the values go out in the order they are in.
    int per_line = 0;
    for (std::size_t i = 0; i < want; ++i) {
        out << density[i];
        out << (++per_line % 3 == 0 ? "\n" : " ");
    }
    if (per_line % 3 != 0) out << "\n";
    out << "attribute \"dep\" string \"positions\"\n"
        << "object \"density\" class field\n"
        << "component \"positions\" value 1\n"
        << "component \"connections\" value 2\n"
        << "component \"data\" value 3\n";
}

void convert_pdb_to_cif(const std::string& pdb_path,
                        const std::string& cif_path,
                        const std::string& probe_id) {
    // Every ATOM/HETATM of the first model (IMP's AllPDBSelector, because dyes
    // often carry HETATM records and unusual residue names and a NonWater
    // selector drops them), read by the core (internal/PdbFrames.h).
    const ProteinFrame atoms =
            internal::read_pdb_frames(pdb_path, internal::PDB_ALL, true, 1)[0];

    std::ofstream out(cif_path.c_str());
    if (!out) IMP_THROW("Cannot write " << cif_path, IOException);

    // The same writer every other category in this module goes through
    // (`internal/Cif.h`), so an atom or component name with a space in
    // it is quoted here as it is there. IMP has no C++ CIF writer to borrow;
    // it has the reader, and that is what reads this back.
    internal::CifWriter w(out);
    w.start_block(probe_id.empty() ? structio::stem_of(pdb_path) : probe_id);

    const std::vector<std::string> cols = {
            "group_PDB", "id", "type_symbol", "label_atom_id", "label_comp_id",
            "label_asym_id", "label_entity_id", "label_seq_id",
            "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv"};
    std::vector<std::vector<std::string> > rows;
    rows.reserve(atoms.get_n_atoms());
    char coord[32];
    for (int i = 0; i < atoms.get_n_atoms(); ++i) {
        const double xyz[3] = {atoms.coords[3 * i], atoms.coords[3 * i + 1],
                               atoms.coords[3 * i + 2]};
        std::string name = structio::trimmed(atoms.atom_types[i]);
        // every ATOM/HETATM line has a residue (`UNK` when the field is empty)
        const std::string comp = atoms.resnames[i];
        std::vector<std::string> row = {
                "HETATM", std::to_string(i + 1),
                structio::element_from_name(name), name, comp, "A", "1", "1"};
        for (int k = 0; k < 3; ++k) {
            std::snprintf(coord, sizeof(coord), "%.3f", xyz[k]);
            row.push_back(coord);
        }
        row.push_back("1.000");
        row.push_back("0.000");
        rows.push_back(row);
    }
    w.write_loop("_atom_site", cols, rows);
}

// --------------------------------------------------------------------------
// the tables beside a structure
// --------------------------------------------------------------------------

std::vector<CrossLink> read_xlink_table(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("Cannot read " << path, IOException);
    std::vector<CrossLink> out;
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (first) { first = false; continue; }  // the header row
        line = structio::trimmed(line);
        std::vector<std::string> parts;
        std::size_t at = 0;
        while (true) {
            const std::size_t comma = line.find(',', at);
            parts.push_back(line.substr(at, comma == std::string::npos
                                                    ? std::string::npos
                                                    : comma - at));
            if (comma == std::string::npos) break;
            at = comma + 1;
        }
        if (parts.size() != 4) continue;
        char* stop_1 = nullptr;
        char* stop_2 = nullptr;
        const std::string r1 = structio::trimmed(parts[1]);
        const std::string r2 = structio::trimmed(parts[3]);
        const long residue_1 = std::strtol(r1.c_str(), &stop_1, 10);
        const long residue_2 = std::strtol(r2.c_str(), &stop_2, 10);
        if (stop_1 == r1.c_str() || *stop_1 != '\0') continue;
        if (stop_2 == r2.c_str() || *stop_2 != '\0') continue;
        CrossLink x;
        x.protein_1 = parts[0];
        x.protein_2 = parts[2];
        x.residue_1 = static_cast<int>(residue_1);
        x.residue_2 = static_cast<int>(residue_2);
        out.push_back(x);
    }
    return out;
}

IMPBFF_END_NAMESPACE
