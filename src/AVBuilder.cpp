/**
 * \file AVBuilder.cpp
 * \brief Building an accessible volume: the two front doors, and the PDB read.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/AVBuilder.h>

#include <IMP/bff/AV.h>
#include <IMP/bff/PathMap.h>
#include <IMP/bff/PathMapHeader.h>
#include <IMP/bff/StripMask.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/json.h>

#include <IMP/algebra/Vector3D.h>
#include <IMP/atom/Atom.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Selection.h>
#include <IMP/atom/pdb.h>
#include <IMP/core/XYZR.h>
#include <IMP/exception.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <sys/stat.h>

IMPBFF_BEGIN_NAMESPACE

const double DEFAULT_ALLOWED_SPHERE_RADIUS = 2.1;

namespace avb {

const double DEFAULT_VDW = 1.70;

std::string upper(const std::string& s) {
    std::string out(s);
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(out[i])));
    }
    return out;
}

std::string trim(const std::string& s) {
    const std::string space = " \t\n\r\f\v";
    const std::size_t a = s.find_first_not_of(space);
    if (a == std::string::npos) return std::string();
    return s.substr(a, s.find_last_not_of(space) - a + 1);
}

const std::map<std::string, int>& element_numbers() {
    static std::map<std::string, int> table;
    if (table.empty()) {
        table["H"] = 1;   table["HE"] = 2;  table["LI"] = 3;  table["BE"] = 4;
        table["B"] = 5;   table["C"] = 6;   table["N"] = 7;   table["O"] = 8;
        table["F"] = 9;   table["MG"] = 12; table["SI"] = 14; table["P"] = 15;
        table["S"] = 16;  table["CL"] = 17; table["K"] = 19;  table["CA"] = 20;
        table["FE"] = 26; table["ZN"] = 30;
    }
    return table;
}

struct CacheKey {
    std::string path;
    long long mtime_ns, size;
    bool operator<(const CacheKey& o) const {
        if (path != o.path) return path < o.path;
        if (mtime_ns != o.mtime_ns) return mtime_ns < o.mtime_ns;
        return size < o.size;
    }
};

bool stat_key(const std::string& path, CacheKey* key) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) return false;
    key->path = path;
#ifdef __APPLE__
    key->mtime_ns = (long long) info.st_mtimespec.tv_sec * 1000000000LL +
                    info.st_mtimespec.tv_nsec;
#else
    key->mtime_ns = (long long) info.st_mtim.tv_sec * 1000000000LL +
                    info.st_mtim.tv_nsec;
#endif
    key->size = (long long) info.st_size;
    return true;
}

std::map<CacheKey, std::vector<PDBAtomRecord> >& record_cache() {
    static std::map<CacheKey, std::vector<PDBAtomRecord> > cache;
    return cache;
}

std::vector<PDBAtomRecord> parse_pdb(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("cannot open " << path, IOException);

    const std::map<int, double> radii = vdw_radii();
    const std::map<std::string, int>& numbers = element_numbers();

    std::vector<PDBAtomRecord> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.compare(0, 6, "ATOM  ") != 0 &&
            line.compare(0, 6, "HETATM") != 0) {
            continue;
        }
        if (line.size() < 54) continue;
        PDBAtomRecord row;
        const std::string res = trim(line.substr(22, 4));
        char* end = 0;
        row.resseq = static_cast<int>(std::strtol(res.c_str(), &end, 10));
        if (end == res.c_str() || *end != '\0') continue;
        const std::string xs = trim(line.substr(30, 8));
        const std::string ys = trim(line.substr(38, 8));
        const std::string zs = trim(line.substr(46, 8));
        row.x = std::strtod(xs.c_str(), &end); if (*end != '\0') continue;
        row.y = std::strtod(ys.c_str(), &end); if (*end != '\0') continue;
        row.z = std::strtod(zs.c_str(), &end); if (*end != '\0') continue;
        row.chain = trim(line.substr(21, 1));
        row.atom_name = trim(line.substr(12, 4));
        row.res_name = trim(line.substr(17, 3));
        // The serial is what CONECT records refer to. A record whose serial
        // does not parse keeps 0 rather than being dropped: the coordinates are
        // what most callers want, and only the bond reader needs the serial.
        const std::string serial = trim(line.substr(6, 5));
        char* serial_end = 0;
        const long value = std::strtol(serial.c_str(), &serial_end, 10);
        if (serial_end != serial.c_str() && *serial_end == '\0') {
            row.serial = static_cast<int>(value);
        }

        row.element = element_symbol_from_pdb_line(line);
        std::map<std::string, int>::const_iterator n = numbers.find(row.element);
        const int atomic_number = n == numbers.end() ? 0 : n->second;
        std::map<int, double>::const_iterator r = radii.find(atomic_number);
        row.vdw_radius = r == radii.end() ? DEFAULT_VDW : r->second;
        rows.push_back(row);
    }
    if (rows.empty()) {
        IMP_THROW("No ATOM/HETATM coordinates found in '" << path << "'",
                  IOException);
    }
    return rows;
}

}  // namespace avb

std::map<int, double> vdw_radii() {
    std::map<int, double> out;
    out[1] = 1.20;  out[2] = 1.40;  out[3] = 1.82;  out[4] = 1.53;
    out[5] = 1.92;  out[6] = 1.70;  out[7] = 1.55;  out[8] = 1.52;
    out[9] = 1.47;  out[12] = 1.73; out[14] = 2.10; out[15] = 1.80;
    out[16] = 1.80; out[17] = 1.75; out[19] = 2.27; out[20] = 1.97;
    out[26] = 1.56; out[30] = 1.39;
    return out;
}

std::string element_symbol_from_pdb_line(const std::string& line) {
    if (line.size() >= 78) {
        const std::string symbol = avb::upper(avb::trim(line.substr(76, 2)));
        if (!symbol.empty()) return symbol;
    }
    if (line.size() < 16) return std::string();

    std::string name_field = avb::upper(line.substr(12, 4));
    while (name_field.size() < 4) name_field += ' ';

    const std::map<std::string, int>& numbers = avb::element_numbers();
    const std::string candidate = avb::trim(name_field.substr(0, 2));
    bool tail_has_digit = false;
    for (std::size_t i = 2; i < name_field.size(); ++i) {
        if (std::isdigit(static_cast<unsigned char>(name_field[i]))) {
            tail_has_digit = true;
            break;
        }
    }
    if (candidate.size() == 2 && numbers.count(candidate) && !tail_has_digit) {
        return candidate;
    }

    std::string letters;
    for (std::size_t i = 0; i < name_field.size(); ++i) {
        if (std::isalpha(static_cast<unsigned char>(name_field[i]))) {
            letters += name_field[i];
        }
    }
    const std::string one = letters.substr(0, 1);
    const std::string two = letters.substr(0, std::min<std::size_t>(2, letters.size()));
    // Left-padding a two-letter element (" ZN ") breaks the column rule, but
    // where the strict reading is not an element at all, take the pair.
    if (!numbers.count(one) && two.size() == 2 && numbers.count(two)) return two;
    return one;
}

std::vector<PDBAtomRecord> read_pdb_records(const std::string& pdb_path) {
    avb::CacheKey key;
    if (avb::stat_key(pdb_path, &key)) {
        std::map<avb::CacheKey, std::vector<PDBAtomRecord> >& cache =
                avb::record_cache();
        std::map<avb::CacheKey, std::vector<PDBAtomRecord> >::const_iterator hit =
                cache.find(key);
        if (hit != cache.end()) return hit->second;
        const std::vector<PDBAtomRecord> rows = avb::parse_pdb(pdb_path);
        cache[key] = rows;
        return rows;
    }
    return avb::parse_pdb(pdb_path);
}

void load_structure_with_vdw(const std::string& pdb_path, double** out_view,
                             int* n_out_view) {
    const std::vector<PDBAtomRecord> rows = read_pdb_records(pdb_path);
    double* out = internal::new_double_view(rows.size() * 4, out_view, n_out_view);
    if (out == NULL) return;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        out[i * 4 + 0] = rows[i].x;
        out[i * 4 + 1] = rows[i].y;
        out[i * 4 + 2] = rows[i].z;
        out[i * 4 + 3] = rows[i].vdw_radius;
    }
}

void find_attachment_point(const std::string& pdb_path, const std::string& chain,
                           int resseq, const std::string& atom_name,
                           double** out_view, int* n_out_view) {
    const std::vector<PDBAtomRecord> rows = read_pdb_records(pdb_path);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        if (rows[i].resseq != resseq || rows[i].atom_name != atom_name) continue;
        if (!chain.empty() && rows[i].chain != chain) continue;
        double* out = internal::new_double_view(3, out_view, n_out_view);
        if (out == NULL) return;
        out[0] = rows[i].x;
        out[1] = rows[i].y;
        out[2] = rows[i].z;
        return;
    }
    // A miss stays a miss: a positional guess yields a plausible and entirely
    // wrong volume.
    internal::new_double_view(0, out_view, n_out_view);
}

AccessibleVolume resample_av(IMP::Model* model, IMP::Particle* source_particle,
                             double linker_length, double linker_width,
                             double r1, double r2, double r3, double disc_step,
                             double allowed_sphere_radius,
                             double contact_volume_thickness,
                             double contact_volume_trapped_fraction,
                             int search_stencil) {
    IMP::Particle* av_particle = new IMP::Particle(model);
    AV::do_setup_particle(model, av_particle->get_index(),
                          source_particle->get_index(), linker_length,
                          IMP::algebra::Vector3D(r1, r2, r3), linker_width,
                          allowed_sphere_radius, contact_volume_thickness,
                          contact_volume_trapped_fraction, disc_step);
    AV av(model, av_particle->get_index());
    if (search_stencil) av.set_search_stencil(search_stencil);
    av.resample();

    PathMap* path_map = av.get_map();
    const IMP::em::DensityHeader* header = path_map->get_header();
    const int nx = header->get_nx(), ny = header->get_ny(), nz = header->get_nz();

    // IMP orders the flat tile values `i = x + nx*y + nx*ny*z` -- *x* fastest --
    // so reading them into an (nx, ny, nz) array in the coordinate axis order
    // means walking z slowest. A straight copy would transpose the volume, and
    // a mirrored volume keeps the right voxel count, bounding box and total
    // volume, so only a voxel-by-voxel comparison against the point cloud
    // catches it.
    const std::vector<float> values = path_map->get_tile_values(
            PM_TILE_ACCESSIBLE_DENSITY,
            std::pair<double, double>(
                    0.0,
                    path_map->get_path_map_header().get_max_path_length()));
    std::vector<double> density(static_cast<std::size_t>(nx) * ny * nz, 0.0);
    for (int iz = 0; iz < nz; ++iz) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                const std::size_t flat =
                        static_cast<std::size_t>(iz) * ny * nx + iy * nx + ix;
                const std::size_t out =
                        (static_cast<std::size_t>(ix) * ny + iy) * nz + iz;
                if (flat < values.size()) density[out] = values[flat];
            }
        }
    }

    // The point cloud comes from IMP directly rather than from the grid above.
    // That is deliberate: it makes the cloud and the density **independent**
    // readings of the same volume, so a test comparing them catches a
    // mis-indexed grid. Deriving the points from the density would make any
    // indexing error self-consistent, and invisible.
    const std::vector<IMP::algebra::Vector4D> xyz_density =
            path_map->get_xyz_density();
    std::vector<double> points;
    points.reserve(xyz_density.size() * 4);
    for (std::size_t i = 0; i < xyz_density.size(); ++i) {
        for (int c = 0; c < 4; ++c) points.push_back(xyz_density[i][c]);
    }

    std::vector<double> origin(3);
    origin[0] = header->get_xorigin();
    origin[1] = header->get_yorigin();
    origin[2] = header->get_zorigin();

    const IMP::algebra::Vector3D source = av.get_source_coordinates();
    std::vector<double> attachment(3);
    attachment[0] = source[0];
    attachment[1] = source[1];
    attachment[2] = source[2];

    return AccessibleVolume(points, density, origin, header->get_spacing(), "",
                            attachment);
}

AccessibleVolume compute_av(double* atoms_xyzr, int n_atoms, int n_cols,
                            const std::vector<double>& source_xyz,
                            double linker_length, double linker_width,
                            double r1, double r2, double r3,
                            double grid_resolution,
                            double allowed_sphere_radius, int search_stencil) {
    if (n_cols != 4) {
        IMP_THROW("obstacles must be (N, 4) of x, y, z, radius, not (" << n_atoms
                          << ", " << n_cols << ")",
                  ValueException);
    }
    if (source_xyz.size() != 3) {
        IMP_THROW("the attachment point is three coordinates, not "
                          << source_xyz.size(),
                  ValueException);
    }

    IMP_NEW(IMP::Model, model, ());
    IMP::atom::Hierarchy root = IMP::atom::Hierarchy::setup_particle(
            new IMP::Particle(model));

    // Every atom is an obstacle, carrying a radius. The attachment site is one
    // of them rather than a separate massless marker: the source has to be a
    // real XYZR particle in the same hierarchy, and `allowed_sphere_radius`
    // handles that it is also an obstacle to itself.
    IMP::Particle* source_particle = NULL;
    for (int i = 0; i < n_atoms; ++i) {
        const double* a = atoms_xyzr + static_cast<std::size_t>(i) * 4;
        IMP::Particle* p = new IMP::Particle(model);
        IMP::core::XYZR xyzr = IMP::core::XYZR::setup_particle(p);
        xyzr.set_coordinates(IMP::algebra::Vector3D(a[0], a[1], a[2]));
        xyzr.set_radius(a[3]);
        IMP::atom::Hierarchy h = IMP::atom::Hierarchy::setup_particle(p);
        std::ostringstream name;
        name << "atom_" << i;
        h->set_name(name.str());
        root.add_child(h);
        if (source_particle == NULL &&
            std::fabs(a[0] - source_xyz[0]) < 1e-8 &&
            std::fabs(a[1] - source_xyz[1]) < 1e-8 &&
            std::fabs(a[2] - source_xyz[2]) < 1e-8) {
            source_particle = p;
        }
    }
    if (source_particle == NULL) {
        // The attachment site is not one of the atoms; add it as a zero-radius
        // particle so it is still a real XYZR to anchor to.
        source_particle = new IMP::Particle(model);
        IMP::core::XYZR xyzr = IMP::core::XYZR::setup_particle(source_particle);
        xyzr.set_coordinates(
                IMP::algebra::Vector3D(source_xyz[0], source_xyz[1], source_xyz[2]));
        xyzr.set_radius(0.0);
        IMP::atom::Hierarchy h =
                IMP::atom::Hierarchy::setup_particle(source_particle);
        h->set_name("source");
        root.add_child(h);
    }

    AccessibleVolume av = resample_av(model, source_particle, linker_length,
                                      linker_width, r1, r2, r3, grid_resolution,
                                      allowed_sphere_radius, 0.0, -1.0,
                                      search_stencil);
    // This door's convention: the attachment coordinate is the one the caller
    // handed in, not the one the solver snapped to a lattice.
    av.set_attachment_point(source_xyz);

    std::map<std::string, std::string> params;
    std::ostringstream s;
    params["backend"] = "imp_bff";
    s.str(""); s << linker_length;    params["linker_length"] = s.str();
    s.str(""); s << linker_width;     params["linker_width"] = s.str();
    s.str(""); s << r1 << ", " << r2 << ", " << r3;
    params["dye_radii"] = s.str();
    s.str(""); s << grid_resolution;  params["grid_resolution"] = s.str();
    s.str(""); s << allowed_sphere_radius;
    params["allowed_sphere_radius"] = s.str();
    s.str(""); s << (search_stencil ? search_stencil : 74);
    params["search_stencil"] = s.str();
    av.set_params(params);
    return av;
}

AccessibleVolume compute_av_from_structure(
        const std::string& pdb_path, const std::string& chain, int resseq,
        const std::string& atom_name, double linker_length, double linker_width,
        double r1, double r2, double r3, double disc_step,
        const std::string& strip_mask, double allowed_sphere_radius,
        double contact_volume_thickness,
        double contact_volume_trapped_fraction) {
    const std::string stripped =
            stripped_pdb_for(pdb_path, chain, resseq, atom_name, strip_mask);

    IMP_NEW(IMP::Model, model, ());
    IMP::atom::Hierarchy hierarchy = IMP::atom::read_pdb(
            stripped, model, new IMP::atom::NonWaterPDBSelector());

    IMP::atom::Selection sel(hierarchy);
    if (!chain.empty()) sel.set_chain_id(chain);
    sel.set_residue_index(resseq);
    sel.set_atom_type(IMP::atom::AtomType(atom_name));
    const IMP::ParticlesTemp particles = sel.get_selected_particles();
    if (particles.empty()) {
        IMP_THROW("Attachment site " << chain << ":" << resseq << ":"
                                     << atom_name << " not found",
                  ValueException);
    }

    // Source clearance. The path search inflates obstacles by half the linker
    // width, so the free sphere around the attachment atom has to clear that
    // inflation (plus a grid step of slack) or the source tile is walled in and
    // the volume comes back empty. The strip above already removes the
    // attachment residue's side chain, which is what lets FPS-calibrated small
    // clearances (`allowed_sphere_radius: 1`) compute a real cloud.
    const double clearance =
            allowed_sphere_radius >= 0.0
                    ? allowed_sphere_radius
                    : std::max(1.5, 0.5 * linker_width + 0.5 * disc_step);

    AccessibleVolume av = resample_av(
            model, particles[0], linker_length, linker_width, r1, r2, r3,
            disc_step, clearance, contact_volume_thickness,
            contact_volume_trapped_fraction, 0);

    // This door's convention: the point weights are forced to one. The array
    // door keeps whatever IMP reported.
    double* points = NULL;
    int n_points = 0;
    av.get_points(&points, &n_points);
    std::vector<double> uniform(points, points + n_points);
    for (int i = 3; i < n_points; i += 4) uniform[i] = 1.0;
    std::free(points);
    av.set_points(uniform);
    return av;
}

namespace {
//! One fps.json `Positions` entry, as the typed call.
/*! Reads the fields the fps dictionary states for a position. A declared
    `simulation_grid_resolution` that disagrees with `disc_step` raises: that
    field is *written into* the particle from `disc_step`, so a caller who
    declares it and omits the step would silently build at the 1.5 A default. */
AccessibleVolume av_from_position(const std::string& pdb_path,
                                  const nlohmann::json& position,
                                  double disc_step) {
    const std::string chain = position.value("chain_identifier", "");
    const int resseq = position.value("residue_seq_number", 0);
    const std::string atom = position.value("atom_name", "CA");
    const double linker_length = position.value("linker_length", 20.0);
    const double linker_width = position.value("linker_width", 1.0);
    const double r1 = position.value("radius1", 3.5);
    const double r2 = position.value("radius2", 0.0);
    const double r3 = position.value("radius3", 0.0);
    const double declared =
            position.contains("simulation_grid_resolution")
                    ? position.at("simulation_grid_resolution").get<double>()
                    : 0.0;
    double step = 1.5;
    if (disc_step > 0.0) {
        step = disc_step;
        if (position.contains("simulation_grid_resolution") &&
            std::abs(declared - step) > 1e-9) {
            IMP_THROW("simulation_grid_resolution=" << declared
                              << " in the position disagrees with disc_step="
                              << step,
                      ValueException);
        }
    } else if (position.contains("simulation_grid_resolution")) {
        IMP_THROW("the position declares simulation_grid_resolution="
                          << declared
                          << " but no disc_step was given, so the AV would be "
                             "built at disc_step=1.5",
                  ValueException);
    }
    return compute_av_from_structure(
            pdb_path, chain, resseq, atom, linker_length, linker_width, r1, r2,
            r3, step, position.value("strip_mask", ""),
            position.contains("allowed_sphere_radius")
                    ? position.at("allowed_sphere_radius").get<double>()
                    : -1.0,
            position.value("contact_volume_thickness", 0.0),
            position.value("contact_volume_trapped_fraction", -1.0));
}
}

AccessibleVolume compute_av_from_structure(
        const std::string& pdb_path, const std::string& position_json,
        double disc_step) {
    return av_from_position(pdb_path, nlohmann::json::parse(position_json),
                            disc_step);
}

std::map<std::string, AccessibleVolume> compute_avs_for_structure(
        const std::string& positions_json, const std::string& pdb_path_or_json,
        double disc_step) {
    const nlohmann::json positions = nlohmann::json::parse(positions_json);
    std::vector<std::string> paths;
    if (!pdb_path_or_json.empty() && pdb_path_or_json[0] == '[') {
        const nlohmann::json arr = nlohmann::json::parse(pdb_path_or_json);
        for (const auto& p : arr) paths.push_back(p.get<std::string>());
    } else {
        // One path, or a comma-separated list (the historical spelling).
        std::string cur;
        for (std::size_t i = 0; i <= pdb_path_or_json.size(); ++i) {
            const char c = i < pdb_path_or_json.size() ? pdb_path_or_json[i] : ',';
            if (c == ',') {
                if (!cur.empty()) paths.push_back(cur);
                cur.clear();
            } else {
                cur += c;
            }
        }
    }

    std::map<std::string, AccessibleVolume> out;
    if (!positions.is_object()) {
        IMP_THROW("positions must be a JSON object keyed by name",
                  ValueException);
    }
    for (auto it = positions.begin(); it != positions.end(); ++it) {
        const std::string& name = it.key();
        const nlohmann::json& p = it.value();
        const int body = p.value("body_id", 0);
        const std::string path =
                paths.empty() ? "" : paths[static_cast<std::size_t>(body) <
                                                   paths.size()
                                           ? body
                                           : 0];
        const double step = disc_step > 0.0
                ? disc_step
                : p.contains("simulation_grid_resolution")
                          ? p.at("simulation_grid_resolution").get<double>()
                          : 1.5;
        double* found = NULL;
        int n_found = 0;
        find_attachment_point(path, p.value("chain_identifier", ""),
                              static_cast<int>(p.value("residue_seq_number", 0)),
                              p.value("atom_name", "CA"), &found, &n_found);
        if (found == NULL || n_found == 0) {
            out[name] = AccessibleVolume(std::vector<double>(),
                                         std::vector<double>(),
                                         std::vector<double>(), step, name);
        } else {
            std::free(found);
            AccessibleVolume av = av_from_position(path, p, step);
            av.set_position_name(name);
            std::map<std::string, std::string> params;
            for (auto it2 = p.begin(); it2 != p.end(); ++it2) {
                params[it2.key()] = it2.value().dump();
            }
            av.set_params(params);
            out[name] = av;
        }
    }
    return out;
}

IMPBFF_END_NAMESPACE
