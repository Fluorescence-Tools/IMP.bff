/**
 * \file AVBuilder.cpp
 * \brief Building an accessible volume without an IMP::Model: the array door,
 *        the lattice search behind it, and the PDB read.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/internal/Text.h>
#include <IMP/bff/AVBuilder.h>

#include <IMP/bff/DensityGrid.h>
#include <IMP/bff/PathMap.h>
#include <IMP/bff/OccupancyGrid.h>
#include <IMP/bff/StripMask.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/json.h>

#include <IMP/algebra/Vector3D.h>
#include <IMP/bff/Base.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <sys/stat.h>

IMPBFF_BEGIN_NAMESPACE


const double DEFAULT_ALLOWED_SPHERE_RADIUS = 2.1;

namespace avb {

using internal::upper;

const double DEFAULT_VDW = 1.70;


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
        const std::string res = internal::trimmed(line.substr(22, 4));
        char* end = 0;
        row.resseq = static_cast<int>(std::strtol(res.c_str(), &end, 10));
        if (end == res.c_str() || *end != '\0') continue;
        const std::string xs = internal::trimmed(line.substr(30, 8));
        const std::string ys = internal::trimmed(line.substr(38, 8));
        const std::string zs = internal::trimmed(line.substr(46, 8));
        row.x = std::strtod(xs.c_str(), &end); if (*end != '\0') continue;
        row.y = std::strtod(ys.c_str(), &end); if (*end != '\0') continue;
        row.z = std::strtod(zs.c_str(), &end); if (*end != '\0') continue;
        row.chain = internal::trimmed(line.substr(21, 1));
        row.atom_name = internal::trimmed(line.substr(12, 4));
        row.res_name = internal::trimmed(line.substr(17, 3));
        // The serial is what CONECT records refer to. A record whose serial
        // does not parse keeps 0 rather than being dropped: the coordinates are
        // what most callers want, and only the bond reader needs the serial.
        const std::string serial = internal::trimmed(line.substr(6, 5));
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

double vdw_radius(const std::string& element) {
    const std::map<std::string, int>& numbers = avb::element_numbers();
    const std::string symbol = avb::upper(internal::trimmed(element));
    std::map<std::string, int>::const_iterator n = numbers.find(symbol);
    if (n == numbers.end()) return 1.7;
    const std::map<int, double> radii = vdw_radii();
    std::map<int, double>::const_iterator r = radii.find(n->second);
    return r == radii.end() ? 1.7 : r->second;
}

std::string element_symbol_from_pdb_line(const std::string& line) {
    if (line.size() >= 78) {
        const std::string symbol = avb::upper(internal::trimmed(line.substr(76, 2)));
        if (!symbol.empty()) return symbol;
    }
    if (line.size() < 16) return std::string();

    std::string name_field = avb::upper(line.substr(12, 4));
    while (name_field.size() < 4) name_field += ' ';

    const std::map<std::string, int>& numbers = avb::element_numbers();
    const std::string candidate = internal::trimmed(name_field.substr(0, 2));
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

void get_attachment_point(const std::string& pdb_path, const std::string& chain,
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


/* The accessible-volume search over spheres: the decorator's lattice path
   (AV::resample_lattice_*), straight through, with none of the per-handle
   caching -- the window, the two occupancy rasters, the attachment-atom
   subtraction, the bounded search, the dye-radius carve, the contact
   weighting, and the read-out. Every step calls the same PathMap and
   OccupancyGrid entry points the decorator calls, in the same order, so the
   two doors compute the same volume; the array-door records in
   test/test_density_grid.py pin that. */
namespace {

//! The lattice search itself: the map every read-out below is taken from.
/*! Split out of get_av_lattice() so that the accessible volume and the
    linker path lengths (#get_linker_path_lengths, LabelLib's
    `minLinkerLength`) are two read-outs of one search rather than two
    searches. */
IMP::Pointer<PathMap> run_lattice_search(
        const std::vector<IMP::algebra::Vector4D>& spheres,
        const IMP::algebra::Vector3D& source, double source_radius,
        double linker_length, double linker_width,
        double r1, double r2, double r3, double h,
        double allowed_sphere_radius, int search_stencil,
        double contact_volume_thickness,
        double contact_volume_trapped_fraction) {
    // the effective parameters, as AV::get_effective_* derive them
    const double ll = linker_length;   // stencil compensation is opt-in and off
    const double allowed = allowed_sphere_radius >= 0.0
            ? allowed_sphere_radius
            : std::max(1.5, 0.5 * linker_width + 0.5 * h);
    std::vector<double> dye_radii;
    dye_radii.push_back(r1);
    if (r2 > 0.0) dye_radii.push_back(r2);
    if (r3 > 0.0) dye_radii.push_back(r3);

    // 0. the window, and the header anchored on it (AV::create_path_map_header)
    int k0[3]; int n;
    lattice_window(source, ll, h, k0, n);
    PathMapHeader header(ll, h);
    header.update_map_dimensions(n, n, n);
    if (search_stencil == 26) {
        header.set_neighbor_radius(std::sqrt(3.0) + 1e-6);
    } else if (search_stencil == 74) {
        header.set_neighbor_radius(std::sqrt(6.0) + 1e-6);
    }
    const IMP::algebra::Vector3D grid_origin(k0[0] * h, k0[1] * h, k0[2] * h);
    header.set_path_origin(source, grid_origin);
    IMP_NEW(PathMap, map, (header, "PathMap%1%"));
    if (search_stencil == 26 || search_stencil == 74) map->set_symmetric_stencil(true);
    map->set_euclidean_search(false);
    map->set_origin_fast(grid_origin);

    // 1. the two occupancy rasters: half the linker width, and the dye radius
    std::shared_ptr<const std::vector<IMP::algebra::Vector4D> > snap =
            std::make_shared<const std::vector<IMP::algebra::Vector4D> >(spheres);
    const double extra1 = linker_width * 0.5;
    IMP_NEW(OccupancyGrid, occ1, (h, extra1, snap));
    occ1->set_was_used(true);
    occ1->set_window(k0[0], k0[1], k0[2], n, n, n);
    occ1->update(true);
    std::vector<IMP::Pointer<OccupancyGrid> > occ_dye;
    for (std::size_t i = 0; i < dye_radii.size(); ++i) {
        IMP_NEW(OccupancyGrid, o, (h, dye_radii[i], snap));
        o->set_was_used(true);
        o->set_window(k0[0], k0[1], k0[2], n, n, n);
        o->update(true);
        occ_dye.push_back(o);
    }

    // 3. obstacles inflated by half the linker width, the attachment atom dropped
    const long nvox = map->get_number_of_voxels();
    std::vector<int32_t> counts(nvox);
    occ1->read_window_counts(k0[0], k0[1], k0[2], n, n, n, counts.data());
    drop_source_obstruction(counts.data(), grid_origin, n, h, source,
                            source_radius, extra1);

    // 4./5./6. block beyond the linker length, open the allowed sphere, search
    const long source_idx = map->get_voxel_by_location(source);
    const float hf = static_cast<float>(h);
    const float llf = static_cast<float>(map->get_path_map_header().get_max_path_length());
    float bound = llf / hf;
    while (bound * hf < llf) bound = std::nextafter(bound, std::numeric_limits<float>::infinity());
    map->search_lattice(source_idx, bound, source, ll, allowed, counts.data());

    // 7. remove tiles closer to obstacles than the dye radius
    occ_dye[0]->read_window_counts(k0[0], k0[1], k0[2], n, n, n, counts.data());
    drop_source_obstruction(counts.data(), grid_origin, n, h, source,
                            source_radius, dye_radii[0]);
    if (dye_radii.size() <= 1) {
        map->carve_lattice(counts.data());
    } else {
        std::vector<std::vector<int32_t> > more(dye_radii.size() - 1,
                                                std::vector<int32_t>(nvox));
        std::vector<const int32_t*> src(dye_radii.size());
        src[0] = counts.data();
        for (std::size_t i = 1; i < dye_radii.size(); ++i) {
            occ_dye[i]->read_window_counts(k0[0], k0[1], k0[2], n, n, n, more[i - 1].data());
            drop_source_obstruction(more[i - 1].data(), grid_origin, n, h, source,
                                    source_radius, dye_radii[i]);
            src[i] = more[i - 1].data();
        }
        map->carve_lattice_fractional(src.data(), static_cast<int>(dye_radii.size()));
    }
    map->apply_contact_weighting(counts.data(), contact_volume_thickness,
                                 contact_volume_trapped_fraction);

    return map;
}

//! The caller's spheres, and where the linker is tied, for the array doors.
/*! The attachment site takes part as one of the obstacles when it coincides
    with an atom -- then its own radius is subtracted from the raster around
    it, as the decorator does for the attachment atom; otherwise it is a
    point with no size, which blocks nothing. This is exactly what the
    IMP::Model these doors used to build expressed with particles. */
IMP::algebra::Vector3D prepare_obstacles(
        double* atoms_xyzr, int n_atoms, const std::vector<double>& source_xyz,
        std::vector<IMP::algebra::Vector4D>* spheres, double* source_radius) {
    spheres->reserve(static_cast<std::size_t>(n_atoms) + 1);
    *source_radius = 0.0;
    bool source_is_atom = false;
    for (int i = 0; i < n_atoms; ++i) {
        const double* a = atoms_xyzr + static_cast<std::size_t>(i) * 4;
        spheres->push_back(IMP::algebra::Vector4D(a[0], a[1], a[2], a[3]));
        if (!source_is_atom &&
            std::fabs(a[0] - source_xyz[0]) < 1e-8 &&
            std::fabs(a[1] - source_xyz[1]) < 1e-8 &&
            std::fabs(a[2] - source_xyz[2]) < 1e-8) {
            source_is_atom = true;
            *source_radius = a[3];
        }
    }
    const IMP::algebra::Vector3D source(source_xyz[0], source_xyz[1], source_xyz[2]);
    if (!source_is_atom) {
        spheres->push_back(IMP::algebra::Vector4D(source[0], source[1], source[2], 0.0));
    }
    return source;
}

}  // namespace

AccessibleVolume get_av_lattice(
        const std::vector<IMP::algebra::Vector4D>& spheres,
        const IMP::algebra::Vector3D& source, double source_radius,
        double linker_length, double linker_width,
        double r1, double r2, double r3, double h,
        double allowed_sphere_radius, int search_stencil,
        double contact_volume_thickness,
        double contact_volume_trapped_fraction) {
    IMP::Pointer<PathMap> map = run_lattice_search(spheres, source, source_radius, linker_length, linker_width,
                              r1, r2, r3, h, allowed_sphere_radius, search_stencil,
                              contact_volume_thickness, contact_volume_trapped_fraction);
    // the read-out: the density cube in (x, y, z) order, the cloud, the frame
    const GridHeader* gh = map->get_header();
    const int nx = gh->get_nx(), ny = gh->get_ny(), nz = gh->get_nz();
    const std::vector<float> values = map->get_tile_values(
            PM_TILE_ACCESSIBLE_DENSITY,
            std::pair<double, double>(0.0, map->get_path_map_header().get_max_path_length()));
    std::vector<double> density(static_cast<std::size_t>(nx) * ny * nz, 0.0);
    for (int iz = 0; iz < nz; ++iz) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                const std::size_t flat = static_cast<std::size_t>(iz) * ny * nx + iy * nx + ix;
                const std::size_t out = (static_cast<std::size_t>(ix) * ny + iy) * nz + iz;
                if (flat < values.size()) density[out] = values[flat];
            }
        }
    }
    const std::vector<IMP::algebra::Vector4D> xyz_density = map->get_xyz_density();
    std::vector<double> points;
    points.reserve(xyz_density.size() * 4);
    for (std::size_t i = 0; i < xyz_density.size(); ++i) {
        for (int c = 0; c < 4; ++c) points.push_back(xyz_density[i][c]);
    }
    std::vector<double> origin(3);
    origin[0] = gh->get_xorigin();
    origin[1] = gh->get_yorigin();
    origin[2] = gh->get_zorigin();
    std::vector<double> attachment(3);
    attachment[0] = source[0];
    attachment[1] = source[1];
    attachment[2] = source[2];
    return AccessibleVolume(points, density, origin, gh->get_spacing(), "", attachment);
}

DensityGrid* get_linker_path_lengths(
        double* atoms_xyzr, int n_atoms, int n_cols,
        const std::vector<double>& source_xyz,
        double linker_length, double linker_width, double dye_radius,
        double grid_resolution, double allowed_sphere_radius, int search_stencil) {
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
    std::vector<IMP::algebra::Vector4D> spheres;
    double source_radius = 0.0;
    const IMP::algebra::Vector3D source =
            prepare_obstacles(atoms_xyzr, n_atoms, source_xyz, &spheres, &source_radius);
    return get_linker_path_lengths(spheres, source, source_radius, linker_length,
                                   linker_width, dye_radius, grid_resolution,
                                   allowed_sphere_radius,
                                   search_stencil ? search_stencil : 74);
}

DensityGrid* get_linker_path_lengths(
        const std::vector<IMP::algebra::Vector4D>& spheres,
        const IMP::algebra::Vector3D& source, double source_radius,
        double linker_length, double linker_width, double dye_radius, double h,
        double allowed_sphere_radius, int search_stencil) {
    IMP::Pointer<PathMap> map = run_lattice_search(
            spheres, source, source_radius, linker_length, linker_width,
            dye_radius, 0.0, 0.0, h, allowed_sphere_radius, search_stencil, 0.0, -1.0);
    const GridHeader* gh = map->get_header();
    const int nx = gh->get_nx(), ny = gh->get_ny(), nz = gh->get_nz();
    const std::vector<float> values = map->get_tile_values(
            PM_TILE_PATH_LENGTH,
            std::pair<double, double>(0.0, map->get_path_map_header().get_max_path_length()));
    IMP_NEW(DensityGrid, out, ("LinkerPathLengths%1%"));
    GridHeader* oh = out->get_header_writable();
    oh->set_spacing(gh->get_spacing());
    oh->update_map_dimensions(nx, ny, nz);
    out->resize(static_cast<long>(nx) * ny * nz);
    out->set_origin(IMP::algebra::Vector3D(gh->get_xorigin(), gh->get_yorigin(),
                                           gh->get_zorigin()));
    // (x, y, z) order, as the volume's density cube is written. A voxel the
    // linker cannot reach is negative, which is LabelLib's convention for the
    // same grid (minLinkerLength); the attachment voxel itself is a reached
    // voxel of length zero.
    const long source_idx = map->get_voxel_by_location(source);
    for (int iz = 0; iz < nz; ++iz) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                const std::size_t flat = static_cast<std::size_t>(iz) * ny * nx + iy * nx + ix;
                const std::size_t o = (static_cast<std::size_t>(ix) * ny + iy) * nz + iz;
                double v = -1.0;
                if (flat < values.size()) {
                    const double raw = values[flat];
                    v = raw > 0.0 ? raw
                        : (static_cast<long>(flat) == source_idx ? 0.0 : -1.0);
                }
                out->set_value(static_cast<long>(o), v);
            }
        }
    }
    out->set_was_used(true);
    return out.release();
}


AccessibleVolume get_av(double* atoms_xyzr, int n_atoms, int n_cols,
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

    std::vector<IMP::algebra::Vector4D> spheres;
    double source_radius = 0.0;
    const IMP::algebra::Vector3D source =
            prepare_obstacles(atoms_xyzr, n_atoms, source_xyz, &spheres, &source_radius);

    AccessibleVolume av = get_av_lattice(
            spheres, source, source_radius, linker_length, linker_width,
            r1, r2, r3, grid_resolution, allowed_sphere_radius,
            search_stencil ? search_stencil : 74, 0.0, -1.0);
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


AccessibleVolume get_av_from_pdb(
        const std::string& pdb_path, const std::string& chain, int resseq,
        const std::string& atom_name, double linker_length, double linker_width,
        double r1, double r2, double r3, double grid_resolution,
        double allowed_sphere_radius, int search_stencil) {
    double* atoms = 0;
    int n_atoms = 0;
    load_structure_with_vdw(pdb_path, &atoms, &n_atoms);
    double* src = 0;
    int n_src = 0;
    get_attachment_point(pdb_path, chain, resseq, atom_name, &src, &n_src);
    if (n_src != 3) {
        std::free(atoms);
        std::free(src);
        IMP_THROW("no attachment atom " << chain << ":" << resseq << ":" << atom_name
                  << " in " << pdb_path, ValueException);
    }
    const std::vector<double> source(src, src + 3);
    AccessibleVolume av = get_av(atoms, n_atoms / 4, 4, source, linker_length, linker_width,
                                 r1, r2, r3, grid_resolution, allowed_sphere_radius,
                                 search_stencil);
    std::free(atoms);
    std::free(src);
    return av;
}

IMPBFF_END_NAMESPACE
