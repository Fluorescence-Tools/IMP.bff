/**
 * \file TrajectoryIO.cpp
 * \brief Reading rotamer-library trajectories from BinaryCIF.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/TrajectoryIO.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/Text.h>
#include <IMP/bff/IMPCompatibility.h>

// The C implementation of ihm, vendored by IMP at
// modules/core/dependency/python-ihm/src/. Only the header is included: the
// symbols come from libimp_atom, which bff already links, and which exports
// them because IMP::atom::read_mmcif compiles the same reader in. Compiling
// ihm_format.c here as well would work and would duplicate it in the library.
#include "ihm_format.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <fcntl.h>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! What the per-row callback accumulates into.
struct Collector {
    std::vector<double> x, y, z;
    ihm_keyword *kx = nullptr, *ky = nullptr, *kz = nullptr;
    bool count_only = false;
    long rows = 0;
};

void on_row(ihm_reader *, int, void *data, ihm_error **) {
    Collector *c = static_cast<Collector *>(data);
    ++c->rows;
    if (c->count_only) return;
    // An omitted coordinate is a malformed trajectory rather than a missing
    // optional, so it is read as zero and left to the caller's own checks --
    // throwing from inside a C callback would cross the C boundary.
    c->x.push_back(c->kx->omitted ? 0.0 : c->kx->data.fval);
    c->y.push_back(c->ky->omitted ? 0.0 : c->ky->data.fval);
    c->z.push_back(c->kz->omitted ? 0.0 : c->kz->data.fval);
}

//! Parse \p path, filling \p c. Throws on any reader error.
void parse(const std::string& path, const std::string& category, Collector& c) {
    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        IMP_THROW("cannot open BinaryCIF trajectory " << path, IMP::IOException);
    }
    ihm_file *fh = ihm_file_new_from_fd(fd);
    ihm_reader *reader = ihm_reader_new(fh, true);        // binary
    ihm_category *cat = ihm_category_new(reader, category.c_str(), on_row,
                                         nullptr, nullptr, &c, nullptr);
    c.kx = ihm_keyword_float_new(cat, "x");
    c.ky = ihm_keyword_float_new(cat, "y");
    c.kz = ihm_keyword_float_new(cat, "z");

    bool more = false;
    ihm_error *err = nullptr;
    const bool ok = ihm_read_file(reader, &more, &err);
    std::string message;
    if (!ok) message = err && err->msg ? err->msg : "unknown error";
    ihm_reader_free(reader);
    if (!ok) {
        IMP_THROW("reading " << path << ": " << message, IMP::IOException);
    }
}

}  // namespace

int bcif_trajectory_rows(const std::string& path, const std::string& category) {
    Collector c;
    c.count_only = true;
    parse(path, category, c);
    return static_cast<int>(c.rows);
}

void read_bcif_trajectory(const std::string& path, int n_atoms,
                          const std::string& category, double** out_view,
                          int* n_out_view) {
    if (n_atoms < 1) {
        IMP_THROW("n_atoms must be positive, not " << n_atoms,
                  IMP::ValueException);
    }
    Collector c;
    parse(path, category, c);
    const std::size_t rows = c.x.size();
    if (rows % static_cast<std::size_t>(n_atoms) != 0) {
        IMP_THROW("BinaryCIF trajectory " << path << " has " << rows
                  << " rows, which is not a multiple of n_atoms = " << n_atoms,
                  IMP::ValueException);
    }
    const std::size_t n_frames = rows / static_cast<std::size_t>(n_atoms);
    double* out = internal::new_double_view(rows * 3, out_view, n_out_view);
    if (out == nullptr) return;
    // Stored atom-major, returned frame-major: every consumer wants
    // (frame, atom, 3), and the storage order exists to make the deltas small
    // rather than to match anyone's indexing.
    for (std::size_t a = 0; a < static_cast<std::size_t>(n_atoms); ++a) {
        for (std::size_t f = 0; f < n_frames; ++f) {
            const std::size_t src = a * n_frames + f;
            const std::size_t dst = (f * n_atoms + a) * 3;
            out[dst + 0] = c.x[src];
            out[dst + 1] = c.y[src];
            out[dst + 2] = c.z[src];
        }
    }
}

// Named, not anonymous: IMP compiles this module as one translation unit.
namespace dcd {
using IMP::bff::internal::ends_with;


//! The leading Fortran record of a DCD header, in bytes.
const int HEADER_RECORD = 84;

//! Read one little- or big-endian 32-bit integer.
int read_i32(const unsigned char* p, bool big_endian) {
    int v = 0;
    if (big_endian) {
        v = (static_cast<int>(p[0]) << 24) | (static_cast<int>(p[1]) << 16) |
            (static_cast<int>(p[2]) << 8) | static_cast<int>(p[3]);
    } else {
        v = (static_cast<int>(p[3]) << 24) | (static_cast<int>(p[2]) << 16) |
            (static_cast<int>(p[1]) << 8) | static_cast<int>(p[0]);
    }
    return v;
}

float read_f32(const unsigned char* p, bool big_endian) {
    const int bits = read_i32(p, big_endian);
    float out;
    std::memcpy(&out, &bits, sizeof(out));
    return out;
}

std::vector<unsigned char> read_bytes(const std::string& path, std::size_t cap) {
    std::ifstream in(path.c_str(), std::ios::binary);
    if (!in) IMP_THROW("Cannot read " << path, IOException);
    std::vector<unsigned char> raw;
    raw.resize(cap);
    in.read(reinterpret_cast<char*>(&raw[0]), static_cast<std::streamsize>(cap));
    raw.resize(static_cast<std::size_t>(in.gcount()));
    return raw;
}

std::size_t file_size(const std::string& path) {
    std::ifstream in(path.c_str(), std::ios::binary | std::ios::ate);
    if (!in) IMP_THROW("Cannot read " << path, IOException);
    return static_cast<std::size_t>(in.tellg());
}

DCDHeader parse_header(const std::vector<unsigned char>& raw,
                       const std::string& path) {
    if (raw.size() < static_cast<std::size_t>(4 + HEADER_RECORD + 8)) {
        IMP_THROW(path << " is too short to be a DCD file", ValueException);
    }
    bool big_endian = false;
    if (read_i32(&raw[0], false) == HEADER_RECORD) {
        big_endian = false;
    } else if (read_i32(&raw[0], true) == HEADER_RECORD) {
        big_endian = true;
    } else {
        IMP_THROW("leading record is not 84 bytes in either byte order",
                  ValueException);
    }
    if (std::memcmp(&raw[4], "CORD", 4) != 0) {
        IMP_THROW(path << " does not start with the CORD magic", ValueException);
    }

    DCDHeader head;
    head.endianness = big_endian ? ">" : "<";
    head.n_frames = read_i32(&raw[8], big_endian);
    head.first_step = read_i32(&raw[12], big_endian);
    head.step_stride = read_i32(&raw[16], big_endian);
    // The float at offset 44 is the timestep; CHARMM writes it single precision.
    head.time_step = read_f32(&raw[44], big_endian);
    head.has_unit_cell = read_i32(&raw[48], big_endian) != 0;
    const bool four_dimensions = read_i32(&raw[48 + 12], big_endian) != 0;
    head.charmm_version = read_i32(&raw[84], big_endian);
    if (four_dimensions) {
        IMP_THROW(path << " is a 4-dimensional trajectory, which is not "
                          "supported",
                  ValueException);
    }

    std::size_t pos = 4 + HEADER_RECORD + 4;  // past the header record's marker
    // Title block: one Fortran record holding a count and that many 80-char
    // lines. Skipped whole; nothing here reads a title.
    const int title_len = read_i32(&raw[pos], big_endian);
    pos += 4 + static_cast<std::size_t>(title_len) + 4;

    if (pos + 8 > raw.size()) {
        IMP_THROW(path << " ends inside its atom-count record", ValueException);
    }
    const int natom_len = read_i32(&raw[pos], big_endian);
    if (natom_len != 4) {
        IMP_THROW(path << " has an unexpected atom-count record (" << natom_len
                       << " bytes)",
                  ValueException);
    }
    head.n_atoms = read_i32(&raw[pos + 4], big_endian);
    pos += 4 + static_cast<std::size_t>(natom_len) + 4;
    head.offset = static_cast<int>(pos);
    return head;
}


}  // namespace dcd

DCDHeader read_dcd_header(const std::string& path) {
    // The header, title block and atom count sit in the first few hundred
    // bytes; reading the whole trajectory to parse them made a header-only call
    // as expensive as a full load. 64 KiB is more than any title block needs.
    return dcd::parse_header(dcd::read_bytes(path, 65536), path);
}

void read_dcd(const std::string& path, int max_frames, double** out_view,
              int* n_out_view) {
    const std::vector<unsigned char> raw =
            dcd::read_bytes(path, dcd::file_size(path));
    const DCDHeader head = dcd::parse_header(raw, path);
    const bool big_endian = head.endianness == ">";

    int n_frames = head.n_frames;
    if (max_frames >= 0 && max_frames < n_frames) n_frames = max_frames;
    const std::size_t axis_bytes = 4 * static_cast<std::size_t>(head.n_atoms);

    double* out = internal::new_double_view(
            static_cast<std::size_t>(n_frames) * head.n_atoms * 3, out_view,
            n_out_view);
    if (out == nullptr) return;

    std::size_t pos = static_cast<std::size_t>(head.offset);
    for (int frame = 0; frame < n_frames; ++frame) {
        if (head.has_unit_cell) {
            // A 48-byte double record ahead of the coordinates.
            const int cell_len = dcd::read_i32(&raw[pos], big_endian);
            pos += 4 + static_cast<std::size_t>(cell_len) + 4;
        }
        for (int axis = 0; axis < 3; ++axis) {
            if (pos + 4 > raw.size()) {
                std::free(out);
                *out_view = nullptr;
                *n_out_view = 0;
                IMP_THROW(path << ": frame " << frame << " ends before its axis "
                               << axis << " record",
                          ValueException);
            }
            const std::size_t rec_len = static_cast<std::size_t>(
                    dcd::read_i32(&raw[pos], big_endian));
            if (rec_len != axis_bytes) {
                std::free(out);
                *out_view = nullptr;
                *n_out_view = 0;
                IMP_THROW(path << ": frame " << frame << " axis " << axis
                               << " record is " << rec_len << " bytes, expected "
                               << axis_bytes << " for " << head.n_atoms
                               << " atoms",
                          ValueException);
            }
            const std::size_t start = pos + 4;
            for (int i = 0; i < head.n_atoms; ++i) {
                out[(static_cast<std::size_t>(frame) * head.n_atoms + i) * 3 +
                    axis] = dcd::read_f32(&raw[start + 4 * i], big_endian);
            }
            pos = start + rec_len + 4;
        }
    }
}

void read_trajectory(const std::string& path, int n_atoms, int max_frames,
                     double** out_view, int* n_out_view) {
    if (dcd::ends_with(path, ".bcif")) {
        if (n_atoms <= 0) {
            IMP_THROW("reading a BinaryCIF trajectory needs n_atoms: it stores "
                      "one row per (atom, frame) and cannot infer where frames "
                      "divide",
                      ValueException);
        }
        double* flat = nullptr;
        int n_flat = 0;
        read_bcif_trajectory(path, n_atoms, "_rotamer_coord", &flat, &n_flat);
        std::size_t keep = static_cast<std::size_t>(n_flat);
        if (max_frames >= 0) {
            const std::size_t capped =
                    static_cast<std::size_t>(max_frames) * n_atoms * 3;
            if (capped < keep) keep = capped;
        }
        double* out = internal::new_double_view(keep, out_view, n_out_view);
        if (out != nullptr) std::memcpy(out, flat, keep * sizeof(double));
        std::free(flat);
        return;
    }
    if (dcd::ends_with(path, ".dcd")) {
        read_dcd(path, max_frames, out_view, n_out_view);
        return;
    }
    IMP_THROW("unsupported trajectory format: " << path, ValueException);
}

IMPBFF_END_NAMESPACE
