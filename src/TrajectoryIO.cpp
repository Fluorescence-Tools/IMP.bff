/**
 * \file TrajectoryIO.cpp
 * \brief Reading rotamer-library trajectories from BinaryCIF.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/TrajectoryIO.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/exception.h>

// The C implementation of ihm, vendored by IMP at
// modules/core/dependency/python-ihm/src/. Only the header is included: the
// symbols come from libimp_atom, which bff already links, and which exports
// them because IMP::atom::read_mmcif compiles the same reader in. Compiling
// ihm_format.c here as well would work and would duplicate it in the library.
#include "ihm_format.h"

#include <cstdio>
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

IMPBFF_END_NAMESPACE
