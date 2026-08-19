/**
 * \file DiffusionSolver.cpp
 * \brief Explicit propagation of an excited-state density on an AV grid.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/DiffusionSolver.h>
#include <IMP/bff/internal/OutputView.h>
// The solver and its adjoint are tttrlib's LatticeDiffusion.h, vendored
// verbatim (test/test_vendored_headers.py keeps the copy identical). This file
// is the IMP-facing shell: the enum names IMP.bff has always exported, the
// numpy out_view returns, nothing numerical.
#define TTTRLIB_LATTICE_NAMESPACE IMP::bff::internal
#include <IMP/bff/internal/LatticeDiffusion.h>

#include <algorithm>


IMPBFF_BEGIN_NAMESPACE

void diffusion_step(
        const std::vector<double>& cur, const std::vector<double>& d,
        const std::vector<double>& decay, const std::vector<double>& bounds,
        int ng, int flux_form, double** out_view, int* n_out_view) {
    std::vector<double> nxt(cur.size(), 0.0);
    internal::lattice_sweep(cur, d, decay, bounds, static_cast<std::size_t>(ng),
                            flux_form == FLUX_SMOLUCHOWSKI, nxt);
    internal::copy_to_view(nxt, out_view, n_out_view);
}

void diffusion_propagate(
        const std::vector<double>& cur, const std::vector<double>& d,
        const std::vector<double>& decay, const std::vector<double>& bounds,
        int ng, int flux_form, int n_steps, int n_out,
        std::vector<double>& fluorescence,
        double** out_view, int* n_out_view) {
    const std::vector<double> final_density = internal::lattice_propagate(
            cur, d, decay, bounds, ng,
            flux_form == FLUX_SMOLUCHOWSKI ? internal::LATTICE_FLUX_SMOLUCHOWSKI : internal::LATTICE_FLUX_ITO,
            n_steps, n_out, fluorescence);
    // One memcpy, against ng^3 x ~35-40 ns of marshalling.
    internal::copy_to_view(final_density, out_view, n_out_view);
}

void diffusion_propagate_adjoint(
        const std::vector<double>& cur, const std::vector<double>& d,
        const std::vector<double>& decay, const std::vector<double>& bounds,
        int ng, int flux_form, int n_steps, int n_out,
        const std::vector<double>& dL_dF,
        std::vector<double>& dL_dd, std::vector<double>& dL_ddecay,
        std::vector<double>& dL_dcur) {
    internal::lattice_propagate_adjoint(
            cur, d, decay, bounds, ng,
            flux_form == FLUX_SMOLUCHOWSKI ? internal::LATTICE_FLUX_SMOLUCHOWSKI : internal::LATTICE_FLUX_ITO,
            n_steps, n_out, dL_dF, dL_dd, dL_ddecay, dL_dcur);
}

void diffusion_propagate_adjoint(
        const std::vector<double>& cur, const std::vector<double>& d,
        const std::vector<double>& decay, const std::vector<double>& bounds,
        int ng, int flux_form, int n_steps, int n_out,
        const std::vector<double>& dL_dF,
        double** out_view, int* n_out_view) {
    std::vector<double> gd, gdec, gcur;
    diffusion_propagate_adjoint(cur, d, decay, bounds, ng, flux_form, n_steps, n_out,
                                dL_dF, gd, gdec, gcur);
    const std::size_t nv = gd.size();
    double* buf = internal::new_double_view(3 * nv, out_view, n_out_view);
    if (buf == nullptr) return;
    std::copy(gd.begin(), gd.end(), buf);
    std::copy(gdec.begin(), gdec.end(), buf + nv);
    std::copy(gcur.begin(), gcur.end(), buf + 2 * nv);
}

IMPBFF_END_NAMESPACE
