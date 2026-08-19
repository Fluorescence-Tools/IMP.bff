/**
 *  \file IMP/bff/GridDiffusionSolver.h
 *  \brief The field picture of a tethered dye: an occupancy density, propagated.
 *
 * The counterpart of #DyeDiffusionSimulation. That one walks a particle and
 * resolves the dye's history; this one propagates the density of where the dye
 * is, which is cheaper and smoother but cannot produce a correlation function
 * because it has thrown the history away. They agree on equilibrium and on the
 * mean squared displacement.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_GRIDDIFFUSIONSOLVER_H
#define IMPBFF_GRIDDIFFUSIONSOLVER_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/DiffusionSolver.h>

#include <IMP/value_macros.h>
#include <IMP/showable_macros.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Propagate an excited-state density on a masked 3-D grid.
/*!
    \param density initial occupancy, flat `ng^3`
    \param bounds non-zero inside the accessible volume
    \param diffusion_map per-voxel `D`, A^2/ns
    \param rate_map per-voxel deactivation rate, 1/ns
    \param dg voxel edge, A
    \param t_step time step, ns
    \param flux_form FLUX_SMOLUCHOWSKI or FLUX_ITO
    \param check_stability refuse a step past the explicit scheme's limit
*/
class IMPBFFEXPORT GridDiffusionSolver {
    std::vector<double> density_;
    std::vector<double> bounds_;
    std::vector<double> diffusion_map_;
    std::vector<double> rate_map_;
    double dg_;
    double t_step_;
    int ng_;
    int flux_form_;
    bool check_stability_;
    int n_iterations_;

    //! \throws IMP::ValueException on a domain touching the outer shell, or a
    //!         time step past `dg^2 / (6 D_max)`. The rate does **not** enter
    //!         that bound: it is carried as `exp(-k dt)`, the exact solution
    //!         over the step, and so constrains nothing.
    void validate() const;
    //! `p b / sum(p b)` -- the normalised start both run() and gradient() use.
    std::vector<double> start_density(const std::vector<double>& p) const;

public:
    GridDiffusionSolver(
            const std::vector<double>& density = std::vector<double>(),
            const std::vector<double>& bounds = std::vector<double>(),
            const std::vector<double>& diffusion_map = std::vector<double>(),
            const std::vector<double>& rate_map = std::vector<double>(),
            double dg = 1.0, double t_step = 0.01,
            int flux_form = FLUX_SMOLUCHOWSKI, bool check_stability = true);

    //! Integrate \p n_steps steps, reporting every \p n_out.
    /*!
        \param[out] out_view,n_out_view the surviving fraction at each reported
                    step -- the donor decay -- followed by the final density.
                    One buffer rather than two returns, because a numpy view is
                    one array; the shim splits it.
    */
    void run(int n_steps, int n_out, double** out_view, int* n_out_view);

    //! Adjoint of run(): the gradient of a loss on the decay, every voxel at once.
    /*!
        For \f$L = \sum_k \bar F_k F(t_k)\f$, returns dL/dD, dL/dk and dL/dp0
        for about four times the cost of one run(). Finite differences need one
        run *per parameter*.

        The gradient is of the *discrete* scheme actually run, so it agrees with
        a central difference of run() to roundoff -- which is also why the time
        step must not depend on the parameters being differentiated.

        \param[in] dL_dF sensitivity of the loss to each reported population
        \param[in] density the start of the run being differentiated; run()
                   replaces the solver's own with the final one, so a completed
                   run has to be handed the density it began with
        \param[out] out_view,n_out_view three concatenated `ng^3` blocks:
                    dL/dD, dL/dk, dL/dp0
    */
    void gradient(const std::vector<double>& dL_dF, int n_steps, int n_out,
                  const std::vector<double>& density,
                  double** out_view, int* n_out_view);

    //! Propagate with no decay until the occupancy stops moving.
    /*! Prefer equilibrium_occupancy(), which is the same answer in closed form.
        This iterates toward it, and on a real site -- where `D` spans orders of
        magnitude through the compounding slow factor -- it can fail to converge
        in any practical number of steps. */
    void equilibrium(int n_steps, double tolerance, int n_check,
                     double** out_view, int* n_out_view);

    void get_density(double** out_view, int* n_out_view) const;
    void set_density(const std::vector<double>& d) { density_ = d; }
    void get_bounds(double** out_view, int* n_out_view) const;
    void get_diffusion_map(double** out_view, int* n_out_view) const;
    void set_diffusion_map(const std::vector<double>& d) { diffusion_map_ = d; }
    void get_rate_map(double** out_view, int* n_out_view) const;
    void set_rate_map(const std::vector<double>& k) { rate_map_ = k; }

    int get_ng() const { return ng_; }
    double get_dg() const { return dg_; }
    double get_t_step() const { return t_step_; }
    int get_n_iterations() const { return n_iterations_; }
    int get_flux_form() const { return flux_form_; }

    IMP_SHOWABLE_INLINE(GridDiffusionSolver,
                        out << "GridDiffusionSolver(ng = " << ng_
                            << ", dg = " << dg_ << " A, dt = " << t_step_
                            << " ns, " << n_iterations_ << " steps taken)");
};
IMP_VALUES(GridDiffusionSolver, GridDiffusionSolvers);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_GRIDDIFFUSIONSOLVER_H
