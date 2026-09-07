/*
 * The cgprobe MD runner is a program -- it reads a directory of force-field
 * systems, propagates each one and writes trajectories -- so it lives in
 * `bin/imp_bff`, beside the command that runs it.
 *
 * What it *computes* does not. The restraint builders
 * (`build_probe_restraints`, `build_steric_restraint`, `build_go_restraints`)
 * and the placement search (`place_guest_by_score`) are C++ in `Scoring.h`,
 * because those are kernels rather than orchestration -- and because a runner
 * with its own copy of the restraint builder disagrees with the shared one
 * about what to do when a bond carries no equilibrium length and about how to
 * score repulsion.
 *
 * There is one repulsion, for dynamics and for Monte Carlo alike:
 * `build_steric_restraint`, soft spheres over every non-excluded pair.
 */
