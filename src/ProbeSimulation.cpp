/**
 * \file ProbeSimulation.cpp
 * \brief The array views of a recorded run.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/ProbeSimulation.h>
#include <IMP/bff/internal/OutputView.h>

IMPBFF_BEGIN_NAMESPACE

void ProbeSimulationTrajectory::get_coordinates(double** out_view, int* n_out_view) const {
    internal::copy_to_view(coordinates, out_view, n_out_view);
}
void ProbeSimulationTrajectory::get_times_fs(double** out_view, int* n_out_view) const {
    internal::copy_to_view(times_fs, out_view, n_out_view);
}
void ProbeSimulationTrajectory::get_potential_energy(double** out_view, int* n_out_view) const {
    internal::copy_to_view(potential_energy, out_view, n_out_view);
}
void ProbeSimulationTrajectory::get_kinetic_energy(double** out_view, int* n_out_view) const {
    internal::copy_to_view(kinetic_energy, out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
