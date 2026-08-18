"""Moved to :mod:`IMP.bff.representation.distance` by PRD-113 stage 4b.

A distance between two labels is a property of what represents them -- an
accessible volume, a rotamer library, a Gaussian -- not of the FRET engine.
Re-exported here so call sites inside this package keep resolving while the rest
of ``fret/`` migrates.
"""

from IMP.bff.representation.distance import *  # noqa: F401,F403
from IMP.bff.representation.distance import (  # noqa: F401
    N_DISTANCE_SAMPLES,
    av_pair_statistics,
    average_distance,
    chi2_score,
    distance_between_mean_positions,
    distance_from_fret_efficiency,
    fit_transfer_polynomial,
    fret_efficiency,
    fret_pair_distribution,
    fret_pair_efficiencies,
    fret_pair_geometry,
    gaussian_rmp_to_rda_mean,
    histogram_rda,
    mean_fret_distance,
    model_distance,
    polynomial_transfer,
    standard_deviation_of_distances,
)
