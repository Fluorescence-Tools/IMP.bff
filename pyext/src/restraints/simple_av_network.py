"""SimpleAVNetworkRestraint — programmatic AV network restraint.

This restraint scores a network of dye-accessible-volume pair distances
against experimental measurements using the asymmetric chi-squared
scoring function.  Unlike ``IMP.bff.AVNetworkRestraint`` it does **not**
require an ``fps.json`` file — AVs and measurements are added
programmatically.

Examples
--------
>>> from IMP.bff.representation.av import BasicAV
>>> from IMP.bff.restraints import SimpleAVNetworkRestraint
>>> import numpy as np
>>>
>>> # Create two AVs from point clouds
>>> pts1 = np.random.randn(200, 4).astype(np.float64)
>>> pts2 = np.random.randn(200, 4).astype(np.float64)
>>> av1 = BasicAV(points=pts1, position_name="donor")
>>> av2 = BasicAV(points=pts2, position_name="acceptor")
>>>
>>> # Build a restraint with one measurement
>>> r = SimpleAVNetworkRestraint()
>>> r.add_av_object("donor", av1)
>>> r.add_av_object("acceptor", av2)
>>> r.add_measurement(
...     av1_name="donor", av2_name="acceptor",
...     distance=55.0, error_neg=3.0, error_pos=5.0,
... )
>>> chi2 = r.evaluate()
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from IMP.bff.representation.av import BasicAV, compute_av
from IMP.bff.representation.distance import chi2_score


@dataclass
class AVMeasurement:
    """An experimental distance measurement between a pair of AVs.

    Parameters
    ----------
    av1_name : str
        Name of the first AV (must match a name passed to
        :meth:`SimpleAVNetworkRestraint.add_av_object`).
    av2_name : str
        Name of the second AV.
    distance : float
        Experimental distance (Å).
    error_neg : float
        Asymmetric negative error (Å); used when ``model < exp``.
    error_pos : float
        Asymmetric positive error (Å); used when ``model >= exp``.
    forster_radius : float
        Förster radius *R*\\ :sub:`0` (Å).  Used for
        ``distance_type="RDAMeanE"``.
    distance_type : str
        Which AV distance metric to compute:

        - ``"RDAMean"`` — mean inter-point distance
          :math:`\\langle R_{DA}\\rangle`
        - ``"RDAMeanE"`` — FRET-averaged distance :math:`R_E`
        - ``"Rmp"`` — mean-position distance
    """

    av1_name: str
    av2_name: str
    distance: float
    error_neg: float
    error_pos: float
    forster_radius: float = 52.0
    distance_type: str = "RDAMean"


class SimpleAVNetworkRestraint:
    """Programmatic AV-network chi-squared restraint.

    Parameters
    ----------
    name : str
        Optional label.

    Attributes
    ----------
    weight : float
        Overall weight applied to the total score.
    """

    def __init__(self, name: str = "SimpleAVNetworkRestraint"):
        self.name = name
        self.weight = 1.0
        self._avs: dict[str, BasicAV] = {}
        self._measurements: list[AVMeasurement] = []

    # ------------------------------------------------------------------
    # Adding AVs
    # ------------------------------------------------------------------
    def add_av_object(self, name: str, av: BasicAV) -> None:
        """Register a pre-computed accessible volume.

        Parameters
        ----------
        name : str
            Unique identifier.
        av : BasicAV
            Computed AV.
        """
        self._avs[name] = av

    def add_av_from_coords(
        self,
        name: str,
        atoms_xyz: np.ndarray,
        atoms_vdw: np.ndarray,
        source_xyz: np.ndarray,
        linker_length: float = 20.0,
        linker_width: float = 0.5,
        dye_radii: tuple[float, float, float] = (3.5, 0.0, 0.0),
        grid_resolution: float = 1.5,
    ) -> "SimpleAVNetworkRestraint":
        """Compute an AV from raw coordinates and register it.

        Parameters
        ----------
        name : str
            Unique identifier.
        atoms_xyz : (N, 3) array
        atoms_vdw : (N,) array
        source_xyz : (3,) array
        linker_length, linker_width, dye_radii, grid_resolution
            Forwarded to :func:`IMP.bff.representation.av.compute_av`.

        Returns
        -------
        self
        """
        result = compute_av(
            atoms_xyz, atoms_vdw, source_xyz,
            linker_length, linker_width, dye_radii, grid_resolution,
        )
        av = BasicAV(
            points=result.points,
            density=result.density,
            grid_origin=result.grid_origin,
            grid_step=result.grid_step,
            position_name=name,
        )
        self._avs[name] = av
        return self

    # ------------------------------------------------------------------
    # Adding measurements
    # ------------------------------------------------------------------
    def add_measurement(
        self,
        av1_name: str,
        av2_name: str,
        distance: float,
        error_neg: float,
        error_pos: float,
        forster_radius: float = 52.0,
        distance_type: str = "RDAMean",
    ) -> "SimpleAVNetworkRestraint":
        """Add an experimental distance measurement.

        Parameters
        ----------
        av1_name, av2_name : str
            AV identifiers (must have been registered).
        distance, error_neg, error_pos, forster_radius, distance_type
            See :class:`AVMeasurement`.

        Returns
        -------
        self
        """
        m = AVMeasurement(
            av1_name=av1_name,
            av2_name=av2_name,
            distance=distance,
            error_neg=error_neg,
            error_pos=error_pos,
            forster_radius=forster_radius,
            distance_type=distance_type,
        )
        self._measurements.append(m)
        return self

    # ------------------------------------------------------------------
    # Evaluation
    # ------------------------------------------------------------------
    def evaluate(self, derivatives: bool = False) -> float:
        """Compute the total chi-squared score.

        Parameters
        ----------
        derivatives : bool
            Currently ignored; reserved for future use.

        Returns
        -------
        float
            Weighted sum of chi-squared contributions.
        """
        total = 0.0
        for m in self._measurements:
            av1 = self._avs[m.av1_name]
            av2 = self._avs[m.av2_name]

            if m.distance_type == "RDAMean":
                model_dist = av1.dRDA(av2)
            elif m.distance_type == "RDAMeanE":
                model_dist = av1.dRDAE(av2, m.forster_radius)
            elif m.distance_type == "Rmp":
                model_dist = av1.dRmp(av2)
            else:
                raise ValueError(f"Unknown distance_type: {m.distance_type!r}")

            total += chi2_score(
                model_dist, m.distance, m.error_neg, m.error_pos,
            )
        return total * self.weight

    def get_model_distances(self) -> dict[str, float]:
        """Return a mapping ``measurement_name → model_distance``.

        The names are auto-generated as ``av1_name_av2_name``.
        """
        d = {}
        for i, m in enumerate(self._measurements):
            av1 = self._avs[m.av1_name]
            av2 = self._avs[m.av2_name]
            if m.distance_type == "RDAMean":
                val = av1.dRDA(av2)
            elif m.distance_type == "RDAMeanE":
                val = av1.dRDAE(av2, m.forster_radius)
            elif m.distance_type == "Rmp":
                val = av1.dRmp(av2)
            else:
                val = 0.0
            key = f"{m.av1_name}_{m.av2_name}"
            d[key] = val
        return d

    def set_weight(self, w: float) -> None:
        """Set the overall restraint weight."""
        self.weight = float(w)

    def get_weight(self) -> float:
        """Get the overall restraint weight."""
        return self.weight

    def __repr__(self) -> str:
        return (
            f"SimpleAVNetworkRestraint(name={self.name!r}, "
            f"n_avs={len(self._avs)}, "
            f"n_measurements={len(self._measurements)})"
        )
