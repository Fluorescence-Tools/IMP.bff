"""DirectLabelingRestraint — fast attachment-atom distance scoring.

This restraint scores distances between **attachment atoms** directly,
**without** computing accessible volumes.  It is much faster than AV-based
restraints and is intended for use in MCMC sampling inner loops where
speed is critical.

The score is an asymmetric chi-squared with the same form used by
:class:`SimpleAVNetworkRestraint` — the difference is only what is
used as the "model distance".

Examples
--------
>>> from IMP.bff.restraints import DirectLabelingRestraint
>>> import numpy as np
>>>
>>> # 3 atoms: residues 2, 5, 7
>>> xyz = np.array([[0.,0.,0.], [10.,0.,0.], [20.,0.,0.]])
>>> r = DirectLabelingRestraint(xyz)
>>> r.add_site(residue_seq_number=0, atom_name="CA",
...            distance=55.0, error_neg=3.0, error_pos=5.0)
>>> r.add_site(residue_seq_number=1, atom_name="CA",
...            distance=65.0, error_neg=4.0, error_pos=6.0)
>>> score = r.evaluate()
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from IMP.bff.distance_metrics import chi2_score


@dataclass
class LabelingSite:
    """A single labeling site with an experimental distance to a partner.

    In practice paired distances are encoded by adding two
    :class:`LabelingSite` instances (one per residue) and letting the
    restraint score all unique pairs.  If a specific pair should be
    omitted, use :class:`AVMeasurement` instead.

    Parameters
    ----------
    residue_seq_number : int
        Residue sequence number (used as index into the coordinate array).
    atom_name : str
        Name of the attachment atom (used when a structured atom array
        is provided).
    distance : float
        Experimental distance (Å) — measured between this site and
        a partner site.
    error_neg, error_pos : float
        Asymmetric experimental errors (Å).
    forster_radius : float
        Förster radius (Å), currently unused by direct labeling.
    """

    residue_seq_number: int
    atom_name: str = "CB"
    distance: float = 0.0
    error_neg: float = 3.0
    error_pos: float = 5.0
    forster_radius: float = 52.0


class DirectLabelingRestraint:
    """Fast labeling restraint using attachment-atom distances.

    Parameters
    ----------
    xyz : (N, 3) ndarray
        Atomic coordinates in the order corresponding to
        *residue_seq_number* in the added sites.
    weight : float
        Overall restraint weight.

    Attributes
    ----------
    sites : list of LabelingSite
    """

    def __init__(
        self,
        xyz: Optional[np.ndarray] = None,
        weight: float = 1.0,
    ):
        self.xyz = np.asarray(xyz, dtype=np.float64) if xyz is not None else None
        self.weight = float(weight)
        self.sites: list[LabelingSite] = []

    def add_site(
        self,
        residue_seq_number: int,
        atom_name: str = "CB",
        distance: float = 0.0,
        error_neg: float = 3.0,
        error_pos: float = 5.0,
        forster_radius: float = 52.0,
    ) -> "DirectLabelingRestraint":
        """Add a labeling site.

        Parameters
        ----------
        residue_seq_number : int
            Row index into the coordinate array.
        atom_name, distance, error_neg, error_pos, forster_radius
            See :class:`LabelingSite`.

        Returns
        -------
        self
        """
        self.sites.append(LabelingSite(
            residue_seq_number=residue_seq_number,
            atom_name=atom_name,
            distance=distance,
            error_neg=error_neg,
            error_pos=error_pos,
            forster_radius=forster_radius,
        ))
        return self

    def evaluate(self, *args, **kwargs) -> float:
        r"""Evaluate the restraint: sum chi-squared over all site pairs.

        For each unique pair of sites :math:`(i, j)` the model distance
        is the Euclidean distance between the corresponding rows in
        *xyz*:

        .. math::

           d_{ij} = \\|\\mathbf{x}_i - \\mathbf{x}_j\\|

        The chi-squared contribution per pair is computed with
        :func:`IMP.bff.distance_metrics.chi2_score`.

        Returns
        -------
        float
            Weighted total.
        """
        if self.xyz is None or len(self.sites) < 2:
            return 0.0

        total = 0.0
        for i in range(len(self.sites)):
            si = self.sites[i]
            ri = si.residue_seq_number
            for j in range(i + 1, len(self.sites)):
                sj = self.sites[j]
                rj = sj.residue_seq_number
                dij = np.linalg.norm(self.xyz[ri] - self.xyz[rj])
                total += chi2_score(
                    dij, si.distance, si.error_neg, si.error_pos,
                )
        return total * self.weight

    def set_weight(self, w: float) -> None:
        """Set the overall restraint weight."""
        self.weight = float(w)

    def get_weight(self) -> float:
        """Get the overall restraint weight."""
        return self.weight

    def __repr__(self) -> str:
        return (
            f"DirectLabelingRestraint("
            f"n_sites={len(self.sites)})"
        )
