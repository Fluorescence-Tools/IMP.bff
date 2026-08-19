"""Explicit all-atom dye modelling under a force field. **Legacy.**

This is molecular mechanics: build a topology for a dye, give it CHARMM36
parameters, attach it to a protein and sample it — Langevin dynamics, RRT,
Boltzmann and mean-field reweighting, kinetic trajectories. It is a real model
and it still runs. It is not what ``IMP.bff`` is *for*.

``IMP.bff`` is a forward-model engine from structure to experiment-neutral
fluorescence observables. Propagating an explicit dye under a potential is
IMP's own territory (``IMP.atom``, ``IMP.core``, ``IMP.isd``), and a
force-field layer specific to this package is a thing to build on IMP rather
than inside a spectroscopy library.

So as of 2026-08-18 this package is **outside the domain layout**:

* it is **not** a domain in :data:`IMP.bff.api.BY_DOMAIN`;
* its names are **not** on the flat ``IMP.bff.<Name>`` surface — 38 of them
  were removed;
* it is imported explicitly, by module path, or not at all.

Nothing in ``IMP.bff`` imports it. It is a leaf.

What was taken out of it first, because it was not molecular mechanics:

===========================================  =============================================
was                                          is
===========================================  =============================================
``cgdye.rotamer``                            ``IMP.bff.representation.rotamer``
``cgdye.sampling.rotamer``                   ``IMP.bff.sampling``
``cgdye.labeling``                           ``IMP.bff.label`` / ``.backbone_frame``
``cgdye.io.*`` (all of it)             ``IMP.bff.io.*``
``cgdye.utils``                              ``IMP.bff.tools``
``cgdye.topology.dye`` LJ table and kernel   ``IMP.bff.scoring``
===========================================  =============================================

A rotamer library is a **representation** — a list of states with weights, the
discrete sibling of an accessible volume — and the design has always said so.
It lived here because this package was the first thing that needed it. The same
is true of attaching a dye to a residue, which is the *system* layer, and of the
Lennard-Jones term, which the rotamer scorer needs and which a core domain must
not reach into a legacy package to get.

Reach what remains by its module path::

    from IMP.bff.cgdye.sampling.langevin import LangevinDyeSampler
    from IMP.bff.cgdye.topology.dye import build_dye_topology
"""

from __future__ import annotations

__all__: list[str] = []
