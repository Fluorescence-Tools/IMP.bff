
 .. Labels:

======
Labels
======

``IMP.bff`` can simulated the positional distribution of biomolecular
labels around their attachment site. Such simulations are required to
compute for a structural model experimental observables that are used
for scoring structures based on fluorescence, epr, and x-link mass
spectrometry experiments.

In accessible volume (AV) computations the positional distribution is
computed for a coarse grained of the label. In the AV computations an
optimal (short) path from the attachment site of the label to other points
on a grid is computed (considering the dimensions of the label in the
construction of the grid).

Path search
-----------

``IMP.bff`` implements different (standard) path search algorithms that
can be used for implementing restraints. The Djikstra algorithm computes
for a start point the distance to all other grid points (till the distance
to the end point is found). The A* algorithm uses additionally a heuristic
to find a path from the start point to the end point. In the AV calculations
no end point is provided and the Djikstra algrothim is used to compute paths.

- Link to wiki for the two algorithm

Accessible volume
-----------------

AVs can be computed for all IMP particles with XYZ decorators. ``IMP.bff``
implements a AV decorator to label particles (for details see XXXX).

Scoring
-------

Networks

Literature
