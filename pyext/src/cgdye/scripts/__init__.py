"""Command-line drivers (``python -m IMP.bff.cgdye.scripts.<name>``).

Only three remain, and only because tests exercise them. The other 24 moved to
``junk/cgdye-scripts/`` in the PRD-113 cleanup (2026-08-18): they were drivers
for one system (``hgbp1_site481``), one-off analyses, or examples of a concept
that has since been retired. Every library function they called is still in
``IMP.bff.cgdye``; the one capability that existed *only* in a script -- writing
mol2 -- was lifted into :mod:`IMP.bff.cgdye.io.mol2` first.

A driver belongs next to the analysis it drove, not in an installed package. If
one of these grows a reusable function, that function belongs in a module.
"""
