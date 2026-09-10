/*
 * The inflated-sphere raster on a global cubic lattice, as an object.
 *
 * `OccupancyGrid` is the core of what `AVOccupancyMap` was (PRD-137 step 5):
 * the same counts, windows and deltas over a plain sphere list. The layer's
 * `AVOccupancyMap` (av.i) derives from it and adds the particle reading, so
 * this must be wrapped first for that class to keep the grid's methods.
 * The raw window readers into caller memory and the snapshot setter (a
 * shared_ptr) stay C++-only; `get_window()` is the Python reader.
 */
%ignore IMP::bff::OccupancyGrid::read_window;
%ignore IMP::bff::OccupancyGrid::read_window_strided;
%ignore IMP::bff::OccupancyGrid::set_coordinate_snapshot;
IMP_SWIG_OBJECT(IMP::bff, OccupancyGrid, OccupancyGrids);
%include "IMP/bff/OccupancyGrid.h"
