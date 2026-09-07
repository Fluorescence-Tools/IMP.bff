/**
 * \file ImpLayer.cpp
 * \brief The connection layer to IMP, compiled as one translation unit.
 *
 * IMP's module tooling compiles `src/*.cpp` and `src/internal/*.cpp` into
 * the module (tools/build/setup_all.py) and leaves any other subdirectory
 * alone -- `src/standalone/` relies on exactly that to stay out of the IMP
 * build. `src/imp/` needs the opposite: in the IMP build, and out of the
 * standalone one. So the layer's sources live in `src/imp/`, where `ls` shows
 * the whole cross-section, and this file is how the IMP build reaches them.
 * It is the only entry in `src/Files.cmake` for them, so a per-cpp build
 * compiles each once as well. The standalone build does not list this file.
 *
 * Order is alphabetical, as in the module's own unity file; nothing here
 * depends on it.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include "imp/AV.cpp"
#include "imp/AVMeanDistanceRestraint.cpp"
#include "imp/AVOccupancyMap.cpp"
#include "imp/Docking.cpp"
#include "imp/EmBridge.cpp"
#include "imp/FPSExport.cpp"
#include "imp/FPSProject.cpp"
#include "imp/HierarchyBridge.cpp"
#include "imp/Potentials.cpp"
#include "imp/ProbeAttachment.cpp"
#include "imp/ProbeDynamics.cpp"
#include "imp/ProbeNetworkRestraint.cpp"
