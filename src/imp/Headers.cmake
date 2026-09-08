# The connection layer's public headers. IMP's module tooling links only
# include/*.h and include/internal/*.h into the build tree (tools/build/
# setup.py), so the layer's headers stay flat in include/ and this list is
# what names them. The standalone (IMP-free) build reads it to leave them
# out; the IMP build ignores it. Sources need no list: src/imp/ is globbed.
set(imp_bff_layer_headers
    "AV.h;AVOccupancyMap.h;ProbeNetworkRestraint.h;AVMeanDistanceRestraint.h;Potentials.h;Docking.h;DyeDynamics.h;ProbeDynamics.h;FPSProject.h;FPSExport.h;ProbeAttachment.h;HierarchyBridge.h;EmBridge.h;AlgebraBridge.h;internal/AVLatticeState.h")
