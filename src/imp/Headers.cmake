# The connection layer's public headers. IMP's module tooling links only
# include/*.h and include/internal/*.h into the build tree (tools/build/
# setup.py), so the layer's headers stay flat in include/ and this list is
# what names them. The standalone (IMP-free) build reads it to leave them
# out; the IMP build ignores it. Sources need no list: src/imp/ is globbed.
set(imp_bff_layer_headers
    "ProbeAccessibleVolumeDecorator.h;StructureReader.h;ProbeAccessibleVolumeOccupancyMap.h;ProbeNetworkRestraint.h;ProbeAccessibleVolumeMeanDistanceRestraint.h;ProbePotentialRestraints.h;Docking.h;MolecularProbeSimulation.h;ProbeDynamics.h;FPSProject.h;FPSExport.h;ProbeAttachment.h;IMPHierarchyBridge.h;IMPEMBridge.h;IMPAlgebraBridge.h;internal/AVLatticeState.h")
