/*
 * The Hierarchy and Particle overloads of the core's functions: every name
 * here existed before PRD-137 step 5 in the core header beside its array
 * twin, and keeps its Python spelling. Wrapped last, after every value type
 * the overloads return (ProteinFrame, StripReport, SelectionAtom, the
 * PathMap) is known.
 */
%include "IMP/bff/HierarchyBridge.h"

// PathMap.set_particles(ps) stays a method in Python; in C++ the lattice knows
// no particle and the layer installs a sphere source on it instead.
%extend IMP::bff::PathMap {
    void set_particles(const IMP::ParticlesTemp& ps) {
        IMP::bff::set_path_map_particles($self, ps);
    }
}
