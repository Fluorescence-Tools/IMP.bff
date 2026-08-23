/**
 *  \file IMP/bff/TopologyBuild.h
 *  \brief The one force-field system builder, from components to the typed
 *         value.
 *
 * `build_forcefield_system` was 410 lines of Python orchestration over
 * already-C++ pieces (the MOL2 reader, the molecular graph, the improper
 * expanders, the template reader). It is one C++ function now: it reads the
 * components, derives sites/bonds/angles/dihedrals/impropers/groups, and
 * assembles the same JSON the Python dict held -- then builds the typed
 * system through #IMP::bff::forcefield_system_from_json, which is the one
 * conversion path and guarantees the same semantics the Python produced.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_TOPOLOGYBUILD_H
#define IMPBFF_TOPOLOGYBUILD_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/DyeForceField.h>

#include <string>

IMPBFF_BEGIN_NAMESPACE

//! The one force-field system builder, from components to the typed value.
/*!
    \param[in] components_json a JSON array of
               `{"name": ..., "mol2": path, "template": path?, "role":
               "fixed"|"mobile"}`. At most one fixed component; a dye alone
               has none.
    \param[in] bond_k,angle_k the harmonic force constants
    \param[in] pi_dihedral_k,linker_dihedral_k the torsion constants (pi:
               both central atoms C/N)
    \param[in] ring_improper_k,pi_improper_k,flat_improper_k,orient_improper_k
               the improper constants, one per template kind
    \param[in] n_steps,write_every,minimize_steps the sampling section
    \param[in] default_radius,default_mass site defaults, A and Da
    \param[in] nonbonded_k,nonbonded_cutoff the excluded-volume term
    \param[in] relative_to a directory the recorded mol2 paths are made
               relative to -- the writer passes the output file's directory.
               Empty records absolute paths.
    \return the typed system. Sites are `component/SITE` with SITE the atom
            name deduplicated by alpha suffix (C, C1 -> C, CA...) so site ids
            stay unique when MOL2 names repeat.
    \throw ValueException for a spec missing name/mol2/role, or for two fixed
           components
    \throw IOException when a MOL2 or template file cannot be read
*/
IMPBFFEXPORT DyeForceFieldSystem build_forcefield_system(
        const std::string& components_json,
        double bond_k = 2000.0, double angle_k = 400.0,
        double pi_dihedral_k = 12.0, double linker_dihedral_k = 1.5,
        double ring_improper_k = 40.0, double pi_improper_k = 180.0,
        double flat_improper_k = 120.0, double orient_improper_k = 220.0,
        int n_steps = 20000, int write_every = 100,
        double default_radius = 1.7, double default_mass = 12.0,
        double nonbonded_k = 5.0, double nonbonded_cutoff = 6.0,
        int minimize_steps = 200, const std::string& relative_to = "");

//! `{serial: site atom name}` -- the MOL2 names deduplicated per component.
/*! Repeated names get an alpha suffix in serial order (`C`, `CA`, `CB`, ...,
    `CZ`, `CAA`), so site ids built from them are unique although MOL2 atom
    names are not (atto655: 83 atoms, 68 names). */
IMPBFFEXPORT std::map<int, std::string> serial_to_site_atom_names(
        const std::map<int, std::string>& serial_to_name);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_TOPOLOGYBUILD_H
