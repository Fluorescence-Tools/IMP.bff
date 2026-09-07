/**
 *  \file IMP/bff/TopologyBuild.h
 *  \brief The one force-field system builder, from components to the typed
 *         value.
 *
 * `create_forcefield_system` is one function over the MOL2 reader, the
 * molecular graph, the improper expanders and the template reader: it reads
 * the components, derives sites/bonds/angles/dihedrals/impropers/groups and
 * assembles the system's JSON, then builds the typed system through
 * #IMP::bff::forcefield_system_from_json -- **the** conversion path, so every
 * route into a system has the same semantics.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_TOPOLOGYBUILD_H
#define IMPBFF_TOPOLOGYBUILD_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/Mol2IO.h>
#include <IMP/bff/ProbeLibrary.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! One component of a system to build: a MOL2, its template, and its role.
/*! The typed form of the `{"name": ..., "mol2": ..., "template": ...,
    "role": ...}` object the JSON builder takes. It exists so a caller says
    what it means in C++ instead of composing JSON by hand. */
struct IMPBFFEXPORT FFComponentSpec {
    std::string name;
    std::string mol2;
    //! A cgprobe component template; empty when the component has none.
    std::string template_path;
    //! `"fixed"` or `"mobile"`.
    std::string role;

    FFComponentSpec() {}
    FFComponentSpec(const std::string& name, const std::string& mol2,
                    const std::string& template_path, const std::string& role)
        : name(name), mol2(mol2), template_path(template_path), role(role) {}

    IMP_SHOWABLE_INLINE(FFComponentSpec,
                        out << "FFComponentSpec(" << name << ", " << role
                            << ")");
};
IMP_VALUES(FFComponentSpec, FFComponentSpecs);

//! The same builder, from typed component specs.
/*! Composes the JSON the other overload takes, so there is one conversion
    path rather than two. */
IMPBFFEXPORT ProbeForceFieldSystem create_forcefield_system(
        const std::vector<FFComponentSpec>& components,
        double bond_k = 2000.0, double angle_k = 400.0,
        double pi_dihedral_k = 12.0, double linker_dihedral_k = 1.5,
        double ring_improper_k = 40.0, double pi_improper_k = 180.0,
        double flat_improper_k = 120.0, double orient_improper_k = 220.0,
        int n_steps = 20000, int write_every = 100,
        double default_radius = 1.7, double default_mass = 12.0,
        double nonbonded_k = 5.0, double nonbonded_cutoff = 6.0,
        int minimize_steps = 200, const std::string& relative_to = "");

//! A system for a protein and a dye, both from MOL2 -- two components.
IMPBFFEXPORT ProbeForceFieldSystem create_probe_protein_system(
        const std::string& protein_mol2, const std::string& dye_mol2,
        const std::string& protein_name, const std::string& probe_name,
        const std::string& protein_template = "",
        const std::string& probe_template = "",
        double default_radius = 1.7, double default_mass = 12.0);

//! A system for one dye alone, with its backbone anchor held fixed.
/*! One mobile component and no protein. What this adds over the plain
    builder: the dye's `N`, `CA`, `C` and `O` sites are collected into a
    `<probe_name>_anchor` group and named in `fixed_groups`, so a sampler can
    hold the labelled residue's backbone still. */
IMPBFFEXPORT ProbeForceFieldSystem probe_forcefield_system(
        const std::string& dye_mol2, const std::string& probe_name = "dye",
        const std::string& probe_template = "", double default_radius = 1.7,
        double default_mass = 12.0);

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
IMPBFFEXPORT ProbeForceFieldSystem create_forcefield_system(
        const std::string& components_json,
        double bond_k = 2000.0, double angle_k = 400.0,
        double pi_dihedral_k = 12.0, double linker_dihedral_k = 1.5,
        double ring_improper_k = 40.0, double pi_improper_k = 180.0,
        double flat_improper_k = 120.0, double orient_improper_k = 220.0,
        int n_steps = 20000, int write_every = 100,
        double default_radius = 1.7, double default_mass = 12.0,
        double nonbonded_k = 5.0, double nonbonded_cutoff = 6.0,
        int minimize_steps = 200, const std::string& relative_to = "");

//! A topology-only system for one MOL2 component: sites, bonds, angles,
//! torsions.
/*!
    What #IntramolecularEnergy needs and no more: the connectivity that
    says which site pairs are 1-2, 1-3 or 1-4 and therefore excluded from the
    dye's own Lennard-Jones sum. No equilibrium lengths, no force constants --
    those come from a template, and a lone dye in vacuum has none.

    Site ids are `dye:<serial>:<atom_name>`.

    \note The sites carry the MOL2's **element**, and must: a site that reaches
    the LJ table elementless is parameterised by whatever the fallback is, so a
    hydrogen silently becomes something else. The element is in the file and in
    `FFSite`.

    \param[in] component the parsed MOL2 (see #read_mol2_component)
*/
IMPBFFEXPORT ProbeForceFieldSystem internal_topology_system(
        const Mol2Component& component);

//! `{serial: site atom name}` -- the MOL2 names deduplicated per component.
/*! Repeated names get an alpha suffix in serial order (`C`, `CA`, `CB`, ...,
    `CZ`, `CAA`), so site ids built from them are unique although MOL2 atom
    names are not (atto655: 83 atoms, 68 names). */
IMPBFFEXPORT std::map<int, std::string> serial_to_site_atom_names(
        const std::map<int, std::string>& serial_to_name);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_TOPOLOGYBUILD_H
