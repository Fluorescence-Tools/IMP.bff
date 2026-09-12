/**
 *  \file IMP/bff/ProbePotentialTables.h
 *  \brief The parameter tables the coarse-grained potentials read, in one file.
 *
 * Every knowledge-based potential is a table and a loop; the loops are in
 * `ProbePotentialRestraints.h` and the tables are here. They arrived as four loose `.npy`
 * files in a chisurf directory -- `mj.npy`, `hb.npy`, `unres.npy`,
 * `rama_ala_pro_gly.npy` -- and they ship as **one PTO container**,
 * `data/potentials.pto`, which is what this module already does for
 * `dyes.drot.pto` and `dyes.mmfdb.pto`.
 *
 * One file rather than four because a potential is not one table: the UNRES
 * centroid term needs its own grid *and* the residue order that indexes it,
 * the Ramachandran map needs to say which of its channels is glycine, and a
 * reader that has to find four files in agreement with each other will one day
 * find three. A container carries them together with a manifest that says what
 * each is, where it came from, and what was done to it on the way in.
 *
 * Two kinds of payload:
 *
 *  - `pot.pmf` -- text in IMP's PMF format, which
 *    #IMP::core::StatisticalPairScore reads directly. Miyazawa-Jernigan and
 *    UNRES are these, so their potential *is* the file and there is no loop
 *    in this module that touches them.
 *  - `pot.grid` -- float64 in C order with the shape in the manifest. The
 *    hydrogen-bond lookup and the Ramachandran map are these, because their
 *    restraints are not pair scores over a distance.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PROBEPOTENTIALTABLES_H
#define IMPBFF_PROBEPOTENTIALTABLES_H

#include <IMP/bff/bff_config.h>

#include <IMP/bff/IMPCompatibility.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! One parameter table out of the container.
struct IMPBFFEXPORT PotentialTable {
    //! What the container calls it: `mj`, `unres`, `hbond`, `ramachandran`.
    std::string name;
    //! `pot.pmf` (text IMP's PMF reader takes) or `pot.grid` (float64).
    std::string kind;
    //! The grid's shape, C order; empty for a `pot.pmf`.
    std::vector<int> shape;
    //! The grid, flat, C order; empty for a `pot.pmf`.
    std::vector<double> values;
    //! The text; empty for a `pot.grid`.
    std::string text;

    PotentialTable() {}
    PotentialTable(std::string name, std::string kind)
        : name(name), kind(kind) {}

    //! The grid as a numpy view; reshape with #shape.
    void get_values(double** out_view, int* n_out_view) const;
    //! How many numbers the shape says there should be.
    int get_size() const;

    IMP_SHOWABLE_INLINE(PotentialTable, out << "PotentialTable(" << name
                                            << ", " << kind << ")");
};
IMP_VALUES(PotentialTable, PotentialTableList);

//! The shipped container, `data/potentials.pto`.
IMPBFFEXPORT std::string get_potential_container_path();

//! Every table in a container, by name.
/*! \param[in] path the container; empty reads the shipped one
    \throw IOException when it cannot be opened */
IMPBFFEXPORT std::vector<std::string> potential_table_names(
        std::string path = "");

//! The container's manifest, as JSON.
/*! What each table is, its shape and units, where it came from, and what the
    converter did to it -- which is the part a loose `.npy` cannot carry. */
IMPBFFEXPORT std::string read_potential_manifest(std::string path = "");

//! One table out of a container.
/*!
    \param[in] name `mj`, `unres`, `hbond` or `ramachandran`
    \param[in] path the container; empty reads the shipped one
    \throw IOException when the container has no table of that name
*/
IMPBFFEXPORT PotentialTable read_potential_table(std::string name,
                                                 std::string path = "");

//! Write a container: the tables, and a manifest that describes them.
/*!
    What `bin/imp_bff_potentials2pto` calls, and what a caller with tables of
    its own calls -- a potential this module does not ship is a file, not a
    patch.

    \param[in] path the container to write
    Each payload uses standard Zstd level 3 compression. Existing containers
    with size-prefixed Brotli payloads remain readable.

    \param[in] tables the tables; each is compressed on its own
    \param[in] manifest_json attached under `manifest.json`
    \throw IOException when the file cannot be written
*/
IMPBFFEXPORT void write_potential_tables(
        const std::string& path, const std::vector<PotentialTable>& tables,
        const std::string& manifest_json = "{}");

IMPBFF_END_NAMESPACE

#endif //IMPBFF_PROBEPOTENTIALTABLES_H
