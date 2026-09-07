/**
 *  \file IMP/bff/DataPaths.h
 *  \brief Where the shipped cgprobe data lives, and how to reach it.
 *
 * Templates, input structures and restraint files are IMP module *data*, not
 * package sources: they live in `imp.bff/data/cgprobe` and are reached through
 * #IMP::bff::get_data_path. Deriving them from a source file's location instead
 * ties them to where the package happens to sit, which is exactly what broke
 * when cgprobe moved out of imp-tricks.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_DATAPATHS_H
#define IMPBFF_DATAPATHS_H

#include <IMP/bff/bff_config.h>

#include <string>

IMPBFF_BEGIN_NAMESPACE

//! The cgprobe data directory.
IMPBFFEXPORT std::string get_cgprobe_data_dir();

//! The component-template directory, optionally with a subpath.
IMPBFFEXPORT std::string get_template_dir(std::string subpath = "");

//! The input-structure directory, optionally with a subpath.
IMPBFFEXPORT std::string get_structure_dir(std::string subpath = "");

//! The output directory, optionally with a subpath.
/*! Relative to the working directory, not to the installation: module data is
    read-only and installed, so nothing may be written beside it. */
IMPBFFEXPORT std::string get_output_dir(std::string subpath = "");

//! Create a directory and every parent it needs, and return it.
/*! An existing directory is not an error. */
IMPBFFEXPORT std::string ensure_dir(std::string path);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_DATAPATHS_H
