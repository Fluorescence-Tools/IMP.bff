/**
 * \file standalone/Config.cpp
 * \brief What IMP's module tooling generates into bff_config.cpp, for the standalone build.
 *
 * The data and example directories: an environment variable first
 * (IMP_BFF_DATA, IMP_BFF_EXAMPLES -- a list, separated as PATH is, ':' or
 * ';' on Windows; the wheel's Python side sets it to the data shipped in
 * the wheel plus the fetch cache), then the directory the build was told
 * about (IMPBFF_STANDALONE_DATA_DIR / _EXAMPLE_DIR, set by CMake from the
 * install layout, "data"/"examples" relative to the working directory when
 * it was told nothing).
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/bff_config.h>
#include <IMP/bff/IMPCompatibility.h>

#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#ifndef IMPBFF_STANDALONE_VERSION
#define IMPBFF_STANDALONE_VERSION "0.0.0"
#endif
#ifndef IMPBFF_STANDALONE_DATA_DIR
#define IMPBFF_STANDALONE_DATA_DIR "data"
#endif
#ifndef IMPBFF_STANDALONE_EXAMPLE_DIR
#define IMPBFF_STANDALONE_EXAMPLE_DIR "examples"
#endif

IMPBFF_BEGIN_NAMESPACE

namespace {
#ifdef _WIN32
const char kPathSep = ';';
#else
const char kPathSep = ':';
#endif

//! The directories to look in: the environment variable's list, then the built-in one.
std::vector<std::string> search_dirs(const char* env, const char* built_in) {
    std::vector<std::string> dirs;
    if (const char* e = std::getenv(env)) {
        std::string item;
        for (const char* c = e;; ++c) {
            if (*c == kPathSep || *c == '\0') {
                if (!item.empty()) dirs.push_back(item);
                item.clear();
                if (*c == '\0') break;
            } else {
                item += *c;
            }
        }
    }
    if (built_in && *built_in) dirs.push_back(built_in);
    return dirs;
}

std::string find_under(const char* env, const char* built_in, const std::string& file_name) {
    const std::vector<std::string> dirs = search_dirs(env, built_in);
    for (std::size_t i = 0; i < dirs.size(); ++i) {
        const std::string path = dirs[i] + "/" + file_name;
        std::ifstream probe(path.c_str());
        if (probe) return path;
    }
    std::string looked;
    for (std::size_t i = 0; i < dirs.size(); ++i) looked += (i ? ", " : "") + dirs[i];
    IMP_THROW("Unable to find data file " << file_name << " under " << looked
              << " (" << env << " is a " << kPathSep << "-separated list of directories)", IOException);
}
}  // namespace

std::string get_module_version() { return IMPBFF_STANDALONE_VERSION; }
std::string get_data_path(std::string file_name) {
    return find_under("IMP_BFF_DATA", IMPBFF_STANDALONE_DATA_DIR, file_name);
}
std::string get_example_path(std::string file_name) {
    return find_under("IMP_BFF_EXAMPLES", IMPBFF_STANDALONE_EXAMPLE_DIR, file_name);
}

IMPBFF_END_NAMESPACE
