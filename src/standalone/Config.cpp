/**
 * \file standalone/Config.cpp
 * \brief What IMP's module tooling generates into bff_config.cpp, for the standalone build.
 *
 * The data and example directories: an environment variable first
 * (IMP_BFF_DATA, IMP_BFF_EXAMPLES), then the directories the build was told
 * about (IMPBFF_STANDALONE_DATA_DIR / _EXAMPLE_DIR, set by CMake from the
 * install layout), then the working directory's data/ and examples/.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/bff_config.h>
#include <IMP/bff/Base.h>

#include <cstdlib>
#include <fstream>
#include <string>

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
std::string find_under(const char* env, const char* built_in, const std::string& file_name) {
    const char* e = std::getenv(env);
    const std::string dirs[2] = {e ? std::string(e) : std::string(), std::string(built_in)};
    for (int i = 0; i < 2; ++i) {
        if (dirs[i].empty()) continue;
        const std::string path = dirs[i] + "/" + file_name;
        std::ifstream probe(path.c_str());
        if (probe) return path;
    }
    IMP_THROW("Unable to find data file " << file_name << " under " << (e ? e : "(unset)")
              << " or " << built_in, IOException);
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
