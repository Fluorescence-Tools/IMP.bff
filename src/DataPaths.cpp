/**
 * \file DataPaths.cpp
 * \brief Where the shipped cgprobe data lives, and how to reach it.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/DataPaths.h>

#include <IMP/bff/Base.h>

#include <cerrno>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

IMPBFF_BEGIN_NAMESPACE

namespace data_paths {

//! Join two path components, tolerating a trailing or leading separator.
std::string join(const std::string& base, const std::string& subpath) {
    if (subpath.empty()) return base;
    if (base.empty()) return subpath;
    const bool base_ends = base[base.size() - 1] == '/';
    const bool sub_starts = subpath[0] == '/';
    if (base_ends && sub_starts) return base + subpath.substr(1);
    if (base_ends || sub_starts) return base + subpath;
    return base + "/" + subpath;
}

}  // namespace data_paths

std::string get_cgprobe_data_dir() { return IMP::bff::get_data_path("cgprobe"); }

std::string get_template_dir(std::string subpath) {
    return data_paths::join(data_paths::join(get_cgprobe_data_dir(), "templates"),
                            subpath);
}

std::string get_structure_dir(std::string subpath) {
    return data_paths::join(
        data_paths::join(get_cgprobe_data_dir(), "inputs/structures"), subpath);
}

std::string get_output_dir(std::string subpath) {
    char buffer[4096];
    const char* cwd = getcwd(buffer, sizeof(buffer));
    if (cwd == NULL) IMP_THROW("cannot read the working directory", IOException);
    return data_paths::join(data_paths::join(std::string(cwd), "output"), subpath);
}

std::string ensure_dir(std::string path) {
    // `mkdir -p`: walk the path making each missing component. An existing
    // directory is the expected case, not an error.
    std::string built;
    std::size_t start = 0;
    if (!path.empty() && path[0] == '/') {
        built = "/";
        start = 1;
    }
    while (start <= path.size()) {
        const std::size_t slash = path.find('/', start);
        const std::size_t end = slash == std::string::npos ? path.size() : slash;
        const std::string component = path.substr(start, end - start);
        if (!component.empty()) {
            built = data_paths::join(built, component);
            if (mkdir(built.c_str(), 0777) != 0 && errno != EEXIST) {
                IMP_THROW("cannot create " << built << ": " << std::strerror(errno),
                          IOException);
            }
        }
        if (slash == std::string::npos) break;
        start = slash + 1;
    }
    return path;
}

IMPBFF_END_NAMESPACE
