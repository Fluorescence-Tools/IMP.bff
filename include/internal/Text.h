/**
 *  \file IMP/bff/internal/Text.h
 *  \brief The string and file helpers shared by this module's readers.
 *
 *  Extension dispatch must not depend on which file the caller sits in, so
 *  `ends_with` is case-insensitive here and everywhere: `LIB.BCIF` is a
 *  trajectory and `X.PDB` is a structure.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_INTERNAL_TEXT_H
#define IMPBFF_INTERNAL_TEXT_H

#include <IMP/bff/bff_config.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <vector>
#include <limits>
#include <string>
#include <sys/stat.h>
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <direct.h>
#else
#include <dirent.h>
#endif

IMPBFF_BEGIN_INTERNAL_NAMESPACE

//! Case-insensitive suffix test, for dispatching on a file extension.
/*! \p tail is compared in lower case, so it must be written that way:
    `ends_with(path, ".pdb")` matches `X.pdb`, `X.PDB` and `X.Pdb`. */
inline bool ends_with(const std::string& s, const std::string& tail) {
    if (s.size() < tail.size()) return false;
    for (std::size_t i = 0; i < tail.size(); ++i) {
        const char c = static_cast<char>(
                std::tolower(static_cast<unsigned char>(
                        s[s.size() - tail.size() + i])));
        if (c != tail[i]) return false;
    }
    return true;
}

//! `s` without leading or trailing whitespace; empty when it is all space.
/*! The whitespace set is the C one (` \t\n\r\f\v`), which is the wider of the
    two that were in use -- the readers that trimmed only ` \t\r\n` cannot see
    a form feed or a vertical tab in the formats they parse. */
inline std::string trimmed(const std::string& s) {
    static const char* space = " \t\n\r\f\v";
    const std::size_t a = s.find_first_not_of(space);
    if (a == std::string::npos) return std::string();
    return s.substr(a, s.find_last_not_of(space) - a + 1);
}

//! Whether the path names something that can be stat'ed.
inline bool file_exists(const std::string& path) {
    struct stat info;
    return stat(path.c_str(), &info) == 0;
}

//! Create a directory (and its parents); no error when it is already there.
inline void make_directory(const std::string& path) {
    if (path.empty() || file_exists(path)) return;
    const std::size_t slash = path.find_last_of('/');
    if (slash != std::string::npos && slash > 0) {
        make_directory(path.substr(0, slash));
    }
#ifdef _WIN32
    _mkdir(path.c_str());
#else
    ::mkdir(path.c_str(), 0755);
#endif
}

//! The files directly in \p path with the given lower-case suffix, sorted.
/*! Empty when \p path is not a directory, which is how a caller tells a
    directory of structures from one structure. \p path may be empty, which
    lists the working directory. */
inline std::vector<std::string> directory_entries(const std::string& path,
                                                  const std::string& suffix) {
    const std::string dir = path.empty() ? std::string(".") : path;
    std::vector<std::string> out;
#ifdef _WIN32
    WIN32_FIND_DATAA data;
    HANDLE handle = FindFirstFileA((dir + "\\*").c_str(), &data);
    if (handle == INVALID_HANDLE_VALUE) return out;
    do {
        const std::string name = data.cFileName;
        if (name == "." || name == "..") continue;
        if (!ends_with(name, suffix)) continue;
        out.push_back(path.empty() ? name : path + "/" + name);
    } while (FindNextFileA(handle, &data));
    FindClose(handle);
#else
    DIR* d = opendir(dir.c_str());
    if (d == NULL) return out;
    for (struct dirent* entry = readdir(d); entry != NULL;
         entry = readdir(d)) {
        const std::string name = entry->d_name;
        if (name == "." || name == "..") continue;
        if (!ends_with(name, suffix)) continue;
        out.push_back(path.empty() ? name : path + "/" + name);
    }
    closedir(d);
#endif
    std::sort(out.begin(), out.end());
    return out;
}

//! Upper-cased, byte by byte. Atom, residue and chain names are compared
//! case-insensitively wherever a file may disagree with a selector.
inline std::string upper(std::string s) {
    for (auto& c : s) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
}

//! A quiet NaN, spelled once.
inline double nan_value() {
    return std::numeric_limits<double>::quiet_NaN();
}

IMPBFF_END_INTERNAL_NAMESPACE

#endif /* IMPBFF_INTERNAL_TEXT_H */
