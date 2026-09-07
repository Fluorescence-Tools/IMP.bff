/**
 * \file StripMask.cpp
 * \brief The fps `strip_mask` dialect: which atoms a selection removes.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/internal/Text.h>
#include <IMP/bff/StripMask.h>


#include <IMP/bff/Base.h>

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <sys/stat.h>

IMPBFF_BEGIN_NAMESPACE

namespace strip {

using internal::upper;

bool all_digits(const std::string& s) {
    if (s.empty()) return false;
    std::size_t i = (s[0] == '-') ? 1 : 0;
    if (i >= s.size()) return false;
    for (; i < s.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
    }
    return true;
}

//! (chain, resseq, name) of an ATOM/HETATM line; false when it is not one.
bool pdb_atom_fields(const std::string& line, std::string* chain, int* resseq,
                     std::string* name) {
    if (line.compare(0, 6, "ATOM  ") != 0 && line.compare(0, 6, "HETATM") != 0) {
        return false;
    }
    if (line.size() < 27) return false;
    const std::string res = internal::trimmed(line.substr(22, 4));
    if (!all_digits(res)) return false;
    *resseq = std::atoi(res.c_str());
    *chain = internal::trimmed(line.substr(21, 1));
    *name = upper(internal::trimmed(line.substr(12, 4)));
    return true;
}

//! Every field of a PDB ATOM/HETATM line a selection can name.
/*! A mask may ask about the residue name, the element, the alternate location
    or a distance, and a reader that only pulled out chain/resid/name would
    answer those with silence. */
bool pdb_selection_atom(const std::string& line, int index, SelectionAtom* out) {
    if (line.compare(0, 6, "ATOM  ") != 0 && line.compare(0, 6, "HETATM") != 0) {
        return false;
    }
    if (line.size() < 27) return false;
    const std::string res = internal::trimmed(line.substr(22, 4));
    if (!all_digits(res)) return false;

    out->index = index;
    out->hetatm = line.compare(0, 6, "HETATM") == 0;
    out->id = std::atoi(internal::trimmed(line.substr(6, 5)).c_str());
    out->name = upper(internal::trimmed(line.substr(12, 4)));
    out->alt = internal::trimmed(line.substr(16, 1));
    out->resn = upper(internal::trimmed(line.substr(17, 3)));
    out->chain = internal::trimmed(line.substr(21, 1));
    out->resi = std::atoi(res.c_str());
    const std::string icode = line.size() > 26
            ? internal::trimmed(line.substr(26, 1)) : std::string();
    if (!icode.empty()) out->resi_text = res + icode;
    if (line.size() >= 54) {
        out->x = std::atof(internal::trimmed(line.substr(30, 8)).c_str());
        out->y = std::atof(internal::trimmed(line.substr(38, 8)).c_str());
        out->z = std::atof(internal::trimmed(line.substr(46, 8)).c_str());
    }
    out->elem = line.size() >= 78
            ? upper(internal::trimmed(line.substr(76, 2))) : std::string();
    if (out->elem.empty()) {
        // No element column: the first letter of the name's leading alphabetic
        // run, which is what the rest of this module reads a MOL2 name as.
        for (std::size_t i = 0; i < out->name.size(); ++i) {
            if (std::isalpha(static_cast<unsigned char>(out->name[i]))) {
                out->elem = std::string(1, out->name[i]);
                break;
            }
        }
    }
    return true;
}

struct CacheKey {
    std::string path;
    long long mtime_ns, size;
    std::string chain, atom_name, mask;
    int resseq;
    bool operator<(const CacheKey& o) const {
        if (path != o.path) return path < o.path;
        if (mtime_ns != o.mtime_ns) return mtime_ns < o.mtime_ns;
        if (size != o.size) return size < o.size;
        if (chain != o.chain) return chain < o.chain;
        if (resseq != o.resseq) return resseq < o.resseq;
        if (atom_name != o.atom_name) return atom_name < o.atom_name;
        return mask < o.mask;
    }
};

std::map<CacheKey, std::string>& stripped_cache() {
    static std::map<CacheKey, std::string> cache;
    return cache;
}

}  // namespace strip

std::vector<std::string> backbone_atom_names() {
    static const char* const names[] = {"N", "CA", "C", "O"};
    return std::vector<std::string>(names, names + 4);
}

SelectionExpression parse_strip_mask(const std::string& mask) {
    return SelectionExpression(internal::trimmed(mask));
}

std::string site_strip_mask(const std::string& chain, int resseq,
                            const std::vector<std::string>& keep_atom_names) {
    std::vector<std::string> keep;
    for (std::size_t i = 0; i < keep_atom_names.size(); ++i) {
        const std::string name = strip::upper(keep_atom_names[i]);
        if (name.empty()) continue;
        bool seen = false;
        for (std::size_t k = 0; k < keep.size(); ++k) {
            if (keep[k] == name) { seen = true; break; }
        }
        if (!seen) keep.push_back(name);
    }

    std::ostringstream out;
    bool first = true;
    if (!chain.empty()) { out << "chain " << chain; first = false; }
    if (!first) out << " and ";
    out << "resid " << resseq;
    if (!keep.empty()) {
        out << " and not name ";
        for (std::size_t i = 0; i < keep.size(); ++i) {
            out << (i ? "+" : "") << keep[i];
        }
    }
    return out.str();
}

std::string default_strip_mask(const std::string& chain, int resseq,
                               const std::string& atom_name) {
    std::vector<std::string> keep = backbone_atom_names();
    keep.push_back(atom_name);
    return site_strip_mask(chain, resseq, keep);
}

std::vector<std::string> strip_pdb_lines(const std::vector<std::string>& lines,
                                         const std::string& mask,
                                         const std::string& keep_chain,
                                         int keep_resseq,
                                         const std::string& keep_atom_name) {
    const SelectionExpression sel = parse_strip_mask(mask);
    const std::string keep_name = strip::upper(internal::trimmed(keep_atom_name));

    // The whole file first: `byres`, `within` and `same chain as` are
    // questions about the structure, not about one line.
    std::vector<SelectionAtom> atoms;
    std::vector<std::size_t> line_of;
    atoms.reserve(lines.size());
    line_of.reserve(lines.size());
    for (std::size_t i = 0; i < lines.size(); ++i) {
        SelectionAtom atom;
        if (strip::pdb_selection_atom(lines[i], static_cast<int>(atoms.size()) + 1,
                                      &atom)) {
            atoms.push_back(atom);
            line_of.push_back(i);
        }
    }
    const std::vector<int> selected = sel.evaluate(atoms);

    std::vector<char> drop(lines.size(), 0);
    for (std::size_t k = 0; k < selected.size(); ++k) {
        if (!selected[k]) continue;
        const SelectionAtom& atom = atoms[k];
        const bool is_attachment =
                !keep_name.empty() && atom.name == keep_name &&
                atom.resi == keep_resseq &&
                (keep_chain.empty() || atom.chain == keep_chain);
        if (!is_attachment) drop[line_of[k]] = 1;
    }

    std::vector<std::string> kept;
    kept.reserve(lines.size());
    for (std::size_t i = 0; i < lines.size(); ++i) {
        if (!drop[i]) kept.push_back(lines[i]);
    }
    return kept;
}

std::string stripped_pdb_for(const std::string& pdb_path,
                             const std::string& chain, int resseq,
                             const std::string& atom_name,
                             const std::string& strip_mask) {
    const std::string mask = internal::trimmed(strip_mask);

    struct stat info;
    strip::CacheKey key;
    key.path = pdb_path;
    key.chain = chain;
    key.resseq = resseq;
    key.atom_name = atom_name;
    key.mask = mask;
    if (stat(pdb_path.c_str(), &info) == 0) {
#ifdef __APPLE__
        key.mtime_ns = (long long) info.st_mtimespec.tv_sec * 1000000000LL +
                       info.st_mtimespec.tv_nsec;
#else
        key.mtime_ns = (long long) info.st_mtim.tv_sec * 1000000000LL +
                       info.st_mtim.tv_nsec;
#endif
        key.size = (long long) info.st_size;
        std::map<strip::CacheKey, std::string>& cache = strip::stripped_cache();
        std::map<strip::CacheKey, std::string>::const_iterator hit =
                cache.find(key);
        if (hit != cache.end()) {
            struct stat unused;
            if (stat(hit->second.c_str(), &unused) == 0) return hit->second;
        }
    }

    std::ifstream source(pdb_path.c_str());
    if (!source) return pdb_path;
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(source, line)) lines.push_back(line + "\n");
    source.close();

    const std::string effective =
            mask.empty() ? default_strip_mask(chain, resseq, atom_name) : mask;
    const std::vector<std::string> kept =
            strip_pdb_lines(lines, effective, chain, resseq, atom_name);

    // `mkstemp` rather than tmpnam: the name is created and opened in one step,
    // so two AVs stripping in parallel cannot be handed the same path.
    std::string tmpl = "/tmp/bff_av_XXXXXX";
    std::vector<char> buffer(tmpl.begin(), tmpl.end());
    buffer.push_back('\0');
    const int fd = mkstemp(&buffer[0]);
    if (fd < 0) return pdb_path;
    const std::string out_path(&buffer[0]);
    close(fd);

    std::ofstream out(out_path.c_str());
    if (!out) return pdb_path;
    for (std::size_t i = 0; i < kept.size(); ++i) out << kept[i];
    out.close();

    strip::stripped_cache()[key] = out_path;
    return out_path;
}

IMPBFF_END_NAMESPACE
