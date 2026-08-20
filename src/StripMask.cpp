/**
 * \file StripMask.cpp
 * \brief The fps `strip_mask` dialect: which atoms a selection removes.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/StripMask.h>

#include <IMP/exception.h>

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <sys/stat.h>

IMPBFF_BEGIN_NAMESPACE

namespace strip {

std::string upper(const std::string& s) {
    std::string out(s);
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(out[i])));
    }
    return out;
}

std::string trim(const std::string& s) {
    const std::string space = " \t\n\r\f\v";
    const std::size_t a = s.find_first_not_of(space);
    if (a == std::string::npos) return std::string();
    return s.substr(a, s.find_last_not_of(space) - a + 1);
}

//! Split on ` and `, case-insensitively, which is the mask's only connective.
std::vector<std::string> split_and(const std::string& mask) {
    std::vector<std::string> out;
    const std::string haystack = upper(mask);
    std::size_t start = 0;
    while (true) {
        std::size_t at = std::string::npos;
        // A separator is the word "and" with whitespace on both sides; the
        // grammar has no other keyword, so this cannot swallow a name.
        for (std::size_t i = start; i + 4 < haystack.size() + 1; ++i) {
            if (haystack.compare(i, 5, " AND ") == 0) { at = i; break; }
        }
        if (at == std::string::npos) {
            out.push_back(trim(mask.substr(start)));
            break;
        }
        out.push_back(trim(mask.substr(start, at - start)));
        start = at + 5;
    }
    return out;
}

//! `keyword rest` when \p term starts with \p keyword, else false.
bool word_arg(const std::string& term, const char* keyword, std::string* rest) {
    const std::string key = upper(keyword);
    const std::string up = upper(term);
    if (up.compare(0, key.size(), key) != 0) return false;
    if (term.size() <= key.size()) return false;
    const char next = term[key.size()];
    if (next != ' ' && next != '\t') return false;
    *rest = trim(term.substr(key.size()));
    return !rest->empty();
}

bool all_digits(const std::string& s) {
    if (s.empty()) return false;
    std::size_t i = (s[0] == '-') ? 1 : 0;
    if (i >= s.size()) return false;
    for (; i < s.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
    }
    return true;
}

std::vector<std::string> split_plus(const std::string& s) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (start <= s.size()) {
        const std::size_t at = s.find('+', start);
        const std::string part =
                s.substr(start, at == std::string::npos ? std::string::npos
                                                        : at - start);
        if (!part.empty()) out.push_back(part);
        if (at == std::string::npos) break;
        start = at + 1;
    }
    return out;
}

//! (chain, resseq, name) of an ATOM/HETATM line; false when it is not one.
bool pdb_atom_fields(const std::string& line, std::string* chain, int* resseq,
                     std::string* name) {
    if (line.compare(0, 6, "ATOM  ") != 0 && line.compare(0, 6, "HETATM") != 0) {
        return false;
    }
    if (line.size() < 27) return false;
    const std::string res = trim(line.substr(22, 4));
    if (!all_digits(res)) return false;
    *resseq = std::atoi(res.c_str());
    *chain = trim(line.substr(21, 1));
    *name = upper(trim(line.substr(12, 4)));
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

StripSelection::StripSelection(const std::string& chain,
                               const std::vector<int>& resids,
                               const std::vector<std::string>& names,
                               bool has_names, bool negate)
    : chain_(chain), resids_(resids), has_names_(has_names), negate_(negate) {
    for (std::size_t i = 0; i < names.size(); ++i) {
        names_.insert(strip::upper(names[i]));
    }
}

std::vector<std::string> StripSelection::get_names() const {
    return std::vector<std::string>(names_.begin(), names_.end());
}

bool StripSelection::matches(const std::string& chain, int resseq,
                             const std::string& name) const {
    if (name.empty()) return false;
    if (!chain_.empty() && chain != chain_) return false;
    if (!resids_.empty()) {
        bool found = false;
        for (std::size_t i = 0; i < resids_.size(); ++i) {
            if (resids_[i] == resseq) { found = true; break; }
        }
        if (!found) return false;
    }
    if (has_names_) {
        const bool in_set = names_.count(strip::upper(name)) > 0;
        if (!in_set != negate_) return false;
    }
    return true;
}

StripSelection parse_strip_mask(const std::string& mask) {
    std::string chain;
    std::vector<int> resids;
    std::vector<std::string> names;
    bool has_names = false, negate = false, seen_name_term = false;

    const std::vector<std::string> terms = strip::split_and(mask);
    for (std::size_t t = 0; t < terms.size(); ++t) {
        const std::string term = terms[t];
        if (term.empty()) continue;
        std::string rest;

        if (strip::word_arg(term, "chain", &rest)) {
            // One token: `chain A or resid 3` is not `chain "A or resid 3"`,
            // it is a mask using a connective the dialect does not have.
            if (rest.find_first_of(" \t") != std::string::npos) {
                IMP_THROW("strip_mask '"
                              << mask
                              << "' is outside the fps dialect (chain <id> and "
                                 "resid <n> and [not] name A+B+...); evaluate "
                                 "the selection externally and hand the "
                                 "stripped structure to the consumer",
                          ValueException);
            }
            if (!chain.empty()) {
                IMP_THROW("strip_mask: repeated 'chain' term in '" << mask << "'",
                          ValueException);
            }
            chain = rest;
            continue;
        }
        if (strip::word_arg(term, "resid", &rest) ||
            strip::word_arg(term, "resi", &rest)) {
            if (rest.find_first_of(" \t") != std::string::npos ||
                !strip::all_digits(rest)) {
                IMP_THROW("strip_mask '"
                              << mask
                              << "' is outside the fps dialect ('resid " << rest
                              << "' is not a plain residue number)",
                          ValueException);
            }
            if (!resids.empty()) {
                IMP_THROW("strip_mask: repeated 'resid' term in '" << mask << "'",
                          ValueException);
            }
            resids.push_back(std::atoi(rest.c_str()));
            continue;
        }
        bool is_not = false;
        std::string name_rest;
        if (strip::word_arg(term, "not", &rest) &&
            strip::word_arg(rest, "name", &name_rest)) {
            is_not = true;
        } else if (!strip::word_arg(term, "name", &name_rest)) {
            IMP_THROW("strip_mask '"
                          << mask
                          << "' is outside the fps dialect (chain <id> and "
                             "resid <n> and [not] name A+B+...); evaluate the "
                             "selection externally and hand the stripped "
                             "structure to the consumer",
                      ValueException);
        }
        if (seen_name_term) {
            IMP_THROW("strip_mask: repeated name term in '" << mask << "'",
                      ValueException);
        }
        // A space-separated name list is a parse error in PyMOL and here: the
        // separator is `+`, and accepting a space would silently select one
        // atom where the document named several.
        if (name_rest.find(' ') != std::string::npos ||
            name_rest.find('\t') != std::string::npos) {
            IMP_THROW("strip_mask '" << mask
                          << "' is outside the fps dialect (the name separator "
                             "is '+', not a space)",
                      ValueException);
        }
        seen_name_term = true;
        negate = is_not;
        names = strip::split_plus(name_rest);
        if (names.empty()) {
            IMP_THROW("strip_mask: empty name list in '" << mask << "'",
                      ValueException);
        }
        has_names = true;
    }
    return StripSelection(chain, resids, names, has_names, negate);
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
    const StripSelection sel = parse_strip_mask(mask);
    const std::string keep_name = strip::upper(strip::trim(keep_atom_name));

    std::vector<std::string> kept;
    kept.reserve(lines.size());
    for (std::size_t i = 0; i < lines.size(); ++i) {
        std::string chain, name;
        int resseq = 0;
        if (strip::pdb_atom_fields(lines[i], &chain, &resseq, &name) &&
            sel.matches(chain, resseq, name)) {
            const bool is_attachment =
                    !keep_name.empty() && name == keep_name &&
                    resseq == keep_resseq &&
                    (keep_chain.empty() || chain == keep_chain);
            if (!is_attachment) continue;
        }
        kept.push_back(lines[i]);
    }
    return kept;
}

std::string stripped_pdb_for(const std::string& pdb_path,
                             const std::string& chain, int resseq,
                             const std::string& atom_name,
                             const std::string& strip_mask) {
    const std::string mask = strip::trim(strip_mask);

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
