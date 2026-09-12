/**
 * \file ProbeRotamerLibrary.cpp
 * \brief The one rotamer library value.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/ProbeRotamerLibrary.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/MMFDBProfile.h>
#include <IMP/bff/internal/ptolib.h>
#include <IMP/bff/ZMatrix.h>
#include <IMP/bff/IMPCompatibility.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

// `out_view` is the module's *managed* out-view name (`types.i`): numpy takes
// the buffer and frees it, so a getter hands over a copy it has malloc'ed and
// never the vector's own storage -- `copy_to_view` is that, and handing out
// `coords.data()` here instead aborts the interpreter the first time numpy
// collects one of these arrays.
void ProbeRotamerLibrary::get_coords(double** out_view, int* n_out_view) const {
    internal::copy_to_view(coords, out_view, n_out_view);
}

void ProbeRotamerLibrary::get_weights(double** out_view, int* n_out_view) const {
    internal::copy_to_view(weights, out_view, n_out_view);
}

void normalize_weights_in_place(std::vector<double>& weights) {
    double total = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i) total += weights[i];
    if (total > 0.0) {
        for (std::size_t i = 0; i < weights.size(); ++i) weights[i] /= total;
    } else if (!weights.empty()) {
        const double w = 1.0 / static_cast<double>(weights.size());
        for (std::size_t i = 0; i < weights.size(); ++i) weights[i] = w;
    }
}

ProbeRotamerLibrary normalize_probe_rotamer_weights(const ProbeRotamerLibrary& lib) {
    ProbeRotamerLibrary out = lib;
    normalize_weights_in_place(out.weights);
    return out;
}

// ===========================================================================
// The `.drot` container: reading.
// ===========================================================================

namespace {

typedef std::map<std::string, std::string> JsonFlat;

// -- layer 1: brotli --------------------------------------------------------

//! Whole-file decompression, streamed: the decoded size is not known up front.
/*! **Not the one-shot API.** The one-shot decoder reports a too-small output
    buffer as an error, not as "needs more room" -- so a reader that guesses
    the size from the compressed length, as this one did with `size * 8`,
    fails outright on anything that compresses better than eight times. Every
    shipped `.drot` is under that ratio, which is why it worked; a library of
    near-identical conformers is not, and would have been the first file to
    be reported as corrupt while being perfectly good. A raw size of 0 sends
    ptolib's codec down the same streaming path: it grows as it goes, so
    there is nothing to guess and nothing to re-decode. */
std::vector<unsigned char> brotli_decompress(const unsigned char* data,
                                             std::size_t size) {
    std::vector<unsigned char> out;
    if (!pto::decompress_bytes("brotli", data, size, 0, out)) {
        IMP_THROW("read_drot: corrupt brotli stream", IOException);
    }
    return out;
}

// -- layer 2: tar -----------------------------------------------------------

//! The ustar subset: one member's (name, data) or an empty name at the end.
struct TarMember {
    std::string name;
    std::vector<unsigned char> data;
};

std::map<std::string, std::vector<unsigned char> > untar(
        const std::vector<unsigned char>& bytes) {
    std::map<std::string, std::vector<unsigned char> > out;
    std::size_t off = 0;
    while (off + 512 <= bytes.size()) {
        const unsigned char* h = &bytes[off];
        if (h[0] == '\0') break;  // end-of-archive block
        std::string name(reinterpret_cast<const char*>(h));
        name = name.substr(0, name.find('\0'));
        std::size_t size = 0;
        for (int i = 124; i < 136 && h[i] != ' ' && h[i] != '\0'; ++i) {
            size = size * 8 + (h[i] - '0');
        }
        const char type = h[156];
        off += 512;
        if (type == '0' || type == '\0') {  // regular file
            if (off + size > bytes.size()) {
                IMP_THROW("read_drot: tar member runs past the archive",
                          IOException);
            }
            out[name] = std::vector<unsigned char>(bytes.begin() + off,
                                                   bytes.begin() + off + size);
        }
        off += (size + 511) / 512 * 512;
    }
    return out;
}

// -- layer 3: minimal JSON --------------------------------------------------

//! Recursive-descent JSON for the writer's grammar (see file comment).
struct Json {
    // Scalars are kept as strings; arrays of numbers are parsed on demand.
    typedef std::map<std::string, Json> Object;

    enum Kind { SCALAR, ARRAY, OBJECT } kind;
    std::string scalar;
    std::vector<Json> array;
    std::map<std::string, Json> object;

    Json() : kind(SCALAR) {}
    explicit Json(Kind k) : kind(k) {}

    double as_double() const {
        return kind == SCALAR ? std::atof(scalar.c_str()) : 0.0;
    }
    int as_int() const {
        return kind == SCALAR ? static_cast<int>(std::atol(scalar.c_str()))
                              : 0;
    }
};

class JsonParser {
    const char* p_;
    const char* end_;

    void skip_ws() {
        while (p_ < end_ && std::isspace(static_cast<unsigned char>(*p_))) {
            ++p_;
        }
    }
    char peek() {
        skip_ws();
        if (p_ >= end_) IMP_THROW("read_drot: JSON ends early", IOException);
        return *p_;
    }
    void expect(char c) {
        if (peek() != c) IMP_THROW("read_drot: unexpected JSON token",
                                   IOException);
        ++p_;
    }

    void parse_string(std::string& out) {
        expect('"');
        out.clear();
        while (p_ < end_ && *p_ != '"') {
            if (*p_ == '\\' && p_ + 1 < end_) {
                ++p_;
                switch (*p_) {
                    case 'n': out += '\n'; break;
                    case 't': out += '\t'; break;
                    default: out += *p_; break;
                }
            } else {
                out += *p_;
            }
            ++p_;
        }
        expect('"');
    }

    void parse_value(Json& v) {
        char c = peek();
        if (c == '{') {
            v.kind = Json::OBJECT;
            ++p_;
            if (peek() == '}') { ++p_; return; }
            while (true) {
                std::string key;
                parse_string(key);
                expect(':');
                parse_value(v.object[key]);
                if (peek() == ',') { ++p_; continue; }
                expect('}');
                break;
            }
        } else if (c == '[') {
            v.kind = Json::ARRAY;
            ++p_;
            if (peek() == ']') { ++p_; return; }
            while (true) {
                v.array.push_back(Json());
                parse_value(v.array.back());
                if (peek() == ',') { ++p_; continue; }
                expect(']');
                break;
            }
        } else if (c == '"') {
            v.kind = Json::SCALAR;
            parse_string(v.scalar);
        } else {
            // number / true / false / null -- kept verbatim
            const char* start = p_;
            while (p_ < end_ && *p_ != ',' && *p_ != '}' && *p_ != ']' &&
                   !std::isspace(static_cast<unsigned char>(*p_))) {
                ++p_;
            }
            if (p_ == start) {
                IMP_THROW("read_drot: empty JSON value", IOException);
            }
            v.kind = Json::SCALAR;
            v.scalar.assign(start, p_);
        }
    }

 public:
    explicit JsonParser(const std::vector<unsigned char>& bytes)
        : p_(reinterpret_cast<const char*>(bytes.data())),
          end_(p_ + bytes.size()) {}
    Json parse() {
        Json v;
        parse_value(v);
        return v;
    }
};

//! Access with a clear error rather than a silent default.
const Json& json_get(const Json& obj, const std::string& key) {
    Json::Object::const_iterator it = obj.object.find(key);
    if (it == obj.object.end()) {
        IMP_THROW("read_drot: drot.json missing '" << key << "'", IOException);
    }
    return it->second;
}

// -- varint weights ---------------------------------------------------------

std::vector<double> decode_varint(const std::vector<unsigned char>& bytes,
                                  std::size_t count) {
    std::vector<double> out(count, 0.0);
    std::size_t i = 0;
    for (std::size_t k = 0; k < count; ++k) {
        unsigned long long v = 0;
        int shift = 0;
        while (i < bytes.size()) {
            const unsigned char b = bytes[i++];
            v |= static_cast<unsigned long long>(b & 0x7F) << shift;
            shift += 7;
            if (!(b & 0x80)) break;
        }
        out[k] = static_cast<double>(v);
    }
    if (i != bytes.size()) {
        IMP_THROW("read_drot: trailing bytes in weights.bin", IOException);
    }
    return out;
}

// -- grids ------------------------------------------------------------------

//! How one grid member is stored: which file, which numbers, which order.
struct GridSpec {
    std::string file;
    std::string dtype;   //!< "f32" or "i16"
    double scale;        //!< counts -> unit, for the integer rungs
    bool shuffle;        //!< byte-plane shuffled (the float rung is)
    bool col_major;      //!< value [dof][rotamer] rather than [rotamer][dof]
    GridSpec() : dtype("i16"), scale(1.0), shuffle(false), col_major(true) {}
};

//! v5/v7/v8 spell a grid `"<file>:<step>:col|row"`; v9 spells it as an object.
GridSpec parse_grid(const Json& grids, const std::string& key, int version) {
    GridSpec g;
    const Json& v = json_get(grids, key);
    if (version >= 9) {
        g.file = json_get(v, "file").scalar;
        g.dtype = json_get(v, "dtype").scalar;
        Json::Object::const_iterator it = v.object.find("scale");
        g.scale = it == v.object.end() ? 1.0 : it->second.as_double();
        it = v.object.find("shuffle");
        g.shuffle = it != v.object.end() && it->second.as_int() != 0;
        it = v.object.find("order");
        g.col_major = it == v.object.end() || it->second.scalar != "row";
        if (g.dtype != "f32" && g.dtype != "i16") {
            IMP_THROW("read_drot: unknown grid dtype '" << g.dtype << "'",
                      IOException);
        }
        return g;
    }
    std::istringstream ss(v.scalar);
    std::string step, order;
    std::getline(ss, g.file, ':');
    std::getline(ss, step, ':');
    std::getline(ss, order, ':');
    g.scale = step.empty() ? 1.0 : std::atof(step.c_str());
    g.col_major = (version == 7 || version == 8);
    return g;
}

//! One grid member as doubles, `count` of them, in the file's own order.
std::vector<double> read_grid(
        const std::map<std::string, std::vector<unsigned char> >& members,
        const GridSpec& g, std::size_t count, const std::string& key) {
    std::map<std::string, std::vector<unsigned char> >::const_iterator it =
            members.find(g.file);
    if (it == members.end()) {
        IMP_THROW("read_drot: missing grid member '" << g.file << "'",
                  IOException);
    }
    const std::vector<unsigned char>& bytes = it->second;
    const std::size_t item = g.dtype == "f32" ? 4 : 2;
    if (bytes.size() != count * item) {
        IMP_THROW("read_drot: " << key << " grid holds " << bytes.size()
                  << " bytes, the header asks for " << count * item,
                  IOException);
    }
    std::vector<unsigned char> plain;
    const unsigned char* raw = bytes.data();
    if (g.shuffle) {  // byte-plane shuffled: all first bytes, then all second
        plain.resize(bytes.size());
        for (std::size_t b = 0; b < item; ++b) {
            for (std::size_t i = 0; i < count; ++i) {
                plain[i * item + b] = bytes[b * count + i];
            }
        }
        raw = plain.data();
    }
    std::vector<double> out(count);
    for (std::size_t i = 0; i < count; ++i) {
        if (item == 4) {
            float v;
            std::memcpy(&v, raw + 4 * i, 4);
            out[i] = v;
        } else {
            const short v = static_cast<short>(static_cast<unsigned short>(
                    raw[2 * i] | (raw[2 * i + 1] << 8)));
            out[i] = v * g.scale;
        }
    }
    return out;
}

}  // namespace

ProbeRotamerLibrary read_probe_rotamer_drot(const std::string& path, const std::string& library) {
    // layer 1+2: the envelope, whichever one this file uses. v10 rides the
    // shared PTO container (one attached object per member, each compressed
    // on its own); v5-v9 are one brotli stream over a tar. Both hand the
    // rest of this function the same map, so only these lines know.
    //
    // Which one it is takes four octets, and asking that way matters: a
    // family container is twenty megabytes holding a hundred libraries, and
    // reading all of it to identify one -- which is what this did while the
    // only envelope was a tar that had to be read whole anyway -- cost more
    // than everything else here put together.
    std::map<std::string, std::vector<unsigned char> > members;
    if (pto::is_pto_file(path)) {
        pto::File pto;
        if (!pto.open(path, false)) {
            IMP_THROW("PTO: cannot open container " << path << ": "
                      << (pto.error().empty() ? "not an EBML document with DocType pto"
                                                : pto.error()), IOException);
        }
        const std::vector<pto::PtoObject>& objects = pto.objects();
        // In a family container every object is named `<library>/<member>`;
        // in a single-library one the member name stands alone. Only the
        // asked-for library's objects are read -- that is the point of the
        // envelope, and of walking it by seeking.
        const std::string prefix = library.empty() ? std::string()
                                                   : library + "/";
        bool found = false;
        for (std::size_t i = 0; i < objects.size(); ++i) {
            const std::string& name = objects[i].name;
            if (!prefix.empty()) {
                if (name.compare(0, prefix.size(), prefix) != 0) continue;
            } else if (name.find('/') != std::string::npos) {
                continue;                     // a bundle, and none was asked
            }
            found = true;
            // ptolib's reader decodes by the object's encoding -- the bytes
            // here are the payload, not the stream. (Since the container
            // moved to ptolib this used to decompress a second time, and
            // every .drot read failed with "corrupt brotli stream".)
            members[name.substr(prefix.size())] = pto.read(objects[i].uid);
        }
        if (!found) {
            std::ostringstream have;
            for (std::size_t i = 0, shown = 0; i < objects.size(); ++i) {
                const std::size_t slash = objects[i].name.find('/');
                if (slash == std::string::npos) continue;
                const std::string lib = objects[i].name.substr(0, slash);
                if (have.str().find(lib) != std::string::npos) continue;
                if (shown++ == 6) { have << ", ..."; break; }
                have << (shown > 1 ? ", " : "") << lib;
            }
            IMP_THROW("read_drot: " << path << " holds no library '" << library
                      << "'" << (have.str().empty() ? std::string()
                                 : std::string("; it holds ") + have.str()),
                      IOException);
        }
    } else {
        // The older envelope is one brotli stream over a tar, so it has to be
        // read whole; nothing shipped is in this shape any more.
        FILE* f = std::fopen(path.c_str(), "rb");
        if (f == nullptr) {
            IMP_THROW("read_drot: cannot open " << path, IOException);
        }
        std::vector<unsigned char> raw;
        char buf[1 << 16];
        std::size_t n;
        while ((n = std::fread(buf, 1, sizeof(buf), f)) > 0) {
            raw.insert(raw.end(), buf, buf + n);
        }
        std::fclose(f);
        const std::vector<unsigned char> payload =
                brotli_decompress(raw.data(), raw.size());
        members = untar(payload);
    }
    if (members.find("drot.json") == members.end()) {
        IMP_THROW("read_drot: no drot.json in " << path
                  << " (not a .drot container?)", IOException);
    }

    // layer 3: header
    const Json header = JsonParser(members.at("drot.json")).parse();
    const int version = json_get(header, "version").as_int();
    if (version != 5 && version != 7 && version != 8 && version != 9 &&
        version != 10) {
        IMP_THROW("read_drot: version " << version
                  << " not supported here (5, 7, 8, 9, 10)", IOException);
    }
    const int n_atoms = json_get(header, "n_atoms").as_int();
    const int n_rot = json_get(header, "n_rotamers").as_int();
    const int n_rows = json_get(header, "n_rows").as_int();
    const int n_base = json_get(header, "n_base").as_int();
    if (n_atoms < 1 || n_rot < 1 || n_rows < 1 || n_base < 3 ||
        n_rows > n_atoms || n_base > n_atoms) {
        IMP_THROW("read_drot: inconsistent shapes in drot.json", IOException);
    }

    // grids: how each member stores its numbers. v5-v8 wrote one string per
    // grid ("<file>:<step>:col|row") and took bond lengths from the template;
    // v9 writes an object per grid and adds the per-conformer bond lengths,
    // which is what makes it lossless (PRD-118).
    const Json& grids = json_get(header, "grids");
    const GridSpec base_grid = parse_grid(grids, "base", version);
    const GridSpec theta_grid = parse_grid(grids, "theta", version);
    const GridSpec phi_grid = parse_grid(grids, "phi", version);
    const bool has_lens = version >= 9 &&
            grids.object.find("lens") != grids.object.end();
    const GridSpec lens_grid = has_lens ? parse_grid(grids, "lens", version)
                                        : GridSpec();

    // base offset
    std::vector<int> base_off;
    {
        const Json& arr = json_get(header, "base_offset");
        for (std::size_t i = 0; i < arr.array.size(); ++i) {
            const int v = arr.array[i].as_int();
            if (v < 0 || v >= n_atoms) {
                IMP_THROW("read_drot: base_offset out of range", IOException);
            }
            base_off.push_back(v);
        }
        if (static_cast<int>(base_off.size()) != n_base) {
            IMP_THROW("read_drot: base_offset count != n_base", IOException);
        }
    }

    // rows: the Z-matrix tree as written by the encoder
    std::vector<int> rows;
    {
        const Json rj = JsonParser(members.at("rows.json")).parse();
        for (std::size_t i = 0; i < rj.array.size(); ++i) {
            rows.push_back(rj.array[i].as_int());
        }
        if (rows.size() != static_cast<std::size_t>(4 * n_rows)) {
            IMP_THROW("read_drot: rows.json size != 4 * n_rows", IOException);
        }
        for (std::size_t i = 0; i < rows.size(); ++i) {
            if (rows[i] < 0 || rows[i] >= n_atoms) {
                IMP_THROW("read_drot: row index out of range", IOException);
            }
        }
    }

    // template: atom name [resname] element x y z
    std::vector<std::string> names, resnames, elements;
    std::vector<double> tpl;
    {
        std::istringstream cif(std::string(
                members.at("template.cif").begin(),
                members.at("template.cif").end()));
        std::string line;
        while (std::getline(cif, line)) {
            if (line.empty() || line[0] == '#' || line[0] == '_' ||
                line.rfind("data_", 0) == 0 || line == "loop_") {
                continue;
            }
            std::istringstream ls(line);
            std::vector<std::string> tok;
            std::string t;
            while (ls >> t) tok.push_back(t);
            if (tok.size() == 6) {
                names.push_back(tok[0]); resnames.push_back(tok[1]);
                elements.push_back(tok[2]);
                tpl.push_back(std::atof(tok[3].c_str()));
                tpl.push_back(std::atof(tok[4].c_str()));
                tpl.push_back(std::atof(tok[5].c_str()));
            } else if (tok.size() == 5) {
                names.push_back(tok[0]); elements.push_back(tok[1]);
                tpl.push_back(std::atof(tok[2].c_str()));
                tpl.push_back(std::atof(tok[3].c_str()));
                tpl.push_back(std::atof(tok[4].c_str()));
            }
        }
        if (static_cast<int>(names.size()) != n_atoms) {
            IMP_THROW("read_drot: template has " << names.size()
                      << " atoms, header says " << n_atoms, IOException);
        }
        if (!resnames.empty() &&
            static_cast<int>(resnames.size()) != n_atoms) {
            resnames.clear();  // partial column: ignore rather than guess
        }
    }

    // layer 4: grids -> exact reconstruction
    ProbeRotamerLibrary out;
    out.n_rotamers = n_rot;
    out.n_atoms = n_atoms;
    out.atom_names = names;
    out.resnames = resnames;
    out.elements = elements;

    // Bond lengths: v9 stores one per row per conformer, which is what makes
    // the store lossless -- MD flexes bonds by ~0.004 A RMS and a torsion-only
    // reconstruction inherits that as 0.016 A of coordinate error, enough to
    // move a FRET efficiency in the third decimal. v5-v8 have no such grid, so
    // there the template's lengths are all there is.
    std::vector<double> tpl_lens(n_rows);
    for (int r = 0; r < n_rows; ++r) {
        const int i = rows[4 * r], p1 = rows[4 * r + 1];
        const double* xi = &tpl[3 * i];
        const double* x1 = &tpl[3 * p1];
        tpl_lens[r] = std::sqrt((xi[0] - x1[0]) * (xi[0] - x1[0]) +
                                (xi[1] - x1[1]) * (xi[1] - x1[1]) +
                                (xi[2] - x1[2]) * (xi[2] - x1[2]));
    }

    const std::size_t base_count =
            static_cast<std::size_t>(n_rot) * n_base * 3;
    const std::size_t grid_count =
            static_cast<std::size_t>(n_rot) * n_rows;
    const std::vector<double> gbase =
            read_grid(members, base_grid, base_count, "base");
    const std::vector<double> gtheta =
            read_grid(members, theta_grid, grid_count, "theta");
    const std::vector<double> gphi =
            read_grid(members, phi_grid, grid_count, "phi");
    std::vector<double> glens;
    if (has_lens) glens = read_grid(members, lens_grid, grid_count, "lens");

    // weights
    {
        const std::string spec = json_get(header, "weights").scalar;
        const std::size_t colon = spec.find(':');
        const std::string kind = spec.substr(0, colon);
        const int count = colon == std::string::npos
                ? n_rot : std::atoi(spec.c_str() + colon + 1);
        const std::vector<unsigned char>& wb = members.at("weights.bin");
        if (kind == "varint") {
            out.weights = decode_varint(wb, count);
        } else {
            if (wb.size() != static_cast<std::size_t>(count) * 4) {
                IMP_THROW("read_drot: weights.bin size disagrees", IOException);
            }
            out.weights.resize(count);
            for (int i = 0; i < count; ++i) {
                float v;
                std::memcpy(&v, &wb[4 * i], 4);
                out.weights[i] = v;
            }
        }
        if (static_cast<int>(out.weights.size()) != n_rot) {
            IMP_THROW("read_drot: weight count != n_rotamers", IOException);
        }
    }

    // Reconstruction, a block of conformers at a time with the rows inside.
    //
    // The obvious order -- one conformer, all its rows -- walks a column-major
    // grid with a stride of `n_rot` doubles, so every row of every conformer
    // is its own cache line: 5.7 kB apart on the shipped libraries, a miss per
    // row. Turning the loops inside out reads each row's values as the run of
    // neighbouring conformers they are stored as, and blocking keeps the
    // block's coordinates (64 conformers, ~127 kB) resident while its rows are
    // placed. The arithmetic is untouched and the order within one conformer
    // is untouched, so the numbers are bit-identical to the plain order.
    //
    // Row `r` may only be placed once its three parents are, and all three are
    // earlier rows or base atoms -- the same reason the plain order works --
    // so placing row `r` for every conformer in the block before row `r + 1`
    // is legal for exactly the same reason.
    out.coords.assign(static_cast<std::size_t>(n_rot) * n_atoms * 3, 0.0);
    // How many conformers share one pass over the rows. Measured against the
    // plain order (a block of 1 *is* the plain order) at 1, 4, 16, 64, 256
    // and 1024: the spread was under 10 % and inside the noise of a loaded
    // machine, so the cache argument for blocking did not survive contact
    // with the measurement. 16 is kept because it was at the fast end and
    // because writing straight into the output, which this order allows,
    // removes a per-conformer frame buffer and its copy.
    const int kBlock = 16;
    const std::size_t stride = static_cast<std::size_t>(n_atoms) * 3;
#pragma omp parallel for schedule(static)
    for (int k0 = 0; k0 < n_rot; k0 += kBlock) {
        const int k1 = std::min(k0 + kBlock, n_rot);
        for (int k = k0; k < k1; ++k) {
            double* f = &out.coords[static_cast<std::size_t>(k) * stride];
            std::memcpy(f, &tpl[0], stride * sizeof(double));
            for (int b = 0; b < n_base; ++b) {
                for (int c = 0; c < 3; ++c) {
                    const std::size_t idx = base_grid.col_major
                            ? (static_cast<std::size_t>(b * 3 + c) * n_rot + k)
                            : (static_cast<std::size_t>(k) * n_base * 3 +
                               b * 3 + c);
                    f[base_off[b] * 3 + c] = gbase[idx];
                }
            }
        }
        for (int r = 0; r < n_rows; ++r) {
            const int i = rows[4 * r], p1 = rows[4 * r + 1],
                      p2 = rows[4 * r + 2], p3 = rows[4 * r + 3];
            const double tl = tpl_lens[r];
            for (int k = k0; k < k1; ++k) {
                double* f = &out.coords[static_cast<std::size_t>(k) * stride];
                const std::size_t idx = theta_grid.col_major
                        ? (static_cast<std::size_t>(r) * n_rot + k)
                        : (static_cast<std::size_t>(k) * n_rows + r);
                const algebra::Vector3D placed = internal2cartesian(
                        algebra::Vector3D(f[3 * p3], f[3 * p3 + 1],
                                          f[3 * p3 + 2]),
                        algebra::Vector3D(f[3 * p2], f[3 * p2 + 1],
                                          f[3 * p2 + 2]),
                        algebra::Vector3D(f[3 * p1], f[3 * p1 + 1],
                                          f[3 * p1 + 2]),
                        has_lens ? glens[idx] : tl, gtheta[idx], gphi[idx]);
                f[3 * i] = placed[0];
                f[3 * i + 1] = placed[1];
                f[3 * i + 2] = placed[2];
            }
        }
    }
    return out;
}

std::vector<std::string> split_probe_rotamer_drot_locator(const std::string& locator) {
    std::vector<std::string> out;
    const std::size_t at = locator.rfind("::");
    if (at == std::string::npos) {
        out.push_back(locator);
        out.push_back(std::string());
    } else {
        out.push_back(locator.substr(0, at));
        out.push_back(locator.substr(at + 2));
    }
    return out;
}

ProbeRotamerLibrary read_probe_rotamer_drot(const std::string& path) {
    const std::vector<std::string> parts = split_probe_rotamer_drot_locator(path);
    return read_probe_rotamer_drot(parts[0], parts[1]);
}

std::vector<std::string> probe_rotamer_drot_catalog(const std::string& path) {
    std::vector<std::string> out;
    if (!pto::is_pto_file(path)) {
        out.push_back(std::string());     // the older envelope: one library
        return out;
    }
    pto::File pto;
    if (!pto.open(path, false)) {
        IMP_THROW("PTO: cannot open container " << path << ": "
                  << (pto.error().empty() ? "not an EBML document with DocType pto"
                                            : pto.error()), IOException);
    }
    const auto matches = pto.find_all("drot.catalog");
    if (!matches.empty()) {
        // The catalog is the one object read: `{"libraries":["a","b",...]}`.
        const std::vector<unsigned char> bytes = pto.read(matches.front());
        const Json listed = JsonParser(bytes).parse();
        const Json& libraries = json_get(listed, "libraries");
        for (std::size_t i = 0; i < libraries.array.size(); ++i) {
            out.push_back(libraries.array[i].scalar);
        }
        return out;
    }
    // No catalog: either one library, or a bundle written without one.
    const std::vector<pto::PtoObject>& objects = pto.objects();
    for (std::size_t i = 0; i < objects.size(); ++i) {
        const std::size_t slash = objects[i].name.find('/');
        if (slash == std::string::npos) continue;
        const std::string lib = objects[i].name.substr(0, slash);
        if (std::find(out.begin(), out.end(), lib) == out.end()) {
            out.push_back(lib);
        }
    }
    if (out.empty()) out.push_back(std::string());
    return out;
}

MfdbTags probe_rotamer_drot_provenance(const std::string& path, const std::string& library) {
    MfdbTags out;
    if (!pto::is_pto_file(path)) return out;   // pre-PTO: not stated
    pto::File pto;
    if (!pto.open(path, false)) {
        IMP_THROW("PTO: cannot open container " << path << ": "
                  << (pto.error().empty() ? "not an EBML document with DocType pto"
                                            : pto.error()), IOException);
    }
    const std::string want = library.empty() ? std::string("provenance.json")
                                             : library + "/provenance.json";
    const auto matches = pto.find_all(want);
    if (matches.empty()) return out;
    const std::vector<unsigned char> bytes = pto.read(matches.front());
    const Json listed = JsonParser(bytes).parse();
    for (std::size_t i = 0; i < listed.array.size(); ++i) {
        const Json& e = listed.array[i];
        out.push_back(MfdbTag(json_get(e, "item").scalar,
                              json_get(e, "value").scalar));
    }
    return out;
}

// ===========================================================================
// The `.drot` container: writing.
// ===========================================================================

namespace {

// -- layer 2: grid serialisation --------------------------------------------

//! Byte-plane shuffle: all first bytes, then all second bytes, ...
std::vector<unsigned char> byte_shuffle(const std::vector<unsigned char>& in,
                                        std::size_t item) {
    const std::size_t n = in.size() / item;
    std::vector<unsigned char> out(in.size());
    for (std::size_t b = 0; b < item; ++b) {
        for (std::size_t i = 0; i < n; ++i) {
            out[b * n + i] = in[i * item + b];
        }
    }
    return out;
}

//! `values[dof][frame]` (column-major) as float32, byte-shuffled.
std::vector<unsigned char> pack_f32(const std::vector<double>& values) {
    std::vector<unsigned char> raw(values.size() * 4);
    for (std::size_t i = 0; i < values.size(); ++i) {
        const float v = static_cast<float>(values[i]);
        std::memcpy(&raw[4 * i], &v, 4);
    }
    return byte_shuffle(raw, 4);
}

//! The same on the compact rung: int16 counts of `scale`, little-endian.
std::vector<unsigned char> pack_i16(const std::vector<double>& values,
                                    double scale, const char* what) {
    std::vector<unsigned char> raw(values.size() * 2);
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double q = values[i] / scale;
        if (!(q > -32768.5 && q < 32767.5)) {
            IMP_THROW("write_drot: " << what << " value " << values[i]
                      << " does not fit an int16 grid of " << scale
                      << " -- use the lossless encoding", ValueException);
        }
        const short s = static_cast<short>(q > 0 ? q + 0.5 : q - 0.5);
        raw[2 * i] = static_cast<unsigned char>(s & 0xFF);
        raw[2 * i + 1] = static_cast<unsigned char>((s >> 8) & 0xFF);
    }
    return raw;
}

//! Weights as varints when they are whole numbers (cluster populations are).
bool weights_are_counts(const std::vector<double>& w) {
    for (std::size_t i = 0; i < w.size(); ++i) {
        if (w[i] < 0.0 || w[i] > 4.0e9) return false;
        const double r = w[i] - static_cast<double>(
                static_cast<unsigned long long>(w[i] + 0.5));
        if (r > 1e-9 || r < -1e-9) return false;
    }
    return true;
}

std::vector<unsigned char> pack_weights(const std::vector<double>& w,
                                        bool as_counts) {
    std::vector<unsigned char> out;
    if (as_counts) {
        for (std::size_t i = 0; i < w.size(); ++i) {
            unsigned long long v = static_cast<unsigned long long>(w[i] + 0.5);
            do {
                unsigned char b = static_cast<unsigned char>(v & 0x7F);
                v >>= 7;
                if (v != 0) b |= 0x80;
                out.push_back(b);
            } while (v != 0);
        }
    } else {
        out.resize(w.size() * 4);
        for (std::size_t i = 0; i < w.size(); ++i) {
            const float v = static_cast<float>(w[i]);
            std::memcpy(&out[4 * i], &v, 4);
        }
    }
    return out;
}

//! Validate output compression before creating a destination file.
void validate_drot_codec(const ProbeRotamerDrotEncoding& encoding) {
    if ((encoding.codec != "zstd" && encoding.codec != "brotli") ||
        !pto::can_compress(encoding.codec)) {
        IMP_THROW("write_drot: unavailable output codec '" << encoding.codec
                  << "'", ValueException);
    }
    if (encoding.codec == "brotli" &&
        (encoding.quality < 0 || encoding.quality > 11)) {
        IMP_THROW("write_drot: brotli quality must be between 0 and 11",
                  ValueException);
    }
}

//! Escape a string for a JSON scalar.
/*! Provenance values are free text -- an algorithm name carries parentheses and
    a settings blob is itself JSON -- so the quotes and backslashes in them must
    not end the string they are being written into. */
std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (std::size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    char buf[7];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    out += buf;
                } else {
                    out += static_cast<char>(c);
                }
        }
    }
    return out;
}

void add_object(pto::File& pto, const std::string& name,
                const std::string& kind, const std::string& encoding_name,
                const std::vector<unsigned char>& payload,
                const ProbeRotamerDrotEncoding& encoding) {
    if (pto.add_coded(kind, encoding_name, encoding.codec, name,
                      payload.data(), payload.size(), encoding.quality) == 0) {
        IMP_THROW("PTO: writing " << pto.filename() << " failed: "
                  << pto.error(), IOException);
    }
}

std::vector<unsigned char> to_bytes(const std::string& s) {
    return std::vector<unsigned char>(s.begin(), s.end());
}

std::string fixed(double v, int digits) {
    std::ostringstream o;
    o.precision(digits);
    o << std::fixed << v;
    return o.str();
}

}  // namespace

void write_probe_rotamer_drot_bundle(const std::vector<std::string>& sources,
                       const std::vector<std::string>& names,
                       const std::string& path) {
    const ProbeRotamerDrotEncoding defaults;
    validate_drot_codec(defaults);
    if (sources.size() != names.size()) {
        IMP_THROW("write_drot_bundle: " << sources.size() << " sources but "
                  << names.size() << " names", ValueException);
    }
    if (sources.empty()) {
        IMP_THROW("write_drot_bundle: nothing to bundle", ValueException);
    }
    for (std::size_t i = 0; i < names.size(); ++i) {
        if (names[i].empty() || names[i].find('/') != std::string::npos) {
            IMP_THROW("write_drot_bundle: '" << names[i]
                      << "' cannot be a library name (empty, or holds a '/')",
                      ValueException);
        }
    }

    pto::File pto;
    if (!pto.create(path, "",
            std::string(pto::kDefaultBanner) + "\nThis container was written by IMP.bff.\nhttps://github.com/tpeulen/IMP.bff\n")) {
        IMP_THROW("PTO: cannot open " << path << " for writing: "
                  << pto.error(), IOException);
    }
    pto.set_writing_app("IMP.bff");

    // The catalog first, so a reader that wants only the listing stops after
    // one object. It is the only thing this function composes; everything
    // else is copied.
    std::ostringstream cj;
    cj << "{\"format\":\"drot.bundle\",\"version\":10"
       << ",\"producer\":\"IMP.bff write_drot_bundle\""
       << ",\"n_libraries\":" << names.size() << ",\"libraries\":[";
    for (std::size_t i = 0; i < names.size(); ++i) {
        if (i) cj << ',';
        cj << '"' << names[i] << '"';
    }
    cj << "]}";
    add_object(pto, "drot.catalog", "drot.catalog", "json",
               to_bytes(cj.str()), defaults);

    for (std::size_t i = 0; i < sources.size(); ++i) {
        pto::File in;
        if (!in.open(sources[i], false)) {
            IMP_THROW("PTO: cannot open container " << sources[i] << ": "
                      << (in.error().empty() ? "not an EBML document with DocType pto"
                                                : in.error()), IOException);
        }
        const std::vector<pto::PtoObject>& objects = in.objects();
        if (objects.empty()) {
            IMP_THROW("write_drot_bundle: " << sources[i] << " holds nothing",
                      IOException);
        }
        for (std::size_t j = 0; j < objects.size(); ++j) {
            if (objects[j].name == "drot.catalog") continue;  // never nested
            // Payload bytes cross unread and unaltered: they are already the
            // encoding their `PtoEncoding` names, and re-compressing them
            // would only cost time and change nothing. read_stored, not
            // data(): data() decodes by the encoding, and decoded bytes
            // labelled `+brotli` are a stream no reader can decode.
            const std::vector<unsigned char> payload =
                    in.read_stored(objects[j].uid);
            if (pto.add(objects[j].kind, objects[j].encoding, names[i] + "/" + objects[j].name, payload.empty() ? NULL : &payload[0], payload.size()) == 0) {
                IMP_THROW("PTO: writing " << pto.filename() << " failed: "
                          << pto.error(), IOException);
            }
        }
    }
    const bool committed = pto.commit();
    const std::string commit_error = pto.error();
    pto.close();
    if (!committed) IMP_THROW("PTO: writing " << path << " failed: "
                             << commit_error, IOException);
}

namespace {
void write_drot_impl(const std::string& path,
                double* rotamer_coords, int n_rotamer_coords,
                const std::vector<std::string>& atom_names,
                const std::vector<std::string>& elements,
                const std::vector<std::string>& resnames,
                double* rotamer_weights, int n_rotamer_weights,
                const ProbeRotamerDrotEncoding& encoding, const MfdbTags& provenance) {
    validate_drot_codec(encoding);
    const int n_atoms = static_cast<int>(atom_names.size());
    if (n_atoms < 4) {
        IMP_THROW("write_drot: need at least four atoms, got " << n_atoms,
                  ValueException);
    }
    if (static_cast<int>(elements.size()) != n_atoms) {
        IMP_THROW("write_drot: " << elements.size() << " elements for "
                  << n_atoms << " atoms", ValueException);
    }
    if (!resnames.empty() && static_cast<int>(resnames.size()) != n_atoms) {
        IMP_THROW("write_drot: " << resnames.size() << " residue names for "
                  << n_atoms << " atoms", ValueException);
    }
    if (rotamer_coords == NULL || n_rotamer_coords <= 0 ||
        n_rotamer_coords % (3 * n_atoms) != 0) {
        IMP_THROW("write_drot: " << n_rotamer_coords
                  << " coordinates is not a whole number of frames of "
                  << n_atoms << " atoms", ValueException);
    }
    const int n_frames = n_rotamer_coords / (3 * n_atoms);
    std::vector<double> weights;
    if (rotamer_weights != NULL && n_rotamer_weights > 0) {
        if (n_rotamer_weights != n_frames) {
            IMP_THROW("write_drot: " << n_rotamer_weights << " weights for "
                      << n_frames << " conformers", ValueException);
        }
        weights.assign(rotamer_weights, rotamer_weights + n_rotamer_weights);
    } else {
        weights.assign(n_frames, 1.0);
    }

    // -- the Z-matrix, from the first conformer --------------------------
    std::vector<algebra::Vector3D> first(n_atoms);
    for (int a = 0; a < n_atoms; ++a) {
        first[a] = algebra::Vector3D(rotamer_coords[3 * a],
                                     rotamer_coords[3 * a + 1],
                                     rotamer_coords[3 * a + 2]);
    }
    int root = 0;
    for (int a = 0; a < n_atoms; ++a) {
        if (atom_names[a] == "CA") { root = a; break; }
    }
    ZMatrix z;
    try {
        z.set_template(first, elements, root);
    } catch (const std::invalid_argument& e) {
        IMP_THROW("write_drot: " << e.what(), ValueException);
    }
    const std::vector<int> rows = z.get_rows();
    const std::vector<int> base = z.get_base();
    const int n_rows = static_cast<int>(rows.size() / 4);
    const int n_base = static_cast<int>(base.size());
    if (n_rows < 1) {
        IMP_THROW("write_drot: the molecule has no Z-matrix rows (only "
                  << n_base << " base atoms)", ValueException);
    }

    // -- measure every conformer, column-major ---------------------------
    std::vector<double> gbase(static_cast<std::size_t>(n_base) * 3 * n_frames);
    std::vector<double> glen(static_cast<std::size_t>(n_rows) * n_frames);
    std::vector<double> gang(glen.size()), gphi(glen.size());
    {
        std::vector<algebra::Vector3D> frame(n_atoms);
        std::vector<double> lens, angs, chi;
        for (int k = 0; k < n_frames; ++k) {
            const double* c = rotamer_coords +
                    static_cast<std::size_t>(k) * n_atoms * 3;
            for (int a = 0; a < n_atoms; ++a) {
                frame[a] = algebra::Vector3D(c[3 * a], c[3 * a + 1],
                                             c[3 * a + 2]);
            }
            z.frame_internals(frame, lens, angs, chi);
            for (int b = 0; b < n_base; ++b) {
                for (int j = 0; j < 3; ++j) {
                    gbase[static_cast<std::size_t>(b * 3 + j) * n_frames + k] =
                            frame[base[b]][j];
                }
            }
            for (int r = 0; r < n_rows; ++r) {
                const std::size_t idx =
                        static_cast<std::size_t>(r) * n_frames + k;
                glen[idx] = lens[r];
                gang[idx] = angs[r];
                gphi[idx] = chi[r];
            }
        }
    }

    // -- members ---------------------------------------------------------
    const bool f32 = encoding.lossless;
    const std::string suffix = f32 ? ".f32" : ".i16";
    std::vector<unsigned char> mbase, mlen, mang, mphi;
    if (f32) {
        mbase = pack_f32(gbase);
        mlen = pack_f32(glen);
        mang = pack_f32(gang);
        mphi = pack_f32(gphi);
    } else {
        mbase = pack_i16(gbase, encoding.grid_a, "base coordinate");
        mlen = pack_i16(glen, encoding.grid_a, "bond length");
        mang = pack_i16(gang, encoding.grid_deg, "bond angle");
        mphi = pack_i16(gphi, encoding.grid_deg, "dihedral");
    }

    std::ostringstream cif;
    cif << "data_drot_template\nloop_\n_atom_site.label_atom_id\n";
    if (!resnames.empty()) cif << "_atom_site.label_comp_id\n";
    cif << "_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n"
        << "_atom_site.Cartn_z\n";
    for (int a = 0; a < n_atoms; ++a) {
        cif << atom_names[a] << ' ';
        if (!resnames.empty()) cif << resnames[a] << ' ';
        cif << elements[a] << ' ' << fixed(first[a][0], 3) << ' '
            << fixed(first[a][1], 3) << ' ' << fixed(first[a][2], 3) << '\n';
    }

    std::ostringstream rj;
    rj << '[';
    for (std::size_t i = 0; i < rows.size(); ++i) {
        if (i) rj << ',';
        rj << rows[i];
    }
    rj << ']';

    const bool counts = weights_are_counts(weights);
    const std::vector<unsigned char> mweights = pack_weights(weights, counts);

    std::ostringstream hj;
    hj << "{\"format\":\"drot\",\"version\":10,\"container\":\"pto\""
       << ",\"producer\":\"IMP.bff write_drot\""
       << ",\"n_atoms\":" << n_atoms
       << ",\"n_rotamers\":" << n_frames
       << ",\"n_rows\":" << n_rows
       << ",\"n_base\":" << n_base
       << ",\"base_offset\":[";
    for (int b = 0; b < n_base; ++b) {
        if (b) hj << ',';
        hj << base[b];
    }
    hj << "],\"grids\":{";
    const char* keys[4] = {"base", "lens", "theta", "phi"};
    for (int g = 0; g < 4; ++g) {
        if (g) hj << ',';
        hj << '"' << keys[g] << "\":{\"file\":\"" << keys[g] << suffix
           << "\",\"dtype\":\"" << (f32 ? "f32" : "i16")
           << "\",\"order\":\"col\",\"shuffle\":" << (f32 ? 1 : 0);
        if (!f32) {
            hj << ",\"scale\":"
               << (g >= 2 ? encoding.grid_deg : encoding.grid_a);
        }
        hj << '}';
    }
    hj << "},\"weights\":\"" << (counts ? "varint:" : "f32:") << n_frames
       << "\"}";

    // -- the PTO envelope --------------------------------------------------
    // independent codec frames preserve selective reads
    pto::File pto;
    if (!pto.create(path, "",
            std::string(pto::kDefaultBanner) + "\nThis container was written by IMP.bff.\nhttps://github.com/tpeulen/IMP.bff\n")) {
        IMP_THROW("PTO: cannot open " << path << " for writing: "
                  << pto.error(), IOException);
    }
    pto.set_writing_app("IMP.bff");
    const std::string enc = std::string(f32 ? "f32.col" : "i16.col");
    add_object(pto, "drot.json", "drot.header", "json",
               to_bytes(hj.str()), encoding);
    add_object(pto, "template.cif", "drot.template", "mmcif",
               to_bytes(cif.str()), encoding);
    add_object(pto, "rows.json", "drot.rows", "json",
               to_bytes(rj.str()), encoding);
    add_object(pto, std::string("base") + suffix, "drot.grid", enc, mbase,
               encoding);
    add_object(pto, std::string("lens") + suffix, "drot.grid", enc, mlen,
               encoding);
    add_object(pto, std::string("theta") + suffix, "drot.grid", enc, mang,
               encoding);
    add_object(pto, std::string("phi") + suffix, "drot.grid", enc, mphi,
               encoding);
    add_object(pto, "weights.bin", "drot.weights",
               std::string(counts ? "varint" : "f32"), mweights,
               encoding);
    if (!provenance.empty()) {
        // mmCIF item/value pairs, as JSON so a reader needs no CIF parser to
        // answer "how was this made". Written last: it is the one object a
        // reader may want without touching a single coordinate.
        std::ostringstream pj;
        pj << '[';
        for (std::size_t i = 0; i < provenance.size(); ++i) {
            if (i) pj << ',';
            pj << "{\"item\":\"" << json_escape(provenance[i].item)
               << "\",\"value\":\"" << json_escape(provenance[i].value) << "\"}";
        }
        pj << ']';
        add_object(pto, "provenance.json", "drot.provenance", "json",
                   to_bytes(pj.str()), encoding);
    }
    const bool committed = pto.commit();
    const std::string commit_error = pto.error();
    pto.close();
    if (!committed) IMP_THROW("PTO: writing " << path << " failed: "
                             << commit_error, IOException);
}
}  // namespace

void write_probe_rotamer_drot(const std::string& path,
                double* rotamer_coords, int n_rotamer_coords,
                const std::vector<std::string>& atom_names,
                const std::vector<std::string>& elements,
                const std::vector<std::string>& resnames,
                double* rotamer_weights, int n_rotamer_weights,
                const ProbeRotamerDrotEncoding& encoding) {
    write_drot_impl(path, rotamer_coords, n_rotamer_coords, atom_names,
                    elements, resnames, rotamer_weights, n_rotamer_weights,
                    encoding, MfdbTags());
}

void write_probe_rotamer_drot_with_provenance(const std::string& path,
                double* rotamer_coords, int n_rotamer_coords,
                const std::vector<std::string>& atom_names,
                const std::vector<std::string>& elements,
                const std::vector<std::string>& resnames,
                double* rotamer_weights, int n_rotamer_weights,
                const ProbeRotamerDrotEncoding& encoding, const MfdbTags& provenance) {
    write_drot_impl(path, rotamer_coords, n_rotamer_coords, atom_names,
                    elements, resnames, rotamer_weights, n_rotamer_weights,
                    encoding, provenance);
}

// --------------------------------------------------------------------------
// Rotamer library IO
// --------------------------------------------------------------------------

namespace {

std::string base_path(const std::string& p) {
    if (p.size() > 4 && (p.substr(p.size() - 4) == ".rmf" || p.substr(p.size() - 4) == ".npy"))
        return p.substr(0, p.size() - 4);
    return p;
}

// Minimal .npy reader for float64 arrays. The numpy .npy format is:
// magic \x93NUMPY, version, header length, header dict, then raw data.
std::vector<double> read_npy(const std::string& path, int& n0, int& n1, int& n2) {
    std::ifstream f(path, std::ios::binary);
    if (!f) IMP_THROW("cannot open " << path, IOException);
    char magic[6];
    f.read(magic, 6);
    if (std::string(magic, 6) != "\x93NUMPY") IMP_THROW("not a .npy file: " << path, IOException);
    char version[2];
    f.read(version, 2);
    int header_len = version[0] == 1 ? 0 : 0;
    if (version[0] == 1) { uint16_t h; f.read(reinterpret_cast<char*>(&h), 2); header_len = h; }
    else if (version[0] == 2) { uint32_t h; f.read(reinterpret_cast<char*>(&h), 4); header_len = h; }
    std::string header(header_len, '\0');
    f.read(&header[0], header_len);
    // Parse shape from header: "'shape': (N, M, 3),"
    // This is a minimal parser; numpy headers are Python dict literals.
    n0 = n1 = n2 = 0;
    auto shape_pos = header.find("'shape'");
    if (shape_pos == std::string::npos) shape_pos = header.find("\"shape\"");
    if (shape_pos != std::string::npos) {
        auto paren = header.find('(', shape_pos);
        if (paren != std::string::npos) {
            std::string s = header.substr(paren + 1);
            // parse comma-separated ints
            std::istringstream ss(s);
            char c;
            std::vector<int> dims;
            int d;
            while (ss >> d) {
                dims.push_back(d);
                ss >> c;
                if (c != ',') break;
            }
            if (dims.size() >= 1) n0 = dims[0];
            if (dims.size() >= 2) n1 = dims[1];
            if (dims.size() >= 3) n2 = dims[2];
        }
    }
    // Read data
    size_t count = 1;
    if (n0 > 0) count *= n0;
    if (n1 > 0) count *= n1;
    if (n2 > 0) count *= n2;
    std::vector<double> data(count);
    f.read(reinterpret_cast<char*>(data.data()), count * sizeof(double));
    return data;
}


}  // namespace

ProbeRotamerLibrary read_probe_rotamer_library(const std::string& path) {
    std::string base = base_path(path);
    ProbeRotamerLibrary lib;
    int n0, n1, n2;
    lib.coords = read_npy(base + "_coords.npy", n0, n1, n2);
    lib.n_rotamers = n0;
    lib.n_atoms = n1;

    // weights
    std::string wpath = base + "_weights.txt";
    std::ifstream wf(wpath);
    if (wf) {
        std::string line;
        while (std::getline(wf, line)) {
            if (!line.empty()) lib.weights.push_back(std::stod(line));
        }
    } else {
        lib.weights.assign(lib.n_rotamers, 1.0 / lib.n_rotamers);
    }

    // atom names
    std::string apath = base + "_atoms.txt";
    std::ifstream af(apath);
    if (af) {
        std::string line;
        while (std::getline(af, line)) {
            if (!line.empty()) lib.atom_names.push_back(line);
        }
    } else {
        for (int i = 0; i < lib.n_atoms; i++)
            lib.atom_names.push_back("AT" + std::to_string(i));
    }

    lib.path = path;
    return lib;
}


IMPBFF_END_NAMESPACE
