/**
 * \file DrotReader.cpp
 * \brief The `.drot` reader: PTO, brotli, minimal JSON, exact reconstruction.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 * Layering, bottom to top, each deliberately small:
 *
 * 1. brotli: the one-shot `BrotliDecoderDecompress` over the vendored
 *    decoder (`src/brotli/`, google/brotli MIT). The decoded size is
 *    unknown up front, so the buffer grows on NEEDS_MORE_OUTPUT.
 * 2. tar: 512-byte headers (name, octal size, typeflag), data padded to
 *    512 -- the ustar subset our writer emits.
 * 3. JSON: a recursive-descent parser for the grammar the writer uses
 *    (objects, arrays, strings, numbers). No third-party dependency.
 * 4. Reconstruction: the embedded `rows.json` IS the Z-matrix tree, so no
 *    tree is rebuilt or guessed -- each row's bond length/angle is measured
 *    from `template.cif` and the conformer is placed by the shared
 *    `internal2cartesian`. The row table is validated (indices in range,
 *    references resolved before use) which is what the cross-check against
 *    a rebuilt tree would catch, without depending on BFS tie-breaking.
 */

#include <IMP/bff/DrotReader.h>
#include <IMP/bff/Pto.h>

#include <IMP/bff/ZMatrix.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/exception.h>

#include "brotli/include/brotli/decode.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <map>
#include <memory>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

namespace {

typedef std::map<std::string, std::string> JsonFlat;

// -- layer 1: brotli --------------------------------------------------------

//! Whole-file decompression, streamed: the decoded size is not known up front.
/*! **Not the one-shot API.** `BrotliDecoderDecompress` reports a too-small
    output buffer as `BROTLI_DECODER_RESULT_ERROR`, not as
    `NEEDS_MORE_OUTPUT` -- the growth branch below it never fired -- so a
    reader that guesses the size from the compressed length, as this one did
    with `size * 8`, fails outright on anything that compresses better than
    eight times. Every shipped `.drot` is under that ratio, which is why it
    worked; a library of near-identical conformers is not, and would have been
    the first file to be reported as corrupt while being perfectly good.

    The streaming decoder is told how much room it has and asks for more when
    it runs out, so there is nothing to guess and nothing to re-decode. */
std::vector<unsigned char> brotli_decompress(const unsigned char* data,
                                             std::size_t size) {
    BrotliDecoderState* state =
            BrotliDecoderCreateInstance(NULL, NULL, NULL);
    if (state == NULL) {
        IMP_THROW("read_drot: cannot start the brotli decoder", IOException);
    }

    std::vector<unsigned char> out;
    std::size_t chunk = std::max<std::size_t>(size * 4, 1u << 16);
    std::size_t available_in = size;
    const unsigned char* next_in = data;
    BrotliDecoderResult r = BROTLI_DECODER_RESULT_NEEDS_MORE_OUTPUT;

    while (r == BROTLI_DECODER_RESULT_NEEDS_MORE_OUTPUT) {
        if (out.size() > (std::size_t)1 << 29) {   // 512 MB sanity cap
            BrotliDecoderDestroyInstance(state);
            IMP_THROW("read_drot: decompressed payload exceeds 512 MB",
                      IOException);
        }
        const std::size_t written = out.size();
        out.resize(written + chunk);
        std::size_t available_out = chunk;
        unsigned char* next_out = &out[written];
        r = BrotliDecoderDecompressStream(state, &available_in, &next_in,
                                          &available_out, &next_out, NULL);
        out.resize(out.size() - available_out);
        chunk *= 2;
    }
    BrotliDecoderDestroyInstance(state);
    if (r != BROTLI_DECODER_RESULT_SUCCESS) {
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

RotamerLibrary read_drot(const std::string& path, const std::string& library) {
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
    if (PtoReader::looks_like_pto(path)) {
        const PtoReader pto(path);
        const std::vector<PtoObject>& objects = pto.objects();
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
            const std::vector<unsigned char> bytes = pto.data(objects[i]);
            const std::string& enc = objects[i].encoding;
            const bool packed = enc.size() >= 7 &&
                    enc.compare(enc.size() - 7, 7, "+brotli") == 0;
            members[name.substr(prefix.size())] = packed
                    ? brotli_decompress(bytes.empty() ? NULL : &bytes[0],
                                        bytes.size())
                    : bytes;
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
    RotamerLibrary out;
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

std::vector<std::string> split_drot_locator(const std::string& locator) {
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

RotamerLibrary read_drot(const std::string& path) {
    const std::vector<std::string> parts = split_drot_locator(path);
    return read_drot(parts[0], parts[1]);
}

std::vector<std::string> drot_catalog(const std::string& path) {
    std::vector<std::string> out;
    if (!PtoReader::looks_like_pto(path)) {
        out.push_back(std::string());     // the older envelope: one library
        return out;
    }
    const PtoReader pto(path);
    const int at = pto.find("drot.catalog");
    if (at >= 0) {
        // The catalog is the one object read: `{"libraries":["a","b",...]}`.
        const std::vector<unsigned char> packed = pto.data(pto.objects()[at]);
        const std::vector<unsigned char> bytes =
                brotli_decompress(packed.empty() ? NULL : &packed[0],
                                  packed.size());
        const Json listed = JsonParser(bytes).parse();
        const Json& libraries = json_get(listed, "libraries");
        for (std::size_t i = 0; i < libraries.array.size(); ++i) {
            out.push_back(libraries.array[i].scalar);
        }
        return out;
    }
    // No catalog: either one library, or a bundle written without one.
    const std::vector<PtoObject>& objects = pto.objects();
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

IMPBFF_END_NAMESPACE