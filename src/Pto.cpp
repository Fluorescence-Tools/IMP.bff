

/**
 * \file Pto.cpp
 * \brief The PTO envelope: EBML variable-width integers and nothing else.
 *
 * EBML is a handful of variable-width integers and twenty element IDs, and
 * the whole argument for the format is that a reader needs exactly that and
 * never a framework -- so this file is deliberately small. The layout it
 * writes is the one in `../tttrlib/okf/specs/pto-binary-decoding.md`:
 *
 *     EBML header (DocType "pto")
 *     Segment
 *       SeekHead      -- reserved at the head, patched at close
 *       Attachments
 *         AttachedFile  FileUID FileName PtoKind PtoEncoding FileData
 *         ...
 *
 * Two rules ride along, both paid for in tttrlib and both cheap here: payload
 * bytes land on an 8-byte boundary (padded with `Void`) so a mapped file
 * hands out slices rather than copies, and any Data Size that may be
 * rewritten in place is written at a fixed width, so patching it cannot shift
 * the bytes that follow.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from Pto.cpp --------

#include <IMP/bff/Pto.h>

#include <IMP/bff/Base.h>

#include <cstdio>
#include <cstring>

IMPBFF_BEGIN_NAMESPACE

namespace {

// -- element IDs, as the specification prints them ---------------------------
const unsigned int kEbmlHeader = 0x1A45DFA3u;
const unsigned int kEbmlVersion = 0x4286u;
const unsigned int kEbmlReadVersion = 0x42F7u;
const unsigned int kEbmlMaxIdLength = 0x42F2u;
const unsigned int kEbmlMaxSizeLength = 0x42F3u;
const unsigned int kDocType = 0x4282u;
const unsigned int kDocTypeVersion = 0x4287u;
const unsigned int kDocTypeReadVersion = 0x4285u;
const unsigned int kSegment = 0x18538067u;
const unsigned int kSeekHead = 0x114D9B74u;
const unsigned int kSeek = 0x4DBBu;
const unsigned int kSeekId = 0x53ABu;
const unsigned int kSeekPosition = 0x53ACu;
const unsigned int kVoid = 0xECu;
const unsigned int kAttachments = 0x1941A469u;
const unsigned int kAttachedFile = 0x61A7u;
const unsigned int kFileUid = 0x46AEu;
const unsigned int kFileName = 0x466Eu;
const unsigned int kFileData = 0x465Cu;
const unsigned int kPtoKind = 0x1E54F001u;
const unsigned int kPtoEncoding = 0x1E54F002u;

//! Payload alignment: 8 keeps float64 columns and 16-byte rows natural.
const std::size_t kAlignment = 8;
//! The width every rewritten Data Size is written at.
const std::size_t kFixedSizeWidth = 8;
//! One Seek entry, fixed by construction so the table patches in place.
const std::size_t kSeekEntryBytes = 19;

typedef std::vector<unsigned char> Bytes;

void put_id(Bytes& out, unsigned int id) {
    // The marker bits are part of the ID -- `0x61A7` is two octets, and
    // `0x1A45DFA3` is four -- which is why the spec prints IDs as numbers.
    int width = 1;
    while (width < 4 && (id >> (8 * width)) != 0) ++width;
    for (int i = width - 1; i >= 0; --i) {
        out.push_back(static_cast<unsigned char>((id >> (8 * i)) & 0xFF));
    }
}

//! A Data Size VINT: `width` octets, or the narrowest that fits.
void put_size(Bytes& out, unsigned long long value, std::size_t width = 0) {
    if (width == 0) {
        width = 1;
        while (width <= 8 && value > ((1ULL << (7 * width)) - 2)) ++width;
    }
    if (width > 8 || value > ((1ULL << (7 * width)) - 2)) {
        IMP_THROW("PTO: " << value << " does not fit a " << width
                          << "-octet Data Size", IOException);
    }
    const unsigned long long raw = (1ULL << (7 * width)) | value;
    for (int i = static_cast<int>(width) - 1; i >= 0; --i) {
        out.push_back(static_cast<unsigned char>((raw >> (8 * i)) & 0xFF));
    }
}

void put_uint(Bytes& out, unsigned int id, unsigned long long value,
              std::size_t octets = 0) {
    if (octets == 0) {
        octets = 1;
        while (octets < 8 && (value >> (8 * octets)) != 0) ++octets;
    }
    put_id(out, id);
    put_size(out, octets);
    for (int i = static_cast<int>(octets) - 1; i >= 0; --i) {
        out.push_back(static_cast<unsigned char>((value >> (8 * i)) & 0xFF));
    }
}

void put_string(Bytes& out, unsigned int id, const std::string& text) {
    put_id(out, id);
    put_size(out, text.size());
    out.insert(out.end(), text.begin(), text.end());
}

//! A `Void` element occupying exactly `total` octets (0, or at least 2).
Bytes void_bytes(std::size_t total) {
    Bytes out;
    if (total == 0) return out;
    if (total < 2) {
        IMP_THROW("PTO: a Void element needs at least two octets", IOException);
    }
    for (std::size_t width = 1; width <= 8; ++width) {
        if (total < 1 + width) continue;
        const std::size_t payload = total - 1 - width;
        if (payload <= ((1ULL << (7 * width)) - 2)) {
            put_id(out, kVoid);
            put_size(out, payload, width);
            out.resize(out.size() + payload, 0);
            return out;
        }
    }
    IMP_THROW("PTO: cannot write a " << total << "-octet Void", IOException);
}

// -- reading ----------------------------------------------------------------

//! Decode an element ID at `at`; returns its width, 0 when there is none.
std::size_t read_id(const unsigned char* buf, std::size_t left,
                    unsigned int* out) {
    if (left == 0 || buf[0] == 0) return 0;
    std::size_t width = 1;
    while (width <= 4 && !(buf[0] & (0x80 >> (width - 1)))) ++width;
    if (width > 4 || width > left) return 0;
    unsigned int value = 0;
    for (std::size_t i = 0; i < width; ++i) {
        value = (value << 8) | buf[i];
    }
    *out = value;
    return width;
}

//! Decode a Data Size at `at`. An unknown size is legal EBML for a live
//! stream and the one thing PTO promises never to write, so it is refused
//! here rather than guessed at.
std::size_t read_size(const unsigned char* buf, std::size_t left,
                      unsigned long long* out) {
    if (left == 0 || buf[0] == 0) return 0;
    std::size_t width = 1;
    while (width <= 8 && !(buf[0] & (0x80 >> (width - 1)))) ++width;
    if (width > 8 || width > left) return 0;
    unsigned long long value = buf[0] & (0xFFu >> width);
    for (std::size_t i = 1; i < width; ++i) value = (value << 8) | buf[i];
    if (value == ((1ULL << (7 * width)) - 1)) return 0;   // unknown size
    *out = value;
    return width;
}

std::string as_string(const unsigned char* at, std::size_t n) {
    return std::string(reinterpret_cast<const char*>(at), n);
}

unsigned long long as_uint(const unsigned char* at, std::size_t n) {
    unsigned long long value = 0;
    for (std::size_t i = 0; i < n; ++i) value = (value << 8) | at[i];
    return value;
}

}  // namespace

// --------------------------------------------------------------------------
// the writer
// --------------------------------------------------------------------------

struct PtoWriter::Impl {
    std::FILE* file;
    std::string path;
    long segment_size_at, segment_data_at;
    long seek_head_data_at;
    long attachments_size_at, attachments_data_at;
    std::size_t seek_head_bytes;
    std::vector<long> seeks;
    bool closed;

    Impl() : file(NULL), segment_size_at(0), segment_data_at(0),
             seek_head_data_at(0), attachments_size_at(0),
             attachments_data_at(0), seek_head_bytes(0), closed(false) {}

    void write(const Bytes& b) {
        if (b.empty()) return;
        if (std::fwrite(&b[0], 1, b.size(), file) != b.size()) {
            IMP_THROW("PTO: writing " << path << " failed", IOException);
        }
    }
};

PtoWriter::PtoWriter(const std::string& path, int doctype_version)
        : impl_(new Impl) {
    impl_->path = path;
    impl_->file = std::fopen(path.c_str(), "wb+");
    if (impl_->file == NULL) {
        delete impl_;
        impl_ = NULL;
        IMP_THROW("PTO: cannot open " << path << " for writing", IOException);
    }

    Bytes header;
    put_uint(header, kEbmlVersion, 1);
    put_uint(header, kEbmlReadVersion, 1);
    put_uint(header, kEbmlMaxIdLength, 4);
    put_uint(header, kEbmlMaxSizeLength, 8);
    put_string(header, kDocType, "pto");
    put_uint(header, kDocTypeVersion, static_cast<unsigned long long>(
                                              doctype_version));
    put_uint(header, kDocTypeReadVersion, 1);

    Bytes out;
    put_id(out, kEbmlHeader);
    put_size(out, header.size());
    out.insert(out.end(), header.begin(), header.end());
    impl_->write(out);

    out.clear();
    put_id(out, kSegment);
    impl_->write(out);
    impl_->segment_size_at = std::ftell(impl_->file);
    out.clear();
    put_size(out, 0, kFixedSizeWidth);
    impl_->write(out);
    impl_->segment_data_at = std::ftell(impl_->file);

    // The seek table is reserved here, at the head of the Segment, and filled
    // at close: a reader then reaches the directory objects without walking
    // every attachment. Reserved as Void, so an interrupted write still
    // leaves a document that walks.
    impl_->seek_head_bytes = kSeekEntryBytes * 8;
    out.clear();
    put_id(out, kSeekHead);
    put_size(out, impl_->seek_head_bytes, kFixedSizeWidth);
    impl_->write(out);
    impl_->seek_head_data_at = std::ftell(impl_->file);
    impl_->write(void_bytes(impl_->seek_head_bytes));

    out.clear();
    put_id(out, kAttachments);
    impl_->write(out);
    impl_->attachments_size_at = std::ftell(impl_->file);
    out.clear();
    put_size(out, 0, kFixedSizeWidth);
    impl_->write(out);
    impl_->attachments_data_at = std::ftell(impl_->file);
}

PtoWriter::~PtoWriter() {
    try {
        close();
    } catch (...) {
        // a destructor is the wrong place to raise; close() reports properly
    }
    delete impl_;
    impl_ = NULL;
}

PtoObject PtoWriter::add(const std::string& name, const std::string& kind,
                         const std::string& encoding, const void* data,
                         std::size_t size) {
    if (impl_ == NULL || impl_->closed) {
        IMP_THROW("PTO: the writer is closed", IOException);
    }
    const unsigned long long uid = objects_.size() + 1;

    Bytes children;
    put_uint(children, kFileUid, uid, 8);
    put_string(children, kFileName, name);
    put_string(children, kPtoKind, kind);
    put_string(children, kPtoEncoding, encoding);

    Bytes data_head;
    put_id(data_head, kFileData);
    put_size(data_head, size, kFixedSizeWidth);

    Bytes head;
    put_id(head, kAttachedFile);
    put_size(head, children.size() + data_head.size() + size, kFixedSizeWidth);
    head.insert(head.end(), children.begin(), children.end());
    head.insert(head.end(), data_head.begin(), data_head.end());

    long at = std::ftell(impl_->file);
    // Pad so the payload -- not the header -- lands on the boundary. A Void
    // cannot be one octet long, so a shortfall of one borrows a whole one.
    std::size_t shortfall =
            (kAlignment - ((at + head.size()) % kAlignment)) % kAlignment;
    if (shortfall == 1) shortfall += kAlignment;
    if (shortfall != 0) {
        impl_->write(void_bytes(shortfall));
        at += static_cast<long>(shortfall);
    }
    impl_->write(head);

    PtoObject object;
    object.uid = uid;
    object.name = name;
    object.kind = kind;
    object.encoding = encoding;
    object.offset = static_cast<std::size_t>(std::ftell(impl_->file));
    object.size = size;
    if (size != 0) {
        if (std::fwrite(data, 1, size, impl_->file) != size) {
            IMP_THROW("PTO: writing " << impl_->path << " failed", IOException);
        }
    }
    // Every object of a `.drot.pto` is a directory object -- there are eight
    // of them, not sixty thousand -- so all of them are seekable.
    if (impl_->seeks.size() * kSeekEntryBytes < impl_->seek_head_bytes) {
        impl_->seeks.push_back(at - impl_->segment_data_at);
    }
    objects_.push_back(object);
    return object;
}

void PtoWriter::close() {
    if (impl_ == NULL || impl_->closed) return;
    const long end = std::ftell(impl_->file);

    Bytes size_bytes;
    put_size(size_bytes, end - impl_->attachments_data_at, kFixedSizeWidth);
    std::fseek(impl_->file, impl_->attachments_size_at, SEEK_SET);
    impl_->write(size_bytes);

    size_bytes.clear();
    put_size(size_bytes, end - impl_->segment_data_at, kFixedSizeWidth);
    std::fseek(impl_->file, impl_->segment_size_at, SEEK_SET);
    impl_->write(size_bytes);

    Bytes table;
    for (std::size_t i = 0; i < impl_->seeks.size(); ++i) {
        Bytes entry_body;
        put_id(entry_body, kSeekId);
        Bytes target;
        put_id(target, kAttachedFile);
        put_size(entry_body, target.size());
        entry_body.insert(entry_body.end(), target.begin(), target.end());
        put_uint(entry_body, kSeekPosition,
                 static_cast<unsigned long long>(impl_->seeks[i]), 8);
        put_id(table, kSeek);
        put_size(table, entry_body.size());
        table.insert(table.end(), entry_body.begin(), entry_body.end());
    }
    const Bytes tail = void_bytes(impl_->seek_head_bytes - table.size());
    table.insert(table.end(), tail.begin(), tail.end());
    std::fseek(impl_->file, impl_->seek_head_data_at, SEEK_SET);
    impl_->write(table);

    std::fclose(impl_->file);
    impl_->file = NULL;
    impl_->closed = true;
}

// --------------------------------------------------------------------------
// the reader
// --------------------------------------------------------------------------


struct PtoReader::Impl {
    //! An element header: its ID, its payload's offset and its payload's size.
    struct Header {
        unsigned int id;
        long body;
        unsigned long long size;
    };

    //! What one buffered read pulls in. An `AttachedFile`'s framing -- its
    //! own header, then uid, name, kind, encoding and the `FileData` header
    //! -- is a couple of hundred octets, so a window this size turns the ten
    //! reads that walking one object would take into one. It is small
    //! deliberately: the payloads it must not drag in are up to a hundred
    //! thousand times larger than the framing between them.
    static const std::size_t kBlock = 512;

    std::FILE* file;
    std::string path;
    long end;
    //! The cached window: `block` holds `block_size` octets from `block_at`.
    std::vector<unsigned char> block;
    long block_at;
    std::size_t block_size;

    Impl() : file(NULL), end(0), block(kBlock), block_at(-1), block_size(0) {}
    ~Impl() { if (file != NULL) std::fclose(file); }

    //! Decode the header at `at`. Reads at most 12 octets -- an ID is up to
    //! four and a Data Size up to eight -- so walking a document of a hundred
    //! objects costs a kilobyte of reads, not the document.
    bool read_header(long at, long limit, Header* out) {
        if (at + 2 > limit) return false;
        unsigned char buf[12];
        std::size_t n = static_cast<std::size_t>(limit - at);
        if (n > sizeof(buf)) n = sizeof(buf);
        read_at(at, buf, n);
        const std::size_t id_len = read_id(buf, n, &out->id);
        if (id_len == 0) return false;
        const std::size_t size_len =
                read_size(buf + id_len, n - id_len, &out->size);
        if (size_len == 0) return false;
        out->body = at + static_cast<long>(id_len + size_len);
        return true;
    }

    //! Read `n` bytes at `at`; the one place this reader touches the disk.
    /*! Small reads come out of the cached window and refill it when they
        miss; a read larger than the window goes straight to the file, so a
        payload is never copied twice. Walking a container of a few hundred
        objects is thousands of these, and the cost of crossing the framing
        was the cost of its syscalls, not of its bytes. */
    void read_at(long at, void* into, std::size_t n) {
        if (n == 0) return;
        if (n <= kBlock) {
            const bool hit = block_at >= 0 && at >= block_at &&
                    at + static_cast<long>(n) <=
                            block_at + static_cast<long>(block_size);
            if (!hit) {
                std::size_t want = kBlock;
                if (end > 0 && at + static_cast<long>(want) > end) {
                    want = static_cast<std::size_t>(end - at);
                }
                if (want < n) want = n;
                if (std::fseek(file, at, SEEK_SET) != 0) {
                    IMP_THROW("PTO: seeking " << path << " to " << at
                              << " failed", IOException);
                }
                block_size = std::fread(&block[0], 1, want, file);
                if (block_size < n) {
                    block_at = -1;
                    IMP_THROW("PTO: reading " << path << " at " << at
                              << " failed", IOException);
                }
                block_at = at;
            }
            std::memcpy(into, &block[at - block_at], n);
            return;
        }
        if (std::fseek(file, at, SEEK_SET) != 0 ||
            std::fread(into, 1, n, file) != n) {
            IMP_THROW("PTO: reading " << path << " at " << at << " failed",
                      IOException);
        }
    }
};

bool PtoReader::looks_like_pto(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (f == NULL) return false;
    unsigned char magic[4] = {0, 0, 0, 0};
    const std::size_t got = std::fread(magic, 1, 4, f);
    std::fclose(f);
    return got == 4 && magic[0] == 0x1A && magic[1] == 0x45 &&
           magic[2] == 0xDF && magic[3] == 0xA3;
}

PtoReader::PtoReader(const std::string& path)
        : impl_(new Impl), doctype_version_(0) {
    impl_->path = path;
    impl_->file = std::fopen(path.c_str(), "rb");
    if (impl_->file == NULL) {
        delete impl_;
        impl_ = NULL;
        IMP_THROW("PTO: cannot open " << path, IOException);
    }
    std::fseek(impl_->file, 0, SEEK_END);
    const long end = std::ftell(impl_->file);
    impl_->end = end;

    unsigned char magic[4] = {0, 0, 0, 0};
    if (end < 4) IMP_THROW("PTO: " << path << " is too short", IOException);
    impl_->read_at(0, magic, 4);
    if (magic[0] != 0x1A || magic[1] != 0x45 || magic[2] != 0xDF ||
        magic[3] != 0xA3) {
        IMP_THROW("PTO: " << path << " does not begin with an EBML header",
                  IOException);
    }

    // Walk the framing. Master elements whose children matter are descended
    // into; everything else is stepped over by its own size, which is what
    // makes an element this reader has never heard of harmless.
    std::string doctype;
    long pos = 0;
    while (pos < end) {
        Impl::Header h;
        if (!impl_->read_header(pos, end, &h)) break;
        if (h.body + static_cast<long>(h.size) > end) {
            IMP_THROW("PTO: " << path << " is truncated", IOException);
        }

        if (h.id == kSegment || h.id == kAttachments) {
            pos = h.body;                                     // descend
            continue;
        }
        if (h.id == kEbmlHeader || h.id == kAttachedFile) {
            const long stop = h.body + static_cast<long>(h.size);
            PtoObject object;
            long at = h.body;
            while (at < stop) {
                Impl::Header c;
                if (!impl_->read_header(at, stop, &c)) break;
                const std::size_t n = static_cast<std::size_t>(c.size);
                if (h.id == kEbmlHeader) {
                    if (c.id == kDocType || c.id == kDocTypeVersion) {
                        std::vector<unsigned char> v(n);
                        impl_->read_at(c.body, n ? &v[0] : NULL, n);
                        if (c.id == kDocType) doctype = as_string(n ? &v[0] : NULL, n);
                        else doctype_version_ = static_cast<int>(as_uint(&v[0], n));
                    }
                } else if (c.id == kFileData) {
                    object.offset = static_cast<std::size_t>(c.body);
                    object.size = n;
                } else if (c.id == kFileUid || c.id == kFileName ||
                           c.id == kPtoKind || c.id == kPtoEncoding) {
                    std::vector<unsigned char> v(n);
                    impl_->read_at(c.body, n ? &v[0] : NULL, n);
                    const unsigned char* p = n ? &v[0] : NULL;
                    if (c.id == kFileUid) object.uid = as_uint(p, n);
                    else if (c.id == kFileName) object.name = as_string(p, n);
                    else if (c.id == kPtoKind) object.kind = as_string(p, n);
                    else object.encoding = as_string(p, n);
                }
                at = c.body + static_cast<long>(c.size);
            }
            if (h.id == kAttachedFile) objects_.push_back(object);
        }
        pos = h.body + static_cast<long>(h.size);
    }

    if (doctype != "pto") {
        IMP_THROW("PTO: " << path << " has DocType '" << doctype
                          << "', not 'pto'", IOException);
    }
}

PtoReader::~PtoReader() {
    delete impl_;
    impl_ = NULL;
}

int PtoReader::find(const std::string& name) const {
    for (std::size_t i = 0; i < objects_.size(); ++i) {
        if (objects_[i].name == name) return static_cast<int>(i);
    }
    return -1;
}

std::vector<unsigned char> PtoReader::data(const PtoObject& object) const {
    std::vector<unsigned char> out(object.size);
    if (object.size != 0) {
        impl_->read_at(static_cast<long>(object.offset), &out[0], object.size);
    }
    return out;
}

IMPBFF_END_NAMESPACE

// -------- from PtoProfile.cpp --------
/**
 * (formerly PtoProfile.cpp, now a section of this file)
 * \brief PTO.MFDB — the vocabulary layer over the PTO container.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/bff_config.h>
#include <IMP/bff/internal/Sha256.h>
#include <IMP/bff/internal/json.h>


#include <algorithm>
#include <map>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

const char* const MFDB_PROFILE = "PTO.MFDB";
const char* const MFDB_PROFILE_VERSION = "1.1";
const char* const MFDB_PROFILE_READ_VERSION = "1";
const char* const MFDB_DICTIONARY_VERSION = "1.8";

namespace {

std::vector<std::string> pp_list(const char* const* v, std::size_t n) {
    return std::vector<std::string>(v, v + n);
}

}  // namespace

// ---------------------------------------------------------------------------
// The controlled vocabularies
// ---------------------------------------------------------------------------
//
// Transcribed from `../mmfdb/src/mmfdb/data/mmfdb_flr_ext.dic`, version 1.8.
// Only the values this package can legitimately write are listed: the
// dictionary's `artifact_kind` carries seventy terms for photon streams,
// imaging and project archives that nothing here produces, and listing them
// would invite a caller to pick one.

const std::vector<std::string>& mfdb_artifact_kinds() {
    static const char* const v[] = {
        "readme", "sample_metadata", "parameter_table", "analysis_result",
        "processed_data", "derived_product", "row_mapping",
        "external_reference", "json_summary", "table_export"};
    static const std::vector<std::string> l = pp_list(v, 10);
    return l;
}

const std::vector<std::string>& mfdb_data_formats() {
    static const char* const v[] = {"json", "csv",  "tsv",  "cif",
                                    "text", "bin",  "npy",  "dstore",
                                    "pto",  "unknown"};
    static const std::vector<std::string> l = pp_list(v, 10);
    return l;
}

const std::vector<std::string>& mfdb_row_grains() {
    static const char* const v[] = {
        "photon", "burst", "dwell",  "segment", "pixel", "voxel",
        "frame",  "line",  "molecule", "track", "spot",  "region",
        "curve_point", "state", "species", "channel", "pair", "file",
        "histogram_bin", "spectrum",
        // Added in dictionary 1.8 for this package: a position on a polymer
        // that a label could be attached to, the thing
        // `_flr_poly_probe_position` addresses. Nothing existing fitted --
        // `species` and `state` count what a molecule is, not a place on one.
        "label_site"};
    static const std::vector<std::string> l = pp_list(v, 21);
    return l;
}

const std::vector<std::string>& mfdb_relationship_types() {
    static const char* const v[] = {
        "included_in", "contains", "derived_from", "supersedes",
        "uses_external_reference", "parameter_depends_on", "parameter_of",
        "linked_to", "project_contains", "grouped_in", "measured_sample",
        "calibrated_by", "maps_rows_of"};
    static const std::vector<std::string> l = pp_list(v, 13);
    return l;
}

const std::vector<std::string>& mfdb_units() {
    static const char* const v[] = {
        "hours", "minutes", "seconds", "milliseconds", "microseconds",
        "nanoseconds", "picoseconds", "femtoseconds", "hertz", "kilohertz",
        "megahertz", "counts", "counts_per_second", "photons", "nanometres",
        "micrometres", "angstroms", "pixels", "degrees", "radians", "celsius",
        "kelvins", "molar", "millimolar", "micromolar", "nanomolar",
        "picomolar", "dimensionless"};
    static const std::vector<std::string> l = pp_list(v, 28);
    return l;
}

const std::vector<std::string>& mfdb_operation_types() {
    // A subset of the dictionary's enumeration: what a writer in *this*
    // package produces. A reader must expect the rest -- most of that
    // enumeration is photon-level and belongs to tttrlib and chisurf.
    //
    // `clustering` is not `ndxplorer_clustering`, which is one tool's
    // burst-selection step. This one covers grouping conformations: a rotamer
    // library built from a trajectory, MSM microstate construction,
    // `cluster_frames_leader`. Added for the .drot libraries, which have to be
    // able to say *which* protocol built them -- dihedral k-means and leader
    // RMSD do not produce comparable libraries, and recording `analysis` for
    // both would be a lie that nothing ever checks.
    static const char* const v[] = {"analysis", "analysis_run", "validation",
                                    "model_fitting", "external_tool",
                                    "import", "calibration", "clustering"};
    static const std::vector<std::string> l = pp_list(v, 8);
    return l;
}

const std::vector<std::string>& mfdb_optical_property_names() {
    static const char* const v[] = {
        "quantum_yield", "extinction_coefficient", "excitation_maximum",
        "emission_maximum", "fluorescence_lifetime", "anisotropy_fundamental",
        "steric_radius", "hydrodynamic_radius"};
    static const std::vector<std::string> l = pp_list(v, 8);
    return l;
}

const std::vector<std::string>& mfdb_optical_property_units() {
    static const char* const v[] = {"dimensionless", "nanometres", "angstroms",
                                    "nanoseconds", "per_molar_per_centimetre"};
    static const std::vector<std::string> l = pp_list(v, 5);
    return l;
}

const std::vector<std::string>& mfdb_spectrum_types() {
    static const char* const v[] = {"absorption", "excitation", "emission",
                                    "transmission", "reflectance",
                                    "quantum_efficiency", "responsivity"};
    static const std::vector<std::string> l = pp_list(v, 7);
    return l;
}

const std::vector<std::string>& mfdb_spectrum_wavelength_units() {
    // `nm`, deliberately: every stored spectrum in this stack is labelled that
    // way and the dictionary follows the data rather than renaming 2896 rows.
    static const char* const v[] = {"nm", "um"};
    static const std::vector<std::string> l = pp_list(v, 2);
    return l;
}

const std::vector<std::string>& mfdb_spectrum_intensity_units() {
    static const char* const v[] = {"normalized", "counts",
                                    "per_molar_per_centimetre",
                                    "dimensionless"};
    static const std::vector<std::string> l = pp_list(v, 4);
    return l;
}

const std::vector<std::string>& mfdb_label_score_types() {
    static const char* const v[] = {
        "conservation", "solvent_exposure", "secondary_structure",
        "charge_environment", "tryptophan_proximity", "cysteine_resemblance",
        "methionine_exclusion", "fret_sensitivity", "measurement", "combined"};
    static const std::vector<std::string> l = pp_list(v, 10);
    return l;
}

const std::vector<std::string>& mfdb_label_score_statuses() {
    static const char* const v[] = {"scored", "excluded", "unresolved",
                                    "unavailable"};
    static const std::vector<std::string> l = pp_list(v, 4);
    return l;
}

const std::vector<std::string>& mfdb_probe_types() {
    static const char* const v[] = {"dye", "fluorescent_protein", "spin_label",
                                    "unspecified"};
    static const std::vector<std::string> l = pp_list(v, 4);
    return l;
}

const std::vector<std::string>& mfdb_label_score_definitions() {
    static const char* const v[] = {"labelizer", "consurf", "dssp", "msms"};
    static const std::vector<std::string> l = pp_list(v, 4);
    return l;
}

namespace {

const std::vector<std::string>* pp_enumeration(const std::string& item) {
    if (item == "_mmfdb_artifact.artifact_kind") return &mfdb_artifact_kinds();
    if (item == "_mmfdb_artifact.data_format") return &mfdb_data_formats();
    if (item == "_mmfdb_artifact.row_grain") return &mfdb_row_grains();
    if (item == "_mmfdb_edge.relationship_type") return &mfdb_relationship_types();
    if (item == "_mmfdb_column.units") return &mfdb_units();
    if (item == "_mmfdb_operation.operation_type") return &mfdb_operation_types();
    if (item == "_mmfdb_label_score.score_type") return &mfdb_label_score_types();
    if (item == "_mmfdb_label_score.status") return &mfdb_label_score_statuses();
    if (item == "_mmfdb_label_score.definition")
        return &mfdb_label_score_definitions();
    if (item == "_mmfdb_optical_property.property_name")
        return &mfdb_optical_property_names();
    if (item == "_mmfdb_optical_property.unit")
        return &mfdb_optical_property_units();
    if (item == "_mmfdb_probe.probe_type") return &mfdb_probe_types();
    if (item == "_mmfdb_spectrum.spectrum_type") return &mfdb_spectrum_types();
    if (item == "_mmfdb_spectrum.wavelength_unit")
        return &mfdb_spectrum_wavelength_units();
    if (item == "_mmfdb_spectrum.intensity_unit")
        return &mfdb_spectrum_intensity_units();
    return 0;
}

}  // namespace

bool mfdb_is_term(const std::string& item, const std::string& value) {
    const std::vector<std::string>* e = pp_enumeration(item);
    if (!e) return false;
    return std::find(e->begin(), e->end(), value) != e->end();
}

void mfdb_check_term(const std::string& item, const std::string& value) {
    const std::vector<std::string>* e = pp_enumeration(item);
    if (!e) {
        IMP_THROW("mfdb_check_term: no enumeration is declared here for "
                  << item << "; add it from mmfdb_flr_ext.dic rather than "
                  << "writing an unchecked value", ValueException);
    }
    if (std::find(e->begin(), e->end(), value) != e->end()) return;
    std::ostringstream allowed;
    for (std::size_t i = 0; i < e->size(); ++i) {
        allowed << (i ? ", " : "") << (*e)[i];
    }
    IMP_THROW("mfdb_check_term: '" << value << "' is not a value of " << item
              << ". Allowed: " << allowed.str()
              << ". A term that is genuinely missing belongs in "
              << "mmfdb_flr_ext.dic first, not here.", ValueException);
}

// ---------------------------------------------------------------------------
// Columns
// ---------------------------------------------------------------------------

std::string mfdb_columns_json(const std::vector<MfdbColumn>& columns) {
    nlohmann::json j = nlohmann::json::array();
    for (std::size_t i = 0; i < columns.size(); ++i) {
        const MfdbColumn& c = columns[i];
        nlohmann::json e;
        e["name"] = c.name;
        // No unit means the unit is unknown, which is a different claim from
        // `dimensionless`; an empty string is therefore omitted, not written.
        if (!c.units.empty()) {
            mfdb_check_term("_mmfdb_column.units", c.units);
            e["units"] = c.units;
        }
        if (!c.item.empty()) e["item"] = c.item;
        if (!c.description.empty()) e["description"] = c.description;
        j.push_back(e);
    }
    return j.dump();
}

// ---------------------------------------------------------------------------
// The tag sets
// ---------------------------------------------------------------------------

std::vector<MfdbTag> mfdb_container_tags() {
    std::vector<MfdbTag> t;
    t.push_back(MfdbTag("_mmfdb_container.profile", MFDB_PROFILE));
    t.push_back(MfdbTag("_mmfdb_container.profile_version", MFDB_PROFILE_VERSION));
    t.push_back(MfdbTag("_mmfdb_container.profile_read_version",
                        MFDB_PROFILE_READ_VERSION));
    t.push_back(MfdbTag("_mmfdb_container.format", "pto"));
    t.push_back(MfdbTag("_mmfdb_container.dictionary_version",
                        MFDB_DICTIONARY_VERSION));
    return t;
}

std::vector<MfdbTag> mfdb_artifact_tags(const std::string& artifact_id,
                                        const std::string& artifact_kind,
                                        const std::string& data_format,
                                        const std::string& row_grain,
                                        long row_count,
                                        const std::string& checksum,
                                        const std::vector<MfdbColumn>& columns) {
    mfdb_check_term("_mmfdb_artifact.artifact_kind", artifact_kind);
    mfdb_check_term("_mmfdb_artifact.data_format", data_format);
    if (!row_grain.empty()) {
        mfdb_check_term("_mmfdb_artifact.row_grain", row_grain);
    }
    std::vector<MfdbTag> t;
    t.push_back(MfdbTag("_mmfdb_artifact.artifact_id", artifact_id));
    t.push_back(MfdbTag("_mmfdb_artifact.artifact_kind", artifact_kind));
    t.push_back(MfdbTag("_mmfdb_artifact.data_format", data_format));
    if (!row_grain.empty()) {
        t.push_back(MfdbTag("_mmfdb_artifact.row_grain", row_grain));
    }
    if (row_count >= 0) {
        std::ostringstream os;
        os << row_count;
        t.push_back(MfdbTag("_mmfdb_artifact.row_count", os.str()));
    }
    if (!checksum.empty()) {
        t.push_back(MfdbTag("_mmfdb_artifact.checksum", checksum));
        t.push_back(MfdbTag("_mmfdb_artifact.checksum_algorithm", "sha256"));
    }
    if (!columns.empty()) {
        t.push_back(MfdbTag("_mmfdb_column", mfdb_columns_json(columns)));
    }
    return t;
}

std::vector<MfdbTag> mfdb_operation_tags(const std::string& operation_type,
                                         const std::string& algorithm,
                                         const std::string& settings_json,
                                         const std::string& software_package,
                                         const std::string& software_version) {
    mfdb_check_term("_mmfdb_operation.operation_type", operation_type);
    std::vector<MfdbTag> t;
    t.push_back(MfdbTag("_mmfdb_operation.operation_type", operation_type));
    if (!algorithm.empty()) {
        t.push_back(MfdbTag("_mmfdb_operation.algorithm", algorithm));
    }
    t.push_back(MfdbTag("_mmfdb_operation.settings_json", settings_json));
    // The identity of a run: same settings, same artifact.
    t.push_back(MfdbTag("_mmfdb_operation.settings_hash",
                        mfdb_checksum(settings_json)));
    t.push_back(MfdbTag("_mmfdb_operation.software_package", software_package));
    t.push_back(MfdbTag("_mmfdb_operation.software_version", software_version));
    t.push_back(MfdbTag("_mmfdb_operation.dictionary_version",
                        MFDB_DICTIONARY_VERSION));
    return t;
}

std::vector<MfdbTag> mfdb_edge_tags(const std::string& source_node_id,
                                    const std::string& target_node_id,
                                    const std::string& relationship_type,
                                    const std::string& source_row_column,
                                    const std::string& target_row_column) {
    mfdb_check_term("_mmfdb_edge.relationship_type", relationship_type);
    std::vector<MfdbTag> t;
    t.push_back(MfdbTag("_mmfdb_edge.source_node_id", source_node_id));
    t.push_back(MfdbTag("_mmfdb_edge.target_node_id", target_node_id));
    t.push_back(MfdbTag("_mmfdb_edge.relationship_type", relationship_type));
    if (!source_row_column.empty()) {
        t.push_back(MfdbTag("_mmfdb_edge.source_row_column", source_row_column));
    }
    if (!target_row_column.empty()) {
        t.push_back(MfdbTag("_mmfdb_edge.target_row_column", target_row_column));
    }
    return t;
}

std::vector<MfdbTag> mfdb_attribution_tags(const MfdbAttribution& a) {
    if (a.license.empty() && a.terms.empty()) {
        IMP_THROW("mfdb_attribution_tags: an artifact needs either a license "
                  "(an SPDX identifier) or terms (the conditions in words). "
                  "Data whose terms are unknown is not the same as data that "
                  "is unrestricted, and the difference matters to whoever "
                  "redistributes it.", ValueException);
    }
    std::vector<MfdbTag> t;
    if (!a.author.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.author", a.author));
    if (!a.citation.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.citation", a.citation));
    if (!a.license.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.license", a.license));
    if (!a.terms.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.terms", a.terms));
    if (!a.terms_url.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.terms_url", a.terms_url));
    if (!a.source.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.source", a.source));
    if (!a.redistributed_via.empty())
        t.push_back(MfdbTag("_mmfdb_artifact.redistributed_via",
                            a.redistributed_via));
    return t;
}

std::string mfdb_checksum(const std::string& bytes) {
    return internal::sha256_hex(bytes);
}

IMPBFF_END_NAMESPACE
