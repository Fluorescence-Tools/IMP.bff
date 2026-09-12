/** \file ProteinSidechainRotamerLibrary.cpp
 * \brief Dunbrack side-chain tables in PTO containers.
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/ProteinSidechainRotamerLibrary.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/ptolib.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

// ===========================================================================
// The backbone-dependent side-chain library: FASPR's binary, in the container
// ===========================================================================


namespace {

//! One record: float32 probability, then chi1..4 and sigma1..4 as int16.
const std::size_t kRecordBytes = 20;
//! 36 phi bins x 36 psi bins per residue.
const std::size_t kBinsPerResidue = 36 * 36;

//! The residues that have rotamers, in the order the binary stores them.
/*! `offset` counts (phi, psi) bin blocks from the start of the file, and
    `n_rotamers` is how many slots one bin holds -- both are FASPR's tables,
    and together they say where a residue's records begin and end. */
struct Residue {
    char code;
    int n_chi;
    int n_rotamers;
    int offset;
};

const Residue kResidues[] = {
    {'R', 4, 75,   0}, {'N', 2, 36,  75}, {'D', 2, 18, 111},
    {'C', 1,  3, 129}, {'Q', 3, 108, 132}, {'E', 3, 54, 240},
    {'H', 2, 36, 294}, {'I', 2,  9, 330}, {'L', 2,  9, 339},
    {'K', 4, 73, 348}, {'M', 3, 27, 421}, {'F', 2, 18, 448},
    {'P', 2,  2, 466}, {'S', 1,  3, 468}, {'T', 1,  3, 471},
    {'W', 2, 36, 474}, {'Y', 2, 18, 510}, {'V', 1,  3, 528},
};
const std::size_t kResidueCount = sizeof(kResidues) / sizeof(kResidues[0]);
//! What the whole library weighs: the last residue's end, in records.
const std::size_t kTotalBins = 531;   // 528 + 3, V's offset plus its rotamers

const Residue* residue(char code) {
    for (std::size_t i = 0; i < kResidueCount; ++i) {
        if (kResidues[i].code == code) return &kResidues[i];
    }
    return NULL;
}

std::size_t residue_bytes(const Residue& r) {
    return kBinsPerResidue * static_cast<std::size_t>(r.n_rotamers) *
           kRecordBytes;
}

std::size_t residue_start(const Residue& r) {
    return static_cast<std::size_t>(r.offset) * kBinsPerResidue * kRecordBytes;
}

//! FASPR's bin: 10 degrees wide, centred, wrapping at 36.
int bin_of(double angle_deg) {
    const double shifted = angle_deg + 180.0;
    int bin = static_cast<int>(std::floor((shifted + 5.0) / 10.0));
    bin %= 36;
    if (bin < 0) bin += 36;
    return bin;
}

std::vector<unsigned char> brotli_pack(const std::vector<unsigned char>& in) {
    std::vector<unsigned char> out;
    if (!pto::compress_bytes("brotli", in.empty() ? NULL : &in[0], in.size(),
                             11, out)) {
        IMP_THROW("dunbrack: brotli compression failed", IOException);
    }
    return out;
}

std::vector<unsigned char> brotli_unpack(const std::vector<unsigned char>& in,
                                         std::size_t expected) {
    std::vector<unsigned char> out;
    if (!pto::decompress_bytes("brotli", in.empty() ? NULL : &in[0], in.size(),
                               expected, out)) {
        IMP_THROW("dunbrack: the payload did not decompress to its stated "
                  << expected << " bytes", IOException);
    }
    return out;
}

std::string object_name(char code) {
    return std::string(1, code) + "/bbdep.records";
}

}  // namespace

void ProteinSidechainDunbrackRotamers::get_chi(double** out_view, int* n_out_view) const {
    internal::copy_to_view(chi, out_view, n_out_view);
}
void ProteinSidechainDunbrackRotamers::get_sigma(double** out_view, int* n_out_view) const {
    internal::copy_to_view(sigma, out_view, n_out_view);
}
void ProteinSidechainDunbrackRotamers::get_probability(double** out_view,
                                       int* n_out_view) const {
    internal::copy_to_view(probability, out_view, n_out_view);
}

void write_protein_sidechain_dunbrack_library(const std::string& bin_path,
                            const std::string& path) {
    std::FILE* f = std::fopen(bin_path.c_str(), "rb");
    if (f == NULL) {
        IMP_THROW("write_dunbrack_library: cannot open " << bin_path,
                  IOException);
    }
    std::fseek(f, 0, SEEK_END);
    const long size = std::ftell(f);
    const std::size_t expected = kTotalBins * kBinsPerResidue * kRecordBytes;
    if (static_cast<std::size_t>(size) != expected) {
        std::fclose(f);
        IMP_THROW("write_dunbrack_library: " << bin_path << " is " << size
                  << " bytes, not the " << expected
                  << " a dun2010bbdep.bin has -- is it the FASPR library?",
                  IOException);
    }

    pto::File pto;
    if (!pto.create(path, "",
            std::string(pto::kDefaultBanner) + "\nThis container was written by IMP.bff.\nhttps://github.com/tpeulen/IMP.bff\n")) {
        IMP_THROW("PTO: cannot open " << path << " for writing: "
                  << pto.error(), IOException);
    }
    pto.set_writing_app("IMP.bff");

    std::ostringstream hj;
    hj << "{\"format\":\"rot.bbdep\",\"version\":1,\"container\":\"pto\""
       << ",\"producer\":\"IMP.bff write_dunbrack_library\""
       << ",\"source\":\"dun2010bbdep.bin\""
       << ",\"library\":\"Dunbrack 2010 backbone-dependent rotamer library\""
       << ",\"citation\":\"Shapovalov MV, Dunbrack RL. Structure "
          "2011;19:844-858\""
       << ",\"redistributed_via\":\"FASPR (Huang X, Pearce R, Zhang Y. "
          "Bioinformatics 2020;36:3758-3765), MIT\""
       << ",\"terms\":\"The rotamer library is free for academic use; this "
          "container carries it unchanged and this notice with it.\""
       << ",\"record\":{\"bytes\":20,\"fields\":[\"f32 probability\","
          "\"i16 chi1..chi4\",\"i16 sigma1..sigma4\"],"
          "\"angle_unit\":\"tenths of a degree\"}"
       << ",\"bins\":{\"phi\":36,\"psi\":36,\"width_deg\":10,"
          "\"index\":\"floor((angle + 185) / 10) mod 36\"}"
       << ",\"residues\":{";
    for (std::size_t i = 0; i < kResidueCount; ++i) {
        if (i) hj << ',';
        hj << '"' << kResidues[i].code << "\":{\"n_chi\":"
           << kResidues[i].n_chi << ",\"n_rotamers\":"
           << kResidues[i].n_rotamers << ",\"bin_offset\":"
           << kResidues[i].offset << '}';
    }
    hj << "}}";
    const std::string header = hj.str();
    std::vector<unsigned char> hb(header.begin(), header.end());
    const std::vector<unsigned char> hp = brotli_pack(hb);
    if (pto.add("rot.bbdep.header", "json+brotli", "bbdep.json", hp.empty() ? NULL : &hp[0], hp.size()) == 0) {
        IMP_THROW("PTO: writing " << pto.filename() << " failed: "
                  << pto.error(), IOException);
    }

    std::ostringstream cj;
    cj << "{\"format\":\"drot.bundle\",\"version\":10"
       << ",\"producer\":\"IMP.bff write_dunbrack_library\""
       << ",\"n_libraries\":" << kResidueCount << ",\"libraries\":[";
    for (std::size_t i = 0; i < kResidueCount; ++i) {
        if (i) cj << ',';
        cj << '"' << kResidues[i].code << '"';
    }
    cj << "]}";
    const std::string catalog = cj.str();
    std::vector<unsigned char> cb(catalog.begin(), catalog.end());
    const std::vector<unsigned char> cp = brotli_pack(cb);
    if (pto.add("drot.catalog", "json+brotli", "drot.catalog", cp.empty() ? NULL : &cp[0], cp.size()) == 0) {
        IMP_THROW("PTO: writing " << pto.filename() << " failed: "
                  << pto.error(), IOException);
    }

    for (std::size_t i = 0; i < kResidueCount; ++i) {
        const Residue& r = kResidues[i];
        const std::size_t n = residue_bytes(r);
        std::vector<unsigned char> records(n);
        if (std::fseek(f, static_cast<long>(residue_start(r)), SEEK_SET) != 0 ||
            std::fread(&records[0], 1, n, f) != n) {
            std::fclose(f);
            IMP_THROW("write_dunbrack_library: reading " << bin_path
                      << " for residue " << r.code << " failed", IOException);
        }
        const std::vector<unsigned char> packed = brotli_pack(records);
        if (pto.add("rot.bbdep.records", "faspr20+brotli", object_name(r.code), packed.empty() ? NULL : &packed[0], packed.size()) == 0) {
            IMP_THROW("PTO: writing " << pto.filename() << " failed: "
                      << pto.error(), IOException);
        }
    }
    std::fclose(f);
    const bool committed = pto.commit();
    const std::string commit_error = pto.error();
    pto.close();
    if (!committed) IMP_THROW("PTO: writing " << path << " failed: "
                             << commit_error, IOException);
}

void write_protein_sidechain_dunbrack_bin(const std::string& path, const std::string& bin_path) {
    pto::File pto;
    if (!pto.open(path, false)) {
        IMP_THROW("PTO: cannot open container " << path << ": "
                  << (pto.error().empty() ? "not an EBML document with DocType pto"
                                            : pto.error()), IOException);
    }
    std::FILE* out = std::fopen(bin_path.c_str(), "wb");
    if (out == NULL) {
        IMP_THROW("write_dunbrack_bin: cannot open " << bin_path
                  << " for writing", IOException);
    }
    for (std::size_t i = 0; i < kResidueCount; ++i) {
        const Residue& r = kResidues[i];
        const auto matches = pto.find_all(object_name(r.code));
        if (matches.empty()) {
            std::fclose(out);
            IMP_THROW("write_dunbrack_bin: " << path << " has no records for "
                      << r.code, IOException);
        }
        // ptolib's reader decodes by the object's encoding.
        const std::vector<unsigned char> records = pto.read(matches.front());
        if (std::fwrite(&records[0], 1, records.size(), out) !=
            records.size()) {
            std::fclose(out);
            IMP_THROW("write_dunbrack_bin: short write to " << bin_path,
                      IOException);
        }
    }
    std::fclose(out);
}

ProteinSidechainDunbrackRotamers read_protein_sidechain_dunbrack_rotamers(const std::string& path, char residue_,
                                        double phi, double psi,
                                        double probability_min,
                                        double probability_accumulated) {
    const Residue* r = residue(residue_);
    if (r == NULL) {
        IMP_THROW("read_dunbrack_rotamers: '" << residue_
                  << "' has no backbone-dependent rotamers (alanine and "
                     "glycine have no chi; is it a one-letter code?)",
                  ValueException);
    }
    pto::File pto;
    if (!pto.open(path, false)) {
        IMP_THROW("PTO: cannot open container " << path << ": "
                  << (pto.error().empty() ? "not an EBML document with DocType pto"
                                            : pto.error()), IOException);
    }
    const auto matches = pto.find_all(object_name(r->code));
    if (matches.empty()) {
        IMP_THROW("read_dunbrack_rotamers: " << path << " has no records for "
                  << r->code, IOException);
    }
    // ptolib's reader decodes by the object's encoding; brotli_unpack's
    // second pass here decompressed plain records and failed.
    const std::vector<unsigned char> records = pto.read(matches.front());
    if (records.size() != static_cast<std::size_t>(r->n_rotamers) *
                                   residue_bytes(*r)) {
        IMP_THROW("read_dunbrack_rotamers: " << path << " holds "
                  << records.size() << " bytes for " << r->code
                  << ", expected " << r->n_rotamers * residue_bytes(*r),
                  IOException);
    }

    // Where this (phi, psi) bin's slots begin, in records.
    const std::size_t bin = static_cast<std::size_t>(36 * bin_of(phi) +
                                                     bin_of(psi));
    const std::size_t first = bin * static_cast<std::size_t>(r->n_rotamers);

    ProteinSidechainDunbrackRotamers out;
    out.n_chi = r->n_chi;
    double accumulated = 0.0;
    for (int j = 0; j < r->n_rotamers; ++j) {
        const unsigned char* rec = &records[(first + j) * kRecordBytes];
        float p;
        std::memcpy(&p, rec, sizeof(p));
        if (p < probability_min) break;
        out.probability.push_back(p);
        for (int k = 0; k < r->n_chi; ++k) {
            short chi_i, sigma_i;
            std::memcpy(&chi_i, rec + 4 + 2 * k, sizeof(chi_i));
            std::memcpy(&sigma_i, rec + 4 + 2 * (4 + k), sizeof(sigma_i));
            out.chi.push_back(chi_i / 10.0);
            out.sigma.push_back(sigma_i / 10.0);
        }
        ++out.n_rotamers;
        accumulated += p;
        if (accumulated > probability_accumulated) break;
    }
    return out;
}

IMPBFF_END_NAMESPACE
