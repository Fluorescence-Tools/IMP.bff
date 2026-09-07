/**
 * \file HierarchyFrame.cpp
 * \brief Reading one frame out of an IMP hierarchy, in one call.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/HierarchyFrame.h>
#include <IMP/bff/internal/PdbFrames.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/Text.h>

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>



IMPBFF_BEGIN_NAMESPACE


namespace {
//! The XYZ leaves, in hierarchy order.
}  // namespace


// --------------------------------------------------------------------------
// ProteinFrame
// --------------------------------------------------------------------------

void ProteinFrame::get_coords(double** out_view, int* n_out_view) const {
    internal::copy_to_view(coords, out_view, n_out_view);
}


// -------- the PDB reader: IMP's rules, no IMP (internal/PdbFrames.h) --------
namespace {

bool starts_with(const std::string& line, const char* rec) {
    return line.compare(0, std::strlen(rec), rec) == 0;
}

std::string field(const std::string& line, std::size_t pos, std::size_t n) {
    return pos < line.size() ? line.substr(pos, n) : std::string();
}

// The name the hierarchy path reads back: the second whitespace token of
// "Atom <type> of residue <n>". Kept as that string operation rather than
// simplified, so a type with a space in it ("HET: S3 ") answers as it did.
std::string name_from_particle_name(const std::string& full) {
    const std::size_t first = full.find_first_not_of(" \t");
    const std::size_t gap = first == std::string::npos ? std::string::npos
                                                       : full.find_first_of(" \t", first);
    if (gap == std::string::npos) {
        return first == std::string::npos ? full : full.substr(first);
    }
    const std::size_t start = full.find_first_not_of(" \t", gap);
    const std::size_t end = start == std::string::npos ? std::string::npos
                                                       : full.find_first_of(" \t", start);
    return start == std::string::npos
                   ? full.substr(first, gap - first)
                   : full.substr(start, end == std::string::npos ? std::string::npos
                                                                 : end - start);
}

}  // namespace

std::vector<ProteinFrame> internal::read_pdb_frames(
        const std::string& path, internal::PdbSelect select,
        bool first_model_only, int max_frames) {
    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("Cannot open " << path, IMP::IOException);
    std::vector<ProteinFrame> frames;
    bool open_frame = false;       // a frame exists for the atoms that follow
    bool first_model_read = false;
    std::string line;
    while (std::getline(in, line)) {
        if (starts_with(line, "MODEL")) {
            if (first_model_read && first_model_only) break;
            open_frame = false;    // the next accepted atom starts a new frame
            first_model_read = true;
            continue;
        }
        const bool is_atom = starts_with(line, "ATOM");
        const bool is_hetatm = starts_with(line, "HETATM");
        if (!is_atom && !is_hetatm) continue;
        if (select == internal::PDB_NON_WATER) {
            const std::string alt = field(line, 16, 1);
            if (!(alt == " " || alt.empty() || alt == "A")) continue;
            const std::string res = field(line, 17, 3);
            if (res == "HOH" || res == "DOD") continue;
        }
        if (!open_frame) {
            if (max_frames >= 0 && static_cast<int>(frames.size()) >= max_frames) break;
            frames.push_back(ProteinFrame());
            open_frame = true;
        }
        ProteinFrame& f = frames.back();
        std::string type = field(line, 12, 4);
        if (is_hetatm) {
            type = "HET:" + type;
        } else {
            type = internal::trimmed(type);
            if (type.empty()) type = "UNK";
        }
        std::string resname = internal::trimmed(field(line, 17, 3));
        if (resname.empty()) resname = "UNK";
        const int resseq = std::atoi(field(line, 22, 4).c_str());
        std::ostringstream pname;
        pname << "Atom " + type << " of residue " << resseq;
        // IMP parses the three coordinate fields into `float` before they
        // reach the particle, so every number downstream -- a rotamer placed
        // in a backbone frame, a Labelizer score -- was computed on values
        // rounded to single precision. Round the same way, on purpose.
        f.coords.push_back(static_cast<double>(static_cast<float>(std::atof(field(line, 30, 8).c_str()))));
        f.coords.push_back(static_cast<double>(static_cast<float>(std::atof(field(line, 38, 8).c_str()))));
        f.coords.push_back(static_cast<double>(static_cast<float>(std::atof(field(line, 46, 8).c_str()))));
        f.residue_indices.push_back(resseq);
        f.atom_names.push_back(name_from_particle_name(pname.str()));
        f.atom_types.push_back(type);
        f.resnames.push_back(resname);
        f.chain_ids.push_back(std::string(1, line.size() > 21 ? line[21] : ' '));
    }
    if (frames.empty()) {
        IMP_THROW("No molecule read from file " << path, IMP::ValueException);
    }
    return frames;
}

std::vector<ProteinFrame> load_protein_frames(const std::string& path,
                                              int max_frames) {
    if (!internal::ends_with(path, ".pdb") && !internal::ends_with(path, ".ent")) {
        IMP_THROW("not a PDB file: " << path
                                     << " (RMF is read through IMP.rmf and "
                                        "protein_frame_from_hierarchy)",
                  IMP::ValueException);
    }
    return internal::read_pdb_frames(path, internal::PDB_NON_WATER, false, max_frames);
}

IMPBFF_END_NAMESPACE
