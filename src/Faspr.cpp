/**
 * \file Faspr.cpp
 * \brief The IMP-facing entry point of the vendored FASPR port.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 * The pipeline is FASPR's main() (junk/FASPR/src/FASPR.cpp) with its CLI
 * removed: read PDB -> load the Dunbrack-2010 library -> build side-chain
 * rotamers -> self and pair energies -> DEE + tree-decomposition search ->
 * write the repacked PDB. FASPR's stdout log is silenced unless \c verbose
 * (the vendored code prints verbatim, so the silencing happens by swapping
 * std::cout's streambuf around the call -- nothing inside the port is
 * touched).
 */
#include <IMP/bff/RotamerLibrary.h>
#include <IMP/bff/Pto.h>

#include "FasprSearch.h"

#include <IMP/bff/Base.h>

#include <fstream>
#include <unistd.h>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

// The two globals FASPR's main() used to define (junk/FASPR/src/FASPR.cpp);
// RotamerBuilder.cpp declares them extern inside the faspr namespace.
namespace faspr {
string PROGRAM_PATH = ".";
string ROTLIB2010 = "dun2010bbdep.bin";
}  // namespace faspr

namespace {

/** A scoped std::cout silencer: FASPR logs to cout verbatim. */
class CoutSilencer {
 public:
  explicit CoutSilencer(bool verbose) : old_(nullptr) {
    if (!verbose) {
      old_ = std::cout.rdbuf(devnull_.rdbuf());
    }
  }
  ~CoutSilencer() { restore(); }
  void restore() {
    if (old_ != nullptr) {
      std::cout.rdbuf(old_);
      old_ = nullptr;
    }
  }

 private:
  std::streambuf* old_;
  std::ostringstream devnull_;
};

//! A `dun2010bbdep.bin` written beside the container, removed on the way out.
/*! FASPR insists on that name -- it is hard-coded -- so the scratch copy
    carries it, in a directory of its own so two calls cannot collide. */
struct TemporaryBin {
  std::string dir, path;

  TemporaryBin() {}
  ~TemporaryBin() {
    if (path.empty()) return;
    std::remove(path.c_str());
    rmdir(dir.c_str());
  }

  void make(const std::string& container) {
    char pattern[] = "/tmp/imp_bff_faspr_XXXXXX";
    const char* made = mkdtemp(pattern);
    if (made == NULL) {
      IMP_THROW("faspr_pack: cannot make a scratch directory for "
                        << container, IOException);
    }
    dir = made;
    path = dir + "/dun2010bbdep.bin";
    write_dunbrack_bin(container, path);
  }

 private:
  TemporaryBin(const TemporaryBin&);
  TemporaryBin& operator=(const TemporaryBin&);
};

}  // namespace

void faspr_pack(const std::string& pdb_in, const std::string& pdb_out,
                const std::string& rotamer_library, bool verbose) {
  if (pdb_in.empty() || pdb_out.empty()) {
    IMP_THROW("faspr_pack: pdb_in and pdb_out must be non-empty",
              IOException);
  }
  {
    std::ifstream probe(rotamer_library.c_str(), std::ios::in | std::ios::binary);
    if (!probe) {
      IMP_THROW("faspr_pack: cannot open rotamer library " + rotamer_library,
                IOException);
    }
  }

  // The shipped side-chain library is `sidechains.drot.pto`, the same
  // container the dye and spin-label rotamers ride in. The vendored engine
  // seeks inside `dun2010bbdep.bin` by name and byte offset, and that is
  // exactly what it should keep doing -- the parity pin depends on it -- so
  // a container is unpacked to a scratch copy first. The conversion is
  // byte-reversible (`write_dunbrack_bin`), so the engine sees the file it
  // was written against, bit for bit.
  std::string library = rotamer_library;
  TemporaryBin scratch;
  if (PtoReader::looks_like_pto(rotamer_library)) {
    scratch.make(rotamer_library);
    library = scratch.path;
  }

  // FASPR opens PROGRAM_PATH + "/" + ROTLIB2010; point them at the library.
  namespace fp = IMP::bff::faspr;
  std::string dir = ".";
  std::string file = library;
  const std::string::size_type slash = library.find_last_of("/\\");
  if (slash != std::string::npos) {
    dir = library.substr(0, slash);
    file = library.substr(slash + 1);
  }
  fp::PROGRAM_PATH = dir;
  fp::ROTLIB2010 = file;

  CoutSilencer quiet(verbose);
  try {
    // FASPR's IO takes non-const string& -- hand it mutable copies.
    std::string in = pdb_in;
    std::string out = pdb_out;
    fp::Solution faspr;
    faspr.ReadPDB(in);
    faspr.LoadSeq();
    faspr.BuildSidechain();
    faspr.CalcSelfEnergy();
    faspr.CalcPairEnergy();
    faspr.Search();
    faspr.WritePDB(out);
  } catch (...) {
    quiet.restore();
    throw;
  }
}

IMPBFF_END_NAMESPACE