/**
 *  \file IMP/bff/Pto.h
 *  \brief PTO -- the EBML container this package's binary artefacts ride in.
 *
 *  A `.pto` file (PhoTon cOntainer, Thomas Peulen; the name winks at
 *  PicoQuant's `.ptu`) is an EBML document (RFC 8794) with `DocType "pto"`
 *  whose payloads are typed attached objects: each carries a name, a
 *  `PtoKind` saying what it is, a `PtoEncoding` saying how its bytes are
 *  coded, and the bytes themselves, 8-byte aligned so a mapped file slices
 *  them without copying. One grammar serves every domain -- photons and
 *  imaging as `.mmfbd.pto` (tttrlib), structures as `.chm.pto` (chimol),
 *  rotamer libraries as `.drot.pto` here. The domain lives in the file's
 *  profile suffix and in each object's kind, never in the container.
 *
 *  **This is a vendored core, not a dependency.** The rule the convergence
 *  record sets is that each repository writes the small walker from the
 *  specification rather than importing another's: imp.bff does not depend on
 *  tttrlib, and this file exists so that it does not have to start. It was
 *  written from `../tttrlib/okf/specs/pto-binary-decoding.md` (§2 the layout,
 *  §3 the reference C99 decoder); `../chimol/chimol/render/pto.py` is the
 *  same exercise in Python and the size a second implementation should be.
 *  Divergence is held down by a shared conformance corpus:
 *  `../tttrlib/test/tools/pto_ebml_check.cpp` walks a container with
 *  **libebml**, the reference implementation, and that -- never this file's
 *  own reader -- is what a container written here is proved against.
 *
 *  What this frames is bytes. What the payloads mean belongs to the profile
 *  and travels in the kind strings, which is the extension point PTO already
 *  provides: a new content pattern is a new kind, so no element ID is
 *  invented here and both profiles keep one grammar.
 *
 *  \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_PTO_H
#define IMPBFF_PTO_H

#include <IMP/bff/bff_config.h>

#include <cstddef>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! One attached object: what it is, and where its bytes are in the file.
struct IMPBFFEXPORT PtoObject {
    //! `FileUID` -- addresses the object independently of its position.
    unsigned long long uid;
    //! `FileName` -- what a reader asks for.
    std::string name;
    //! `PtoKind` -- what the payload is, e.g. `"drot.grid"`.
    std::string kind;
    //! `PtoEncoding` -- how the bytes are coded, e.g. `"f32.col+brotli"`.
    std::string encoding;
    //! Absolute offset of the payload in the file, and its length. The pair
    //! is what a mapped container slices and what an HTTP `Range` asks for.
    std::size_t offset, size;

    PtoObject() : uid(0), offset(0), size(0) {}
};

//! Write a PTO document: objects in, one file out.
/*!
    Objects are appended in the order they are added, each payload padded to
    the alignment boundary with `Void`. Sizes that close over the whole file
    (`Segment`, `Attachments`) are written as fixed-octet placeholders and
    patched at `close()` -- fixed-octet because a Data Size that narrows when
    its value shrinks shifts every byte after it, which is how a `FileUID`
    once corrupted a container; the lesson is cheap to keep and expensive to
    relearn.

    The file is complete only after `close()`, which the destructor calls.
*/
class IMPBFFEXPORT PtoWriter {
public:
    //! \param[in] path the file to write
    //! \param[in] doctype_version what the document claims to be
    explicit PtoWriter(const std::string& path, int doctype_version = 2);
    ~PtoWriter();

    //! Attach one payload; returns where it landed.
    /*! \param[in] name `FileName`
        \param[in] kind `PtoKind`
        \param[in] encoding `PtoEncoding`; empty means the bytes are plain
        \param[in] data,size the payload
        \throw IOException when the file cannot be written */
    PtoObject add(const std::string& name, const std::string& kind,
                  const std::string& encoding, const void* data,
                  std::size_t size);

    //! Patch the deferred sizes and close. Idempotent.
    void close();

    //! The objects written so far.
    const std::vector<PtoObject>& objects() const { return objects_; }

private:
    PtoWriter(const PtoWriter&);
    PtoWriter& operator=(const PtoWriter&);

    struct Impl;
    Impl* impl_;
    std::vector<PtoObject> objects_;
};

//! Read a PTO document: its objects, and any one payload on demand.
/*!
    Opening walks the framing only -- element headers, which are a few bytes
    each -- and reads no payload. `data()` then seeks to the one object asked
    for. That is what makes a container holding a hundred libraries cheap to
    open: pulling one 25 KB library out of a 19 MB file touches 25 KB, not
    19 MB. Seeking rather than mapping keeps this portable; the offsets in
    `PtoObject` are absolute, so a mapping reader can be dropped in later
    without changing a caller.
*/
class IMPBFFEXPORT PtoReader {
public:
    //! \throw IOException when the file is missing, is not a PTO document, or
    //!        its framing does not decode
    explicit PtoReader(const std::string& path);
    ~PtoReader();

    //! Every attached object, in file order.
    const std::vector<PtoObject>& objects() const { return objects_; }

    //! The object with this `FileName`, or -1.
    int find(const std::string& name) const;

    //! The payload bytes of an object of this reader's.
    std::vector<unsigned char> data(const PtoObject& object) const;

    //! What the document's `DocTypeVersion` says.
    int doctype_version() const { return doctype_version_; }

    //! Does this file begin with an EBML header?
    static bool looks_like_pto(const std::string& path);

private:
    PtoReader(const PtoReader&);
    PtoReader& operator=(const PtoReader&);

    struct Impl;
    Impl* impl_;
    std::vector<PtoObject> objects_;
    int doctype_version_;
};

IMPBFF_END_NAMESPACE

#endif /* IMPBFF_PTO_H */
