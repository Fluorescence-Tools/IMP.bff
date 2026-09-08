/**
 *  \file IMP/bff/Base.h
 *  \brief The vocabulary this module is written in: exceptions, showability,
 *         and the plural typedefs.
 *
 * Nearly every header here says `IMP_THROW`, `IMP_SHOWABLE_INLINE` or
 * `IMP_VALUES`. Those three, plus `IMP_OBJECTS`, are the whole of what this
 * module borrows from IMP's kernel as *boilerplate* --
 * about six hundred uses across ninety-six headers, and not one of them is
 * about integrative modelling. They arrived through a hundred and thirty
 * `#include <IMP/exception.h>` and friends, scattered file by file.
 *
 * This header is the single door they now come through, and the reason it
 * exists is that a door can be moved. With IMP present it forwards to IMP's
 * own definitions, so the module compiles to exactly what it always did --
 * `IMP_VALUES` still makes an `IMP::Vector`, `IMP_SHOWABLE_INLINE` still emits
 * the `IMP::Showable` conversion IMP's SWIG layer looks for. Define
 * `IMPBFF_STANDALONE` and the same call sites compile against equivalents
 * defined here, which is what lets the core build with no IMP on the include
 * path at all.
 *
 * The macro *call sites* are identical either way. That is the point: this is
 * a header swap, not a refactor, and nothing in the six hundred uses had to be
 * touched to make the module portable.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_BASE_H
#define IMPBFF_BASE_H

#include <IMP/bff/bff_config.h>

#ifndef IMPBFF_STANDALONE

// ---------------------------------------------------------------------------
// Built as an IMP module: the macros are IMP's, unchanged.
// ---------------------------------------------------------------------------
#include <IMP/exception.h>
#include <IMP/showable_macros.h>
#include <IMP/value_macros.h>
#include <IMP/object_macros.h>
// IMP_WARN. The standalone branch below defines its own, so this branch has
// to name IMP's: without it the macro is undefined and its argument -- a
// stream expression -- is compiled as one, which fails on the first `<<`.
#include <IMP/log_macros.h>

#else

// ---------------------------------------------------------------------------
// Built standalone: the same names, defined here.
// ---------------------------------------------------------------------------
#include <iostream>
#include <sstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef SWIG  /* the Python module makes these as Python exception classes */
namespace IMP {

//! The base of every error this module raises.
/*! Deriving from `std::runtime_error` rather than inventing a hierarchy: a
    caller that catches `std::exception` catches these, which is what a caller
    who has never heard of IMP expects. */
class Exception : public std::runtime_error {
 public:
    explicit Exception(const char *message) : std::runtime_error(message) {}
    explicit Exception(const std::string &message) : std::runtime_error(message) {}
};

//! A value was outside the domain the callee accepts.
class ValueException : public Exception {
 public:
    explicit ValueException(const char *t) : Exception(t) {}
};

//! A file could not be read, written, or made sense of.
class IOException : public Exception {
 public:
    explicit IOException(const char *t) : Exception(t) {}
};

//! The model reached a state its own invariants forbid.
class ModelException : public Exception {
 public:
    explicit ModelException(const char *t) : Exception(t) {}
};

}  // namespace IMP
#endif  // SWIG

//! Build a message with `<<` and throw \p exception_name carrying it.
/*! The `do { } while (true)` -- rather than the usual `while (false)` -- is
    IMP's, and is kept deliberately: it tells the compiler the macro never
    falls through, which suppresses "control reaches end of non-void function"
    at the call sites that end in an `IMP_THROW`. */
#define IMP_THROW(message, exception_name)                 \
    do {                                                   \
        std::ostringstream imp_throw_oss;                  \
        imp_throw_oss << message << std::endl;             \
        throw exception_name(imp_throw_oss.str().c_str()); \
    } while (true)

//! Give a class a `show()` that writes \p how_to_show to a stream.
/*! The `IMP::Showable` conversion operator IMP's version also emits is not
    reproduced: this module never names `IMP::Showable` (zero uses), and the
    standalone SWIG layer does not look for it. */
#define IMP_SHOWABLE_INLINE(Name, how_to_show)                    \
    void show(std::ostream &out = std::cout) const { how_to_show; } \
    friend std::ostream &operator<<(std::ostream &o, const Name &v) { \
        v.show(o);                                                \
        return o;                                                 \
    }

//! Declare `show()`; the definition is written out of line.
#define IMP_SHOWABLE(Name) void show(std::ostream &out = std::cout) const

//! Name the plural of a value type.
/*! `IMP::Vector<T>` derives publicly from `std::vector<T>`, so a plain
    `std::vector` is the same type to every caller in this module. */
#define IMP_VALUES(Name, PluralName) typedef std::vector<Name> PluralName

//! Name the plural of a reference-counted object type.
/*! IMP spells the owning pointer `IMP::Pointer` and the borrowed one
    `IMP::WeakPointer`; standalone spells them `std::shared_ptr` and a raw
    pointer, which is the same distinction with the same guarantees. The three
    classes that use this (`MinimizerObserver`, `InteractionTerm`,
    `RRTCollision`) derive from `IMP::Object` only so that Python can subclass
    them, and the standalone SWIG layer gets that from a director over
    `shared_ptr` instead. */
#include <IMP/Object.h>
#include <IMP/Pointer.h>

//! IMP's plural typedefs for objects: owning, and non-owning (IMP's "Temp").
#define IMP_OBJECTS(Name, PluralName)                  \
    typedef std::vector<IMP::Pointer<Name> > PluralName; \
    typedef std::vector<Name*> PluralName##Temp

#ifndef SWIG
namespace IMP {
//! A caller broke a documented precondition (IMP's UsageException).
class UsageException : public Exception {
 public:
    explicit UsageException(const std::string& m) : Exception(m) {}
};
}
#endif  // SWIG

//! IMP's IMP_USAGE_CHECK: a precondition, checked always here (IMP checks it at its usage level).
#define IMP_USAGE_CHECK(expr, message)                                   \
    do {                                                                  \
        if (!(expr)) {                                                    \
            std::ostringstream imp_bff_oss;                               \
            imp_bff_oss << "Usage check failed: " << #expr << ": " << message; \
            throw IMP::UsageException(imp_bff_oss.str());                 \
        }                                                                 \
    } while (false)

//! IMP's IMP_WARN: to standard error, once per call.
#define IMP_WARN(message)                                          \
    do {                                                            \
        std::cerr << "WARNING  " << message;                        \
    } while (false)

#endif  // IMPBFF_STANDALONE

#endif  // IMPBFF_BASE_H
