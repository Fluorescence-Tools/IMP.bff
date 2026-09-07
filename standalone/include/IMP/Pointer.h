/**
 *  \file IMP/Pointer.h
 *  \brief IMP::Pointer for the standalone build: an intrusive smart pointer over IMP::Object.
 */
#ifndef IMPBFF_STANDALONE_POINTER_H
#define IMPBFF_STANDALONE_POINTER_H

#include <IMP/Object.h>
#include <cstddef>
#include <utility>

namespace IMP {

template <class O>
class Pointer {
    O* ptr_ = nullptr;
    void take(O* p) {
        if (p) p->ref();
        O* old = ptr_;
        ptr_ = p;
        if (old) old->unref();
    }

 public:
    Pointer() {}
    Pointer(O* p) { take(p); }  // NOLINT: implicit, as IMP's is
    Pointer(const Pointer& o) { take(o.ptr_); }
    template <class P>
    Pointer(const Pointer<P>& o) { take(o.get()); }
    Pointer& operator=(const Pointer& o) { take(o.ptr_); return *this; }
    Pointer& operator=(O* p) { take(p); return *this; }
    ~Pointer() { if (ptr_) ptr_->unref(); }

    O* get() const { return ptr_; }
    O* operator->() const { return ptr_; }
    O& operator*() const { return *ptr_; }
    operator O*() const { return ptr_; }  // NOLINT: as IMP's
    explicit operator bool() const { return ptr_ != nullptr; }
    bool operator!() const { return ptr_ == nullptr; }
    bool operator==(const Pointer& o) const { return ptr_ == o.ptr_; }
    bool operator!=(const Pointer& o) const { return ptr_ != o.ptr_; }
    bool operator==(const O* p) const { return ptr_ == p; }
    bool operator!=(const O* p) const { return ptr_ != p; }
    void reset() { take(nullptr); }
    //! Hand the object over without deleting it (IMP's semantics): the count
    //! drops, the object lives on for whoever takes the raw pointer.
    O* release() {
        O* p = ptr_;
        ptr_ = nullptr;
        if (p) p->release_ref();
        return p;
    }
};

template <class O>
using PointerMember = Pointer<O>;
template <class O>
using WeakPointer = O*;

}  // namespace IMP

//! IMP's IMP_NEW: a Pointer to a freshly constructed object.
#define IMP_NEW(Typename, varname, args) IMP::Pointer<Typename> varname(new Typename args)

#endif  // IMPBFF_STANDALONE_POINTER_H
