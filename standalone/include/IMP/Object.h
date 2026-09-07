/**
 *  \file IMP/Object.h
 *  \brief IMP::Object for the standalone build: a named, reference-counted base.
 *
 * The subset of IMP's Object the core uses: a name (with IMP's "%1%"
 * pattern replaced by a running number), set_was_used(), show(), and the
 * intrusive reference count IMP::Pointer and the Python binding drive
 * through ref()/unref(). Nothing else of IMP's Object -- no Model, no
 * logging levels, no check levels.
 */
#ifndef IMPBFF_STANDALONE_OBJECT_H
#define IMPBFF_STANDALONE_OBJECT_H

#include <atomic>
#include <cstddef>
#include <iostream>
#include <string>

namespace IMP {

class Object {
    std::string name_;
    mutable std::atomic<int> refs_{0};
    bool was_used_ = false;

    static std::string expand(std::string name) {
        static std::atomic<unsigned long> counter{0};
        const std::size_t at = name.find("%1%");
        if (at != std::string::npos) {
            name.replace(at, 3, std::to_string(++counter));
        }
        return name;
    }

 public:
    explicit Object(const std::string& name) : name_(expand(name)) {}
    virtual ~Object() {}
    Object(const Object&) = delete;
    Object& operator=(const Object&) = delete;

    const std::string& get_name() const { return name_; }
    void set_name(const std::string& name) { name_ = name; }
    void set_was_used(bool tf) { was_used_ = tf; }
    bool get_was_used() const { return was_used_; }
    virtual std::string get_type_name() const { return "Object"; }
    virtual void show(std::ostream& out = std::cout) const { out << "\"" << get_name() << "\""; }

    //! Intrusive reference count -- what IMP::Pointer and SWIG's ref/unref features drive.
    void ref() const { ++refs_; }
    void unref() const {
        if (--refs_ == 0) delete this;
    }
    int get_ref_count() const { return refs_.load(); }
    //! Drop one count without the delete unref() does at zero: Pointer::release().
    void release_ref() const { --refs_; }
};

inline std::ostream& operator<<(std::ostream& out, const Object& o) {
    o.show(out);
    return out;
}

}  // namespace IMP

//! The methods IMP's IMP_OBJECT_METHODS declares, reduced to what a standalone object needs.
#define IMP_OBJECT_METHODS(Name)                                        \
    virtual std::string get_type_name() const override { return #Name; } \
    virtual ~Name() {}

#endif  // IMPBFF_STANDALONE_OBJECT_H
