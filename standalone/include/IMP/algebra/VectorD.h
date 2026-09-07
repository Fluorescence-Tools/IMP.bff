/**
 *  \file IMP/algebra/VectorD.h
 *  \brief IMP::algebra::VectorD for the standalone build: a fixed-size vector of doubles.
 *
 * The subset the core uses -- indexing, the arithmetic operators, the
 * magnitude, the scalar and vector products -- with IMP's names, so a core
 * source compiles against either. Nothing of IMP's geometry beyond that.
 */
#ifndef IMPBFF_STANDALONE_ALGEBRA_VECTORD_H
#define IMPBFF_STANDALONE_ALGEBRA_VECTORD_H

#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace IMP {
namespace algebra {

template <int D>
class VectorD {
    double v_[D];

 public:
    VectorD() { for (int i = 0; i < D; ++i) v_[i] = 0.0; }
    VectorD(double x, double y, double z) {
        static_assert(D == 3, "three values make a Vector3D");
        v_[0] = x; v_[1] = y; v_[2] = z;
    }
    VectorD(double x, double y, double z, double w) {
        static_assert(D == 4, "four values make a Vector4D");
        v_[0] = x; v_[1] = y; v_[2] = z; v_[3] = w;
    }
    explicit VectorD(const std::vector<double>& v) {
        if (static_cast<int>(v.size()) != D) throw std::invalid_argument("VectorD: wrong length");
        for (int i = 0; i < D; ++i) v_[i] = v[i];
    }
    static unsigned int get_dimension() { return D; }
    double operator[](unsigned int i) const { return v_[i]; }
    double& operator[](unsigned int i) { return v_[i]; }
    const double* get_data() const { return v_; }
    double get_squared_magnitude() const {
        double s = 0.0;
        for (int i = 0; i < D; ++i) s += v_[i] * v_[i];
        return s;
    }
    double get_magnitude() const { return std::sqrt(get_squared_magnitude()); }
    VectorD get_unit_vector() const {
        const double m = get_magnitude();
        return m > 0.0 ? *this / m : *this;
    }
    double get_scalar_product(const VectorD& o) const {
        double s = 0.0;
        for (int i = 0; i < D; ++i) s += v_[i] * o.v_[i];
        return s;
    }
    VectorD operator+(const VectorD& o) const { VectorD r; for (int i = 0; i < D; ++i) r.v_[i] = v_[i] + o.v_[i]; return r; }
    VectorD operator-(const VectorD& o) const { VectorD r; for (int i = 0; i < D; ++i) r.v_[i] = v_[i] - o.v_[i]; return r; }
    VectorD operator-() const { VectorD r; for (int i = 0; i < D; ++i) r.v_[i] = -v_[i]; return r; }
    VectorD operator*(double s) const { VectorD r; for (int i = 0; i < D; ++i) r.v_[i] = v_[i] * s; return r; }
    VectorD operator/(double s) const { VectorD r; for (int i = 0; i < D; ++i) r.v_[i] = v_[i] / s; return r; }
    double operator*(const VectorD& o) const { return get_scalar_product(o); }
    VectorD& operator+=(const VectorD& o) { for (int i = 0; i < D; ++i) v_[i] += o.v_[i]; return *this; }
    VectorD& operator-=(const VectorD& o) { for (int i = 0; i < D; ++i) v_[i] -= o.v_[i]; return *this; }
    VectorD& operator*=(double s) { for (int i = 0; i < D; ++i) v_[i] *= s; return *this; }
    VectorD& operator/=(double s) { for (int i = 0; i < D; ++i) v_[i] /= s; return *this; }
    void show(std::ostream& out = std::cout) const {
        out << "(";
        for (int i = 0; i < D; ++i) out << (i ? ", " : "") << v_[i];
        out << ")";
    }
};

template <int D>
inline VectorD<D> operator*(double s, const VectorD<D>& v) { return v * s; }
template <int D>
inline std::ostream& operator<<(std::ostream& out, const VectorD<D>& v) { v.show(out); return out; }
template <int D>
inline double get_squared_distance(const VectorD<D>& a, const VectorD<D>& b) { return (a - b).get_squared_magnitude(); }
template <int D>
inline double get_distance(const VectorD<D>& a, const VectorD<D>& b) { return (a - b).get_magnitude(); }
template <int D>
inline VectorD<D> get_zero_vector_d() { return VectorD<D>(); }

typedef VectorD<3> Vector3D;
typedef VectorD<4> Vector4D;
typedef std::vector<Vector3D> Vector3Ds;
typedef std::vector<Vector4D> Vector4Ds;

inline Vector3D get_vector_product(const Vector3D& a, const Vector3D& b) {
    return Vector3D(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

}  // namespace algebra
}  // namespace IMP

#endif  // IMPBFF_STANDALONE_ALGEBRA_VECTORD_H
