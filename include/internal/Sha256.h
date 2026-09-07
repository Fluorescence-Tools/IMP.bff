/**
 *  \file IMP/bff/internal/Sha256.h
 *  \brief SHA-256, header-only, for container checksums.
 *
 * FIPS 180-4. Written here rather than pulled in because it is the only hash
 * this package needs and the alternative is a dependency on OpenSSL for sixty
 * lines of arithmetic. It is not a security primitive in this use: it answers
 * "are these the bytes that were written", which is what
 * `_mmfdb_artifact.checksum` records.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_INTERNAL_SHA256_H
#define IMPBFF_INTERNAL_SHA256_H

#include <IMP/bff/bff_config.h>

#include <cstddef>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

IMPBFF_BEGIN_INTERNAL_NAMESPACE

namespace sha256_detail {

inline unsigned int rotr(unsigned int x, unsigned int n) {
    return (x >> n) | (x << (32 - n));
}

inline const unsigned int* k_table() {
    static const unsigned int k[64] = {
        0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu,
        0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u,
        0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u,
        0xc19bf174u, 0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
        0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau, 0x983e5152u,
        0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
        0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu,
        0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
        0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u,
        0xd6990624u, 0xf40e3585u, 0x106aa070u, 0x19a4c116u, 0x1e376c08u,
        0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu,
        0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
        0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};
    return k;
}

}  // namespace sha256_detail

//! The SHA-256 of a byte range, lowercase hexadecimal.
inline std::string sha256_hex(const unsigned char* data, std::size_t n) {
    using namespace sha256_detail;
    unsigned int h[8] = {0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
                         0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u};

    // Message, padded to a multiple of 64 with the bit length at the end.
    std::vector<unsigned char> m(data, data + n);
    const unsigned long long bits = static_cast<unsigned long long>(n) * 8ull;
    m.push_back(0x80);
    while (m.size() % 64 != 56) m.push_back(0x00);
    for (int i = 7; i >= 0; --i) {
        m.push_back(static_cast<unsigned char>((bits >> (i * 8)) & 0xff));
    }

    const unsigned int* k = k_table();
    for (std::size_t off = 0; off < m.size(); off += 64) {
        unsigned int w[64];
        for (int i = 0; i < 16; ++i) {
            w[i] = (static_cast<unsigned int>(m[off + 4 * i + 0]) << 24) |
                   (static_cast<unsigned int>(m[off + 4 * i + 1]) << 16) |
                   (static_cast<unsigned int>(m[off + 4 * i + 2]) << 8) |
                   (static_cast<unsigned int>(m[off + 4 * i + 3]));
        }
        for (int i = 16; i < 64; ++i) {
            const unsigned int s0 =
                    rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const unsigned int s1 =
                    rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        unsigned int a = h[0], b = h[1], c = h[2], d = h[3];
        unsigned int e = h[4], f = h[5], g = h[6], hh = h[7];
        for (int i = 0; i < 64; ++i) {
            const unsigned int S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
            const unsigned int ch = (e & f) ^ ((~e) & g);
            const unsigned int t1 = hh + S1 + ch + k[i] + w[i];
            const unsigned int S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
            const unsigned int maj = (a & b) ^ (a & c) ^ (b & c);
            const unsigned int t2 = S0 + maj;
            hh = g; g = f; f = e; e = d + t1;
            d = c; c = b; b = a; a = t1 + t2;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d;
        h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
    }

    char out[65];
    for (int i = 0; i < 8; ++i) std::sprintf(out + 8 * i, "%08x", h[i]);
    out[64] = '\0';
    return std::string(out, 64);
}

//! The SHA-256 of a string, lowercase hexadecimal.
inline std::string sha256_hex(const std::string& s) {
    return sha256_hex(reinterpret_cast<const unsigned char*>(s.data()),
                      s.size());
}

IMPBFF_END_INTERNAL_NAMESPACE

#endif /* IMPBFF_INTERNAL_SHA256_H */
