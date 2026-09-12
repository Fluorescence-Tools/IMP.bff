/**
 * \file DensityGrid.cpp
 * \brief The MRC writer for the module's lattice.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/DensityGrid.h>
#include <IMP/bff/IMPCompatibility.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

namespace {

void put_i32(std::ofstream& s, std::int32_t v) { s.write(reinterpret_cast<const char*>(&v), 4); }
void put_f32(std::ofstream& s, float v) { s.write(reinterpret_cast<const char*>(&v), 4); }

}  // namespace

void write_mrc(const std::string& path, const GridHeader& header,
               const float* values, std::size_t n) {
    const std::int32_t nx = header.get_nx(), ny = header.get_ny(), nz = header.get_nz();
    IMP_USAGE_CHECK(n == static_cast<std::size_t>(nx) * ny * nz,
                    "write_mrc: " << n << " values for a " << nx << "x" << ny << "x" << nz << " grid");
    std::ofstream s(path.c_str(), std::ios::binary);
    if (!s) IMP_THROW("Cannot write " << path, IOException);
    // the statistics the header carries, from the values themselves
    double dmin = n ? values[0] : 0.0, dmax = n ? values[0] : 0.0, sum = 0.0, sum2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double v = values[i];
        if (v < dmin) dmin = v;
        if (v > dmax) dmax = v;
        sum += v;
        sum2 += v * v;
    }
    const double mean = n ? sum / n : 0.0;
    const double var = n ? sum2 / n - mean * mean : 0.0;
    const float rms = static_cast<float>(var > 0.0 ? std::sqrt(var) : 0.0);
    const float spacing = header.get_spacing();
    // MRC2014, word by word (all 32-bit)
    put_i32(s, nx); put_i32(s, ny); put_i32(s, nz);       // 1-3   NX NY NZ
    put_i32(s, 2);                                        // 4     MODE 2: 32-bit float
    put_i32(s, 0); put_i32(s, 0); put_i32(s, 0);          // 5-7   NXSTART NYSTART NZSTART
    put_i32(s, nx); put_i32(s, ny); put_i32(s, nz);       // 8-10  MX MY MZ
    put_f32(s, nx * spacing); put_f32(s, ny * spacing); put_f32(s, nz * spacing);  // 11-13 CELLA
    put_f32(s, 90.0f); put_f32(s, 90.0f); put_f32(s, 90.0f);                      // 14-16 CELLB
    put_i32(s, 1); put_i32(s, 2); put_i32(s, 3);          // 17-19 MAPC MAPR MAPS
    put_f32(s, static_cast<float>(dmin)); put_f32(s, static_cast<float>(dmax));
    put_f32(s, static_cast<float>(mean));                 // 20-22 DMIN DMAX DMEAN
    put_i32(s, 1);                                        // 23    ISPG: 1, a volume
    put_i32(s, 0);                                        // 24    NSYMBT
    for (int i = 25; i <= 27; ++i) put_i32(s, 0);         // 25-27 EXTRA (27 = EXTTYP, none)
    put_i32(s, 20140);                                    // 28    NVERSION
    for (int i = 29; i <= 49; ++i) put_i32(s, 0);         // 29-49 EXTRA
    put_f32(s, header.get_xorigin()); put_f32(s, header.get_yorigin()); put_f32(s, header.get_zorigin());  // 50-52 ORIGIN
    s.write("MAP ", 4);                                   // 53    MAP
    {                                                     // 54    MACHST: this machine's byte order
        const std::uint16_t probe = 0x0102;
        const bool little = *reinterpret_cast<const unsigned char*>(&probe) == 0x02;
        const unsigned char stamp[4] = {static_cast<unsigned char>(little ? 0x44 : 0x11),
                                        static_cast<unsigned char>(little ? 0x44 : 0x11), 0, 0};
        s.write(reinterpret_cast<const char*>(stamp), 4);
    }
    put_f32(s, rms);                                      // 55    RMS
    put_i32(s, 1);                                        // 56    NLABL
    char labels[10 * 80];
    std::memset(labels, 0, sizeof(labels));
    std::memcpy(labels, "IMP.bff write_mrc", 17);        // 57-  LABEL 1..10
    s.write(labels, sizeof(labels));
    s.write(reinterpret_cast<const char*>(values), static_cast<std::streamsize>(sizeof(float) * n));
    if (!s) IMP_THROW("Error writing " << path, IOException);
}

void DensityGrid::write_mrc(const std::string& path) const {
    std::vector<float> values(data_.begin(), data_.end());
    IMP::bff::write_mrc(path, header_, values.data(), values.size());
}

IMPBFF_END_NAMESPACE
