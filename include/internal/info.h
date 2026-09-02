// SPDX-License-Identifier: BSD-3-Clause
//
// A deliberately reduced copy of tttrlib's `modules/util/include/info.h`,
// carrying only the three environment-variable helpers that the vendored
// `Random.h` (in this same directory) needs via its unqualified
// `#include "info.h"`. The full file also drives tttrlib's AVX/NEON runtime
// dispatch, which depends on tttrlib's own CMake-defined macros
// (TTTRLIB_WITH_AVX, TTTRLIB_COMPILE_AVX, ...) and does not compile outside
// that build -- so unlike `Random.h` and `DecayConvolution.h`, this one is
// NOT pinned byte-identical against the upstream file. The three functions
// below ARE copied verbatim (same names, same bodies, same namespace) so
// that a diff against tttrlib's info.h for just those symbols is still
// meaningful.
//
// tttrlib owns the originals, at `modules/util/include/info.h` lines
// ~136-166 (`safe_getenv`, `is_false_value`, `is_feature_enabled_by_env`).
#ifndef TTTRLIB_INFO_H
#define TTTRLIB_INFO_H

#include <cstdlib>
#include <cstring>

namespace tttrlib {
namespace cpu_features {

    // Cross-platform helper to safely get environment variable
    inline const char* safe_getenv(const char* env_var_name, char* buffer, size_t buffer_size) {
#ifdef _WIN32
        size_t required_size = 0;
        errno_t err = getenv_s(&required_size, buffer, buffer_size, env_var_name);
        return (err == 0 && required_size > 0) ? buffer : nullptr;
#else
        (void)buffer;       // Unused on non-Windows
        (void)buffer_size;  // Unused on non-Windows
        return std::getenv(env_var_name);
#endif
    }

    // Helper to check if string represents a false value
    inline bool is_false_value(const char* value) {
        return std::strcmp(value, "0") == 0 ||
               std::strcmp(value, "false") == 0 ||
               std::strcmp(value, "FALSE") == 0 ||
               std::strcmp(value, "off") == 0 ||
               std::strcmp(value, "OFF") == 0;
    }

    // Helper to check environment variable for feature override
    inline bool is_feature_enabled_by_env(const char* env_var_name, bool default_value) {
        char env_buffer[256];
        const char* env_val = safe_getenv(env_var_name, env_buffer, sizeof(env_buffer));

        if (env_val != nullptr) {
            return !is_false_value(env_val);
        }
        return default_value;
    }

} // namespace cpu_features
} // namespace tttrlib

#endif // TTTRLIB_INFO_H
