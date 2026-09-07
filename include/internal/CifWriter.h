/**
 *  \file IMP/bff/internal/CifWriter.h
 *  \brief The one place this module turns values into mmCIF text.
 *
 * Reading CIF is IMP's job and is done through IMP's parser: `ihm_format.h`,
 * the C reader IMP vendors at `modules/core/dependency/python-ihm/src/` and
 * drives its own `IMP::atom::read_mmcif` with. This module adds no parser.
 *
 * Writing has no IMP facility to borrow. The vendored ihm library is
 * read-only, `IMP::atom` has `read_mmcif`/`read_bcif` and no writer, and
 * `IMP.mmcif` is Python (python-ihm) building deposition systems -- a C++
 * module cannot call it, and it has no spelling for the `_ff_*`, `_cgprobe_*`
 * or `_flr_*` tables written here. So the writer lives here, once: every
 * category this module emits goes through this class, and the quoting rule
 * below is stated in exactly one place.
 *
 * Internal. Nothing public takes a CifWriter.
 */
#ifndef IMPBFF_INTERNAL_CIFWRITER_H
#define IMPBFF_INTERNAL_CIFWRITER_H

#include <IMP/bff/bff_config.h>

#include <cctype>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

IMPBFF_BEGIN_INTERNAL_NAMESPACE

//! Writes mmCIF categories and loops to a stream.
class CifWriter {
public:
    explicit CifWriter(std::ostream& out) : out_(out) {}

    //! `data_<name>`, the block header.
    void start_block(const std::string& name) { out_ << "data_" << name << "\n"; }

    //! A key-value category: one `_name.key value` line per pair.
    void write_category(const std::string& name,
                        const std::vector<std::pair<std::string, std::string>>& kv) {
        out_ << "#" << "\n";
        for (auto& [k, v] : kv) {
            out_ << name << "." << k << " " << quote(v) << "\n";
        }
        out_ << "#" << "\n";
    }

    //! A `loop_` with \p cols columns and one line per row.
    void write_loop(const std::string& name,
                    const std::vector<std::string>& cols,
                    const std::vector<std::vector<std::string>>& rows) {
        out_ << "#" << "\n" << "loop_" << "\n";
        for (auto& c : cols) out_ << name << "." << c << "\n";
        for (auto& row : rows) {
            bool after_text_field = false;
            for (size_t i = 0; i < row.size(); i++) {
                const std::string v = quote(row[i]);
                // A text field opens and closes with a `;` at the start of a
                // line, so it brings its own leading newline and the value
                // after it has to start on a line of its own too.
                const bool text_field = v.compare(0, 2, "\n;") == 0;
                if (i) out_ << (after_text_field ? "\n" : (text_field ? "" : " "));
                out_ << v;
                after_text_field = text_field;
            }
            out_ << "\n";
        }
        out_ << "#" << "\n";
    }

    //! A value the ihm reader will tokenize as one token.
    /*!
        Quoting here, and not at each call site, is what makes a metadata value
        like `N1 and resname A48, N2 and resname A48` survive a round trip --
        the reader splits rows on whitespace, so an unquoted multi-word value
        is read as several and the loop comes up short. Already-quoted values
        and the CIF nulls (`.`/`?`) pass through untouched.

        The delimiter is chosen, not escaped. CIF has no escape: a quoted
        string ends at its delimiter *followed by whitespace or end of line*,
        which means an embedded quote is ordinary text unless whitespace
        follows it, and no doubling makes it otherwise. (Writing `""` for a
        literal quote -- the CSV convention -- reads back as two quotes.) So a
        value that cannot be closed with `"` is written with `'`, and one that
        can be closed with neither, or that contains a newline, is written as a
        semicolon text field, whose delimiters must each start a line.
    */
    static std::string quote(const std::string& s) {
        if (s.empty() || s == "." || s == "?") return s;
        if (s[0] == '"' || s[0] == '\'') return s;
        const bool plain = s.find_first_of(" \t\r\n,#'\"") == std::string::npos
                && s[0] != '_' && s[0] != ';';
        if (plain) return s;
        if (s.find('\n') == std::string::npos && s.find('\r') == std::string::npos) {
            if (closable_with(s, '"')) return "\"" + s + "\"";
            if (closable_with(s, '\'')) return "'" + s + "'";
        }
        // The one form with no forbidden character: the delimiters are lines
        // of their own, so nothing inside can end it early.
        return "\n;" + s + "\n;";
    }

private:
    //! Whether \p delim closes \p s: no `delim`+whitespace inside, and \p s
    //! does not itself end with one (which the closing delimiter would follow).
    static bool closable_with(const std::string& s, char delim) {
        if (s[s.size() - 1] == delim) return false;
        for (std::size_t i = 0; i + 1 < s.size(); ++i) {
            if (s[i] == delim && std::isspace(static_cast<unsigned char>(s[i + 1]))) {
                return false;
            }
        }
        return true;
    }

    std::ostream& out_;
};

//! An empty string is the CIF null; everything else is quoted by the rule above.
inline std::string cif_val(const std::string& s) {
    return s.empty() ? "." : CifWriter::quote(s);
}
inline std::string cif_val(const char* s) { return cif_val(std::string(s)); }
inline std::string cif_val(double v) {
    std::ostringstream ss;
    ss << v;
    return ss.str();
}
inline std::string cif_val(int v) { return std::to_string(v); }
inline std::string cif_val(bool v) { return v ? "YES" : "NO"; }

//! An empty string is the CIF null. Kept as its own name where a caller means
//! "omitted" rather than "empty"; the text is the same.
inline std::string cif_val_or_omit(const std::string& s) { return cif_val(s); }

IMPBFF_END_INTERNAL_NAMESPACE

#endif  // IMPBFF_INTERNAL_CIFWRITER_H
