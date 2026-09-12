/**
 *  \file SelectionExpression.cpp
 *  \brief The selection language: tokenize, parse, evaluate, compile.
 */

#include <IMP/bff/SelectionExpression.h>

#include <IMP/bff/internal/Text.h>

#include <IMP/bff/Base.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <map>
#include <set>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

using internal::trimmed;
using internal::upper;

namespace {

// --------------------------------------------------------------------------
// The keyword table
// --------------------------------------------------------------------------

//! What a keyword does.
enum KeywordKind {
    KW_PROPERTY,   //!< takes one value: `name CA`
    KW_FLAG,       //!< takes none: `hydro`
    KW_UNARY,      //!< `not`, `byres`
    KW_BINARY,     //!< `and`, `or`, `-`
    KW_DISTANCE,   //!< `within R of SEL`, `beyond`, `around`
    KW_SAME,       //!< `same residue as SEL`
    KW_UNSUPPORTED //!< valid notation, but not answerable from a structure
};

//! Which atom field a property keyword reads.
enum Property {
    P_NAME, P_RESN, P_RESI, P_CHAIN, P_SEGI, P_ELEM, P_ALT, P_INDEX, P_ID
};

//! Which class of atom a flag keyword selects.
enum Flag {
    F_ALL, F_NONE, F_HYDRO, F_NOH, F_HETATM, F_POLYMER, F_PROTEIN, F_NUCLEIC,
    F_SOLVENT, F_ORGANIC, F_INORGANIC, F_BACKBONE, F_SIDECHAIN, F_GUIDE
};

//! The distance operators.
enum Distance { D_WITHIN, D_BEYOND, D_AROUND };

//! Operators, with their precedences: `not` 0x70, `and` 0x60, `or` 0x40,
//! `byres` 0x20 -- so `byres` binds loosest and `not` tightest.
enum Op { OP_NOT, OP_AND, OP_OR, OP_ANDNOT, OP_BYRES };

struct Keyword {
    KeywordKind kind;
    int what;      //!< a Property, a Flag or an Op
    int precedence;
};

const std::map<std::string, Keyword>& keyword_table() {
    static const std::map<std::string, Keyword> t = {
        // operators
        {"not", {KW_UNARY, OP_NOT, 0x70}},  {"!", {KW_UNARY, OP_NOT, 0x70}},
        {"and", {KW_BINARY, OP_AND, 0x60}}, {"&", {KW_BINARY, OP_AND, 0x60}},
        {"or", {KW_BINARY, OP_OR, 0x40}},   {"|", {KW_BINARY, OP_OR, 0x40}},
        {"+", {KW_BINARY, OP_OR, 0x40}},    {"-", {KW_BINARY, OP_ANDNOT, 0x40}},
        {"byres", {KW_UNARY, OP_BYRES, 0x20}},
        {"byresi", {KW_UNARY, OP_BYRES, 0x20}},
        {"byresidue", {KW_UNARY, OP_BYRES, 0x20}},
        {"br.", {KW_UNARY, OP_BYRES, 0x20}},
        // properties
        {"name", {KW_PROPERTY, P_NAME, 0x80}}, {"n.", {KW_PROPERTY, P_NAME, 0x80}},
        {"resn", {KW_PROPERTY, P_RESN, 0x80}},
        {"resname", {KW_PROPERTY, P_RESN, 0x80}},
        {"r.", {KW_PROPERTY, P_RESN, 0x80}},
        {"resi", {KW_PROPERTY, P_RESI, 0x80}},
        {"resid", {KW_PROPERTY, P_RESI, 0x80}},
        {"residue", {KW_PROPERTY, P_RESI, 0x80}},
        {"i.", {KW_PROPERTY, P_RESI, 0x80}},
        {"chain", {KW_PROPERTY, P_CHAIN, 0x80}},
        {"c.", {KW_PROPERTY, P_CHAIN, 0x80}},
        {"segi", {KW_PROPERTY, P_SEGI, 0x80}},
        {"segname", {KW_PROPERTY, P_SEGI, 0x80}},
        {"segid", {KW_PROPERTY, P_SEGI, 0x80}},
        {"segment", {KW_PROPERTY, P_SEGI, 0x80}},
        {"s.", {KW_PROPERTY, P_SEGI, 0x80}},
        {"elem", {KW_PROPERTY, P_ELEM, 0x80}},
        {"element", {KW_PROPERTY, P_ELEM, 0x80}},
        {"symbol", {KW_PROPERTY, P_ELEM, 0x80}},
        {"e.", {KW_PROPERTY, P_ELEM, 0x80}},
        {"alt", {KW_PROPERTY, P_ALT, 0x80}},
        {"altloc", {KW_PROPERTY, P_ALT, 0x80}},
        {"index", {KW_PROPERTY, P_INDEX, 0x80}},
        {"idx.", {KW_PROPERTY, P_INDEX, 0x80}},
        {"id", {KW_PROPERTY, P_ID, 0x80}},
        {"serial", {KW_PROPERTY, P_ID, 0x80}},
        // flags
        {"all", {KW_FLAG, F_ALL, 0x90}},   {"*", {KW_FLAG, F_ALL, 0x90}},
        {"none", {KW_FLAG, F_NONE, 0x90}},
        {"hydro", {KW_FLAG, F_HYDRO, 0x90}},
        {"hydrogens", {KW_FLAG, F_HYDRO, 0x90}},
        {"h.", {KW_FLAG, F_HYDRO, 0x90}},
        {"hetatm", {KW_FLAG, F_HETATM, 0x80}},
        {"het", {KW_FLAG, F_HETATM, 0x80}},
        {"polymer", {KW_FLAG, F_POLYMER, 0x90}},
        {"pol.", {KW_FLAG, F_POLYMER, 0x90}},
        {"protein", {KW_FLAG, F_PROTEIN, 0x90}},
        {"polymer.protein", {KW_FLAG, F_PROTEIN, 0x90}},
        {"pro.", {KW_FLAG, F_PROTEIN, 0x90}},
        {"nucleic", {KW_FLAG, F_NUCLEIC, 0x90}},
        {"polymer.nucleic", {KW_FLAG, F_NUCLEIC, 0x90}},
        {"nuc.", {KW_FLAG, F_NUCLEIC, 0x90}},
        {"solvent", {KW_FLAG, F_SOLVENT, 0x90}},
        {"water", {KW_FLAG, F_SOLVENT, 0x90}},
        {"waters", {KW_FLAG, F_SOLVENT, 0x90}},
        {"noh", {KW_FLAG, F_NOH, 0x90}},
        {"sol.", {KW_FLAG, F_SOLVENT, 0x90}},
        {"organic", {KW_FLAG, F_ORGANIC, 0x90}},
        {"org.", {KW_FLAG, F_ORGANIC, 0x90}},
        {"inorganic", {KW_FLAG, F_INORGANIC, 0x90}},
        {"ino.", {KW_FLAG, F_INORGANIC, 0x90}},
        {"backbone", {KW_FLAG, F_BACKBONE, 0x90}},
        {"bb.", {KW_FLAG, F_BACKBONE, 0x90}},
        {"sidechain", {KW_FLAG, F_SIDECHAIN, 0x90}},
        {"sc.", {KW_FLAG, F_SIDECHAIN, 0x90}},
        {"guide", {KW_FLAG, F_GUIDE, 0x90}},
        // valid notation this layer cannot answer
        {"within", {KW_DISTANCE, D_WITHIN, 0x30}},
        {"w.", {KW_DISTANCE, D_WITHIN, 0x30}},
        {"beyond", {KW_DISTANCE, D_BEYOND, 0x30}},
        {"be.", {KW_DISTANCE, D_BEYOND, 0x30}},
        {"around", {KW_DISTANCE, D_AROUND, 0x30}},
        {"a.", {KW_DISTANCE, D_AROUND, 0x30}},
        {"same", {KW_SAME, 0, 0x20}},
        {"expand", {KW_UNSUPPORTED, 0, 0}}, {"x.", {KW_UNSUPPORTED, 0, 0}},
        {"extend", {KW_UNSUPPORTED, 0, 0}}, {"xt.", {KW_UNSUPPORTED, 0, 0}},
        {"near_to", {KW_UNSUPPORTED, 0, 0}}, {"nto.", {KW_UNSUPPORTED, 0, 0}},
        {"neighbor", {KW_UNSUPPORTED, 0, 0}}, {"nbr.", {KW_UNSUPPORTED, 0, 0}},
        {"bound_to", {KW_UNSUPPORTED, 0, 0}}, {"bto.", {KW_UNSUPPORTED, 0, 0}},
        {"bymolecule", {KW_UNSUPPORTED, 0, 0}}, {"bymol", {KW_UNSUPPORTED, 0, 0}},
        {"byobject", {KW_UNSUPPORTED, 0, 0}}, {"byobj", {KW_UNSUPPORTED, 0, 0}},
        {"model", {KW_UNSUPPORTED, 0, 0}}, {"m.", {KW_UNSUPPORTED, 0, 0}},
        {"state", {KW_UNSUPPORTED, 0, 0}}, {"rep", {KW_UNSUPPORTED, 0, 0}},
        {"color", {KW_UNSUPPORTED, 0, 0}}, {"visible", {KW_UNSUPPORTED, 0, 0}},
        {"enabled", {KW_UNSUPPORTED, 0, 0}}, {"flag", {KW_UNSUPPORTED, 0, 0}},
        {"b", {KW_UNSUPPORTED, 0, 0}}, {"q", {KW_UNSUPPORTED, 0, 0}},
        {"partial_charge", {KW_UNSUPPORTED, 0, 0}},
        {"formal_charge", {KW_UNSUPPORTED, 0, 0}},
        {"ss", {KW_UNSUPPORTED, 0, 0}}, {"in", {KW_UNSUPPORTED, 0, 0}},
        {"like", {KW_UNSUPPORTED, 0, 0}}, {"pepseq", {KW_UNSUPPORTED, 0, 0}},
    };
    return t;
}

std::string lower(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::tolower(
                              static_cast<unsigned char>(c)));
    return s;
}

// --------------------------------------------------------------------------
// Residue classification, read off the residue name -- the only thing a PDB
// or mmCIF guarantees. A viewer that sets per-atom flags at load time can be
// more precise about, say, a modified residue.
// --------------------------------------------------------------------------

const std::set<std::string>& amino_acids() {
    static const std::set<std::string> s = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "MSE", "SEC", "PYL", "ASX", "GLX", "UNK", "HID", "HIE", "HIP", "CYX",
        "HSD", "HSE", "HSP", "ACE", "NME", "NMA"};
    return s;
}

const std::set<std::string>& nucleotides() {
    static const std::set<std::string> s = {
        "A", "C", "G", "T", "U", "DA", "DC", "DG", "DT", "DU",
        "ADE", "CYT", "GUA", "THY", "URA", "RA", "RC", "RG", "RU"};
    return s;
}

const std::set<std::string>& solvents() {
    static const std::set<std::string> s = {
        "HOH", "WAT", "H2O", "DOD", "SOL", "TIP", "TIP3", "TIP4", "SPC"};
    return s;
}

//! The backbone, in the order a chain is read: N, CA, C, O, and a terminal
//! OXT. The order is part of it -- a mask written from this list reads
//! `not name N+CA+C+O+OXT`, which is how a person writes it.
const std::vector<std::string>& protein_backbone_order() {
    static const std::vector<std::string> v = {"N", "CA", "C", "O", "OXT"};
    return v;
}

const std::set<std::string>& protein_backbone() {
    static const std::set<std::string> s(protein_backbone_order().begin(),
                                         protein_backbone_order().end());
    return s;
}

const std::set<std::string>& nucleic_backbone() {
    static const std::set<std::string> s = {
        "P", "OP1", "OP2", "O5'", "C5'", "C4'", "C3'", "O3'", "O1P", "O2P"};
    return s;
}

bool is_protein(const SelectionAtom& a) {
    return amino_acids().count(upper(a.resn)) > 0;
}
bool is_nucleic(const SelectionAtom& a) {
    return nucleotides().count(upper(a.resn)) > 0;
}
bool is_solvent(const SelectionAtom& a) {
    return solvents().count(upper(a.resn)) > 0;
}

// --------------------------------------------------------------------------
// The value matcher: `+` and `,` separate a list, `-` gives an inclusive
// numeric range, `*` and `?` are wildcards, `\` escapes the next character.
// --------------------------------------------------------------------------

bool wildcard_match(const std::string& pattern, const std::string& text) {
    // `*` any run, `?` one character. Iterative, with a backtrack point, so a
    // pattern of several stars cannot go quadratic on a long name.
    std::size_t p = 0, t = 0, star = std::string::npos, mark = 0;
    while (t < text.size()) {
        if (p < pattern.size() &&
            (pattern[p] == '?' || pattern[p] == text[t])) {
            ++p; ++t;
        } else if (p < pattern.size() && pattern[p] == '*') {
            star = p++; mark = t;
        } else if (star != std::string::npos) {
            p = star + 1; t = ++mark;
        } else {
            return false;
        }
    }
    while (p < pattern.size() && pattern[p] == '*') ++p;
    return p == pattern.size();
}

//! Whether \p token is an integer, and its value.
bool as_int(const std::string& s, int& out) {
    if (s.empty()) return false;
    char* end = nullptr;
    const long v = std::strtol(s.c_str(), &end, 10);
    if (end == s.c_str() || *end != '\0') return false;
    out = static_cast<int>(v);
    return true;
}

//! One alternative of a value: a literal/wildcard, or a numeric range.
struct ValueTerm {
    std::string text;      //!< as written, minus quotes
    bool is_range;
    int lo, hi;
    bool has_wildcard;
    ValueTerm() : is_range(false), lo(0), hi(0), has_wildcard(false) {}
};

std::string unquote(const std::string& s) {
    if (s.size() >= 2 && (s[0] == '"' || s[0] == '\'') && s[s.size() - 1] == s[0]) {
        return s.substr(1, s.size() - 2);
    }
    return s;
}

//! Split a value into its alternatives.
std::vector<ValueTerm> parse_value(const std::string& raw) {
    std::vector<ValueTerm> out;
    const std::string text = unquote(raw);
    std::vector<std::string> parts;
    std::string cur;
    bool escape = false;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if (escape) { cur += c; escape = false; continue; }
        if (c == '\\') { escape = true; continue; }
        if (c == '+' || c == ',') { parts.push_back(cur); cur.clear(); continue; }
        cur += c;
    }
    parts.push_back(cur);

    for (const std::string& part : parts) {
        if (part.empty()) continue;
        ValueTerm t;
        t.text = part;
        t.has_wildcard = part.find('*') != std::string::npos ||
                         part.find('?') != std::string::npos;
        // A range is `lo-hi` with both ends integers. A leading `-` is a sign,
        // so the separator is searched from position 1.
        const std::size_t dash = part.find('-', 1);
        if (dash != std::string::npos && !t.has_wildcard) {
            int lo = 0, hi = 0;
            if (as_int(part.substr(0, dash), lo) &&
                as_int(part.substr(dash + 1), hi)) {
                t.is_range = true; t.lo = lo; t.hi = hi;
            }
        }
        out.push_back(t);
    }
    return out;
}

//! Whether one of \p terms matches \p text (and \p value, for a range).
bool value_matches(const std::vector<ValueTerm>& terms, const std::string& text,
                   bool numeric, int value, bool ignore_case) {
    const std::string t = ignore_case ? upper(text) : text;
    for (const ValueTerm& term : terms) {
        if (numeric && term.is_range) {
            if (value >= term.lo && value <= term.hi) return true;
            continue;
        }
        const std::string p = ignore_case ? upper(term.text) : term.text;
        if (term.has_wildcard) {
            if (wildcard_match(p, t)) return true;
        } else if (p == t) {
            return true;
        } else if (numeric) {
            // `resi 5` must match the atom whose number is 5 whatever the
            // text spelling is (`5`, ` 5`), so compare numerically too.
            int v = 0;
            if (as_int(term.text, v) && v == value) return true;
        }
    }
    return false;
}

}  // namespace

// --------------------------------------------------------------------------
// The tokenizer
// --------------------------------------------------------------------------

std::vector<std::string> selection_tokens(const std::string& expression) {
    std::vector<std::string> out;
    std::string cur;
    bool in_word = false, in_quote = false;
    char quote_char = '"';
    for (std::size_t i = 0; i < expression.size(); ++i) {
        const char c = expression[i];
        if (in_quote) {
            cur += c;
            if (c == quote_char) in_quote = false;
            continue;
        }
        switch (c) {
            case '"':
            case '\'':
                in_quote = true; quote_char = c; in_word = true; cur += c;
                break;
            case ' ':
            case '\t':
            case '\n':
                if (in_word) { out.push_back(cur); cur.clear(); in_word = false; }
                break;
            case ';':  // a word terminator that stays part of the word
                cur += c; out.push_back(cur); cur.clear(); in_word = false;
                break;
            case '!':
            case '&':
            case '|':
            case '(':
            case ')':
            case '%':
                if (in_word) { out.push_back(cur); cur.clear(); in_word = false; }
                out.push_back(std::string(1, c));
                break;
            default:
                cur += c; in_word = true;
                break;
        }
    }
    if (in_word && !cur.empty()) out.push_back(cur);
    return out;
}

// --------------------------------------------------------------------------
// Which spelling an expression is written in
// --------------------------------------------------------------------------

SelectionDialect selection_dialect(const std::string& expression) {
    const std::vector<std::string> tokens = selection_tokens(expression);
    bool plus = false, space = false;

    // Keywords only one vocabulary has. `resid`, `chain`, `name`, `index`,
    // `within` and the logical operators are in both and count for neither.
    static const std::set<std::string> plus_words = {
        "resi", "resn", "segi", "byres", "byresi", "byresidue", "br.",
        "hetatm", "het", "solvent", "sol.", "alt", "altloc", "idx.",
        "n.", "i.", "c.", "e.", "r.", "s.", "h.", "bb.", "sc.", "pol.",
        "pro.", "nuc.", "org.", "ino.", "w.", "be.", "a.", "id"};
    static const std::set<std::string> space_words = {
        "resname", "segname", "serial", "same", "water", "waters", "noh",
        "noself", "pbc", "to", "of"};

    bool after_property = false;
    for (std::size_t i = 0; i < tokens.size(); ++i) {
        const std::string w = lower(tokens[i]);
        if (plus_words.count(w)) plus = true;
        if (space_words.count(w)) space = true;
        // `+` inside a value is the one vocabulary's list separator; a second
        // bare word after a value is the other's.
        const Keyword* kw = keyword_table().count(w) ? &keyword_table().at(w)
                                                     : nullptr;
        if (kw != nullptr && kw->kind == KW_PROPERTY) {
            after_property = true;
            continue;
        }
        if (after_property) {
            if (w.find('+') != std::string::npos) plus = true;
            after_property = false;
            // a second value word, with no operator between: the space list
            if (i + 1 < tokens.size()) {
                const std::string next = lower(tokens[i + 1]);
                if (!keyword_table().count(next) && next != "(" && next != ")") {
                    space = true;
                }
            }
            continue;
        }
    }
    // `of` alone is not a marker -- both spell `within R of SEL` -- so it only
    // counts when nothing else does.
    if (plus && space) return SELECTION_DIALECT_MIXED;
    if (plus) return SELECTION_DIALECT_PLUS_LIST;
    if (space) return SELECTION_DIALECT_SPACE_LIST;
    return SELECTION_DIALECT_UNMARKED;
}

// --------------------------------------------------------------------------
// The parse tree
// --------------------------------------------------------------------------

struct SelectionExpression::Node {
    enum Type { PROPERTY, FLAG, NOT, AND, OR, ANDNOT, BYRES, BYCHAIN,
                DISTANCE } type;
    int what;                       //!< a Property, a Flag or a Distance
    std::vector<ValueTerm> values;  //!< for PROPERTY
    double radius;                  //!< for DISTANCE
    bool noself;                    //!< for DISTANCE
    std::shared_ptr<Node> a, b;
    Node() : type(FLAG), what(F_NONE), radius(0.0), noself(false) {}
};

namespace {

// Not `Node`: this is an unnamed namespace, whose members are visible at
// IMP::bff scope in the unity build, and IMP::bff::GraphNode is the reactive
// chinet node (Port.h/Node.h) -- a bare `Node` here made every later file
// in bff_all.cpp ambiguous.
typedef SelectionExpression::Node ParseNode;
typedef std::shared_ptr<ParseNode> NodePtr;

//! A recursive-descent parser over the token list.
class Parser {
public:
    Parser(const std::vector<std::string>& tokens, const std::string& source)
        : t_(tokens), at_(0), source_(source) {}

    NodePtr parse() {
        if (t_.empty()) return NodePtr();
        NodePtr n = parse_expression(0);
        if (at_ < t_.size()) {
            fail("unexpected `" + t_[at_] + "`");
        }
        return n;
    }

private:
    const std::vector<std::string>& t_;
    std::size_t at_;
    std::string source_;

    void fail(const std::string& why) const {
        IMP_THROW("selection `" << source_ << "`: " << why, ValueException);
    }

    bool done() const { return at_ >= t_.size(); }
    const std::string& peek() const { return t_[at_]; }

    static const Keyword* lookup(const std::string& word) {
        const auto& table = keyword_table();
        auto it = table.find(lower(word));
        if (it != table.end()) return &it->second;
        // `id` is also spelled `ID`; everything else is lower case.
        it = table.find(word);
        return it == table.end() ? nullptr : &it->second;
    }

    //! Everything at precedence \p min_precedence or tighter.
    NodePtr parse_expression(int min_precedence) {
        NodePtr left = parse_unary();
        while (!done()) {
            const Keyword* kw = lookup(peek());
            if (kw != nullptr && kw->kind == KW_DISTANCE &&
                0x30 >= min_precedence) {
                left = parse_distance(left);
                continue;
            }
            if (kw == nullptr || kw->kind != KW_BINARY) break;
            if (kw->precedence < min_precedence) break;
            const int op = kw->what;
            const int prec = kw->precedence;
            ++at_;
            if (done()) fail("`" + t_[at_ - 1] + "` needs a right-hand side");
            NodePtr right = parse_expression(prec + 1);
            NodePtr n(new ParseNode());
            n->type = op == OP_AND ? ParseNode::AND
                                   : (op == OP_OR ? ParseNode::OR : ParseNode::ANDNOT);
            n->a = left; n->b = right;
            left = n;
        }
        return left;
    }

    NodePtr parse_unary() {
        if (done()) fail("expression ends where a selection was expected");
        const Keyword* kw = lookup(peek());
        if (kw != nullptr && kw->kind == KW_SAME) {
            // `same residue as SEL`, `same chain as SEL`
            ++at_;
            if (done()) fail("`same` needs a property and `as`");
            const std::string prop = lower(t_[at_++]);
            if (done() || lower(t_[at_]) != "as") {
                fail("`same " + prop + "` needs `as` before its selection");
            }
            ++at_;
            NodePtr operand = parse_expression(0x20);
            NodePtr n(new ParseNode());
            if (prop == "residue" || prop == "resid" || prop == "resi") {
                n->type = ParseNode::BYRES;
            } else if (prop == "chain") {
                n->type = ParseNode::BYCHAIN;
            } else {
                fail("`same " + prop + " as` is not supported -- `residue` and"
                     " `chain` are");
            }
            n->a = operand;
            return n;
        }
        if (kw != nullptr && kw->kind == KW_DISTANCE) {
            return parse_distance(NodePtr());
        }
        if (kw != nullptr && kw->kind == KW_UNARY) {
            const int op = kw->what;
            const int prec = kw->precedence;
            ++at_;
            // `not` binds tightly, `byres` binds loosest -- both parse their
            // operand at their own precedence, which is what makes
            // `byres chain A and name CA` extend the whole conjunction.
            NodePtr operand = parse_expression(prec);
            NodePtr n(new ParseNode());
            n->type = op == OP_NOT ? ParseNode::NOT : ParseNode::BYRES;
            n->a = operand;
            return n;
        }
        return parse_primary();
    }

    //! `within R [pbc] [noself] of SEL`, `beyond R of SEL`, `SEL around R`.
    /*! \p left is the left-hand selection of the infix form, or null. */
    NodePtr parse_distance(NodePtr left) {
        const Keyword* kw = lookup(t_[at_]);
        const std::string word = t_[at_];
        const int what = kw->what;
        ++at_;
        if (done()) fail("`" + word + "` needs a distance");
        double radius = 0.0;
        {
            char* end = nullptr;
            radius = std::strtod(t_[at_].c_str(), &end);
            if (end == t_[at_].c_str() || *end != '\0') {
                fail("`" + word + "` needs a number, not `" + t_[at_] + "`");
            }
            ++at_;
        }
        bool noself = what == D_AROUND;
        while (!done()) {
            const std::string q = lower(t_[at_]);
            if (q == "noself") { noself = true; ++at_; continue; }
            if (q == "pbc") {
                fail("`pbc` asks for periodic wrapping, and there is no cell"
                     " here to wrap in");
            }
            break;
        }
        NodePtr n(new ParseNode());
        n->type = ParseNode::DISTANCE;
        n->what = what;
        n->radius = radius;
        n->noself = noself;

        if (what == D_AROUND) {
            // `SEL around R` -- the operand is what came before.
            if (!left) fail("`around` needs a selection before it");
            n->a = left;
            return n;
        }
        if (done() || lower(t_[at_]) != "of") {
            fail("`" + word + " " + std::to_string(radius) + "` needs `of`"
                 " and a selection");
        }
        ++at_;
        n->a = parse_expression(0x60);
        if (!left) return n;
        NodePtr both(new ParseNode());
        both->type = ParseNode::AND;
        both->a = left;
        both->b = n;
        return both;
    }

    NodePtr parse_primary() {
        if (done()) fail("expression ends where a selection was expected");
        const std::string word = t_[at_];

        if (word == "(") {
            ++at_;
            NodePtr n = parse_expression(0);
            if (done() || t_[at_] != ")") fail("unbalanced `(`");
            ++at_;
            return n;
        }
        if (word == ")") fail("unbalanced `)`");

        const Keyword* kw = lookup(word);
        if (kw == nullptr) {
            IMP_THROW("selection `" << source_ << "`: `" << word
                      << "` is not a selection keyword. This layer selects on"
                         " a structure's own fields (name, resn, resi, chain,"
                         " segi, elem, alt, index, id) and knows no named"
                         " objects or selections.", ValueException);
        }
        if (kw->kind == KW_UNSUPPORTED) {
            IMP_THROW("selection `" << source_ << "`: `" << word
                      << "` is valid notation this layer cannot answer -- it"
                         " needs coordinates, bonds, or a viewer's object"
                         " model, none of which a mask carries.",
                      ValueException);
        }
        ++at_;

        if (kw->kind == KW_FLAG) {
            NodePtr n(new ParseNode());
            n->type = ParseNode::FLAG;
            n->what = kw->what;
            return n;
        }

        // A property keyword: one value, plus any bare words that follow it.
        if (done()) fail("`" + word + "` needs a value");
        NodePtr n(new ParseNode());
        n->type = ParseNode::PROPERTY;
        n->what = kw->what;
        n->values = parse_value(t_[at_]);
        ++at_;
        while (!done()) {
            const std::string& next = t_[at_];
            if (next == "(" || next == ")") break;
            // `resid 50 to 60` -- the other spelling of `resid 50-60`.
            if (lower(next) == "to" && at_ + 1 < t_.size() &&
                !n->values.empty()) {
                int lo = 0, hi = 0;
                if (as_int(n->values.back().text, lo) &&
                    as_int(t_[at_ + 1], hi)) {
                    ValueTerm range;
                    range.is_range = true; range.lo = lo; range.hi = hi;
                    range.text = n->values.back().text + "-" + t_[at_ + 1];
                    n->values.back() = range;
                    at_ += 2;
                    continue;
                }
                fail("`to` needs a number on each side");
            }
            if (lookup(next) != nullptr) break;
            // A bare word extends the value list: `resn HOH SOL WAT`.
            const std::vector<ValueTerm> more = parse_value(next);
            n->values.insert(n->values.end(), more.begin(), more.end());
            ++at_;
        }
        if (n->values.empty()) fail("`" + word + "` needs a value");
        return n;
    }
};

//! Whether the node tree contains a `byres`.
//! `index` from 0 becomes `index` from 1, once, at parse time.
void shift_index_base(const NodePtr& n) {
    if (!n) return;
    if (n->type == ParseNode::PROPERTY && n->what == P_INDEX) {
        for (ValueTerm& v : n->values) {
            if (v.is_range) { v.lo += 1; v.hi += 1; }
            int one = 0;
            if (as_int(v.text, one)) v.text = std::to_string(one + 1);
        }
    }
    shift_index_base(n->a);
    shift_index_base(n->b);
}

bool needs_whole_structure(const NodePtr& n) {
    if (!n) return false;
    if (n->type == ParseNode::BYRES || n->type == ParseNode::BYCHAIN ||
        n->type == ParseNode::DISTANCE) {
        return true;
    }
    return needs_whole_structure(n->a) || needs_whole_structure(n->b);
}

bool flag_matches(int flag, const SelectionAtom& a) {
    switch (flag) {
        case F_ALL: return true;
        case F_NONE: return false;
        case F_HYDRO: {
            const std::string e = upper(a.elem);
            return e == "H" || e == "D";
        }
        case F_NOH: {
            const std::string e = upper(a.elem);
            return !(e == "H" || e == "D");
        }
        case F_HETATM: return a.hetatm;
        case F_PROTEIN: return is_protein(a);
        case F_NUCLEIC: return is_nucleic(a);
        case F_POLYMER: return is_protein(a) || is_nucleic(a);
        case F_SOLVENT: return is_solvent(a);
        case F_ORGANIC:
            return !is_protein(a) && !is_nucleic(a) && !is_solvent(a) &&
                   upper(a.elem) == "C";
        case F_INORGANIC:
            return !is_protein(a) && !is_nucleic(a) && !is_solvent(a) &&
                   upper(a.elem) != "C";
        case F_BACKBONE:
            return (is_protein(a) && protein_backbone().count(upper(a.name))) ||
                   (is_nucleic(a) && nucleic_backbone().count(upper(a.name)));
        case F_SIDECHAIN:
            return (is_protein(a) || is_nucleic(a)) && !flag_matches(F_BACKBONE, a);
        case F_GUIDE: return is_protein(a) && upper(a.name) == "CA";
        default: return false;
    }
}

bool property_matches(const ParseNode& n, const SelectionAtom& a) {
    switch (n.what) {
        case P_NAME: return value_matches(n.values, a.name, false, 0, true);
        case P_RESN: return value_matches(n.values, a.resn, false, 0, true);
        case P_ELEM: return value_matches(n.values, a.elem, false, 0, true);
        case P_ALT: return value_matches(n.values, a.alt, false, 0, true);
        // chain and segi are case-sensitive: `a` and `A` are different
        // chains.
        case P_CHAIN: return value_matches(n.values, a.chain, false, 0, false);
        case P_SEGI: return value_matches(n.values, a.segi, false, 0, false);
        case P_RESI: {
            const std::string text = a.resi_text.empty()
                    ? std::to_string(a.resi) : a.resi_text;
            return value_matches(n.values, text, true, a.resi, true);
        }
        case P_INDEX:
            return value_matches(n.values, std::to_string(a.index), true,
                                 a.index, true);
        case P_ID:
            return value_matches(n.values, std::to_string(a.id), true, a.id,
                                 true);
        default: return false;
    }
}

void evaluate_node(const NodePtr& n, const std::vector<SelectionAtom>& atoms,
                   std::vector<char>& out) {
    out.assign(atoms.size(), 0);
    if (!n) return;
    switch (n->type) {
        case ParseNode::PROPERTY:
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                out[i] = property_matches(*n, atoms[i]) ? 1 : 0;
            }
            break;
        case ParseNode::FLAG:
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                out[i] = flag_matches(n->what, atoms[i]) ? 1 : 0;
            }
            break;
        case ParseNode::NOT: {
            std::vector<char> a;
            evaluate_node(n->a, atoms, a);
            for (std::size_t i = 0; i < atoms.size(); ++i) out[i] = a[i] ? 0 : 1;
            break;
        }
        case ParseNode::AND:
        case ParseNode::OR:
        case ParseNode::ANDNOT: {
            std::vector<char> a, b;
            evaluate_node(n->a, atoms, a);
            evaluate_node(n->b, atoms, b);
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                out[i] = n->type == ParseNode::AND ? (a[i] && b[i])
                       : n->type == ParseNode::OR  ? (a[i] || b[i])
                                              : (a[i] && !b[i]);
            }
            break;
        }
        case ParseNode::DISTANCE: {
            std::vector<char> a;
            evaluate_node(n->a, atoms, a);
            const double r2 = n->radius * n->radius;
            // Brute force over the reference set. A structure's atom count is
            // in the thousands and a reference set is usually far smaller; a
            // grid would pay for itself only well past that.
            std::vector<std::size_t> reference;
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                if (a[i]) reference.push_back(i);
            }
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                bool is_near = false;   // `near` is a word MSVC keeps (see RRT.cpp)
                for (std::size_t k = 0; k < reference.size() && !is_near; ++k) {
                    const std::size_t j = reference[k];
                    if (n->noself && j == i) continue;
                    const double dx = atoms[i].x - atoms[j].x;
                    const double dy = atoms[i].y - atoms[j].y;
                    const double dz = atoms[i].z - atoms[j].z;
                    is_near = (dx * dx + dy * dy + dz * dz) <= r2;
                }
                if (n->what == D_BEYOND) {
                    out[i] = is_near ? 0 : 1;
                } else {
                    out[i] = is_near ? 1 : 0;
                }
            }
            break;
        }
        case ParseNode::BYCHAIN: {
            std::vector<char> a;
            evaluate_node(n->a, atoms, a);
            std::set<std::string> hit;
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                if (a[i]) hit.insert(atoms[i].chain + "\x1f" + atoms[i].segi);
            }
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                out[i] = hit.count(atoms[i].chain + "\x1f" + atoms[i].segi)
                                 ? 1 : 0;
            }
            break;
        }
        case ParseNode::BYRES: {
            std::vector<char> a;
            evaluate_node(n->a, atoms, a);
            std::set<std::string> hit;
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                if (!a[i]) continue;
                hit.insert(atoms[i].chain + "\x1f" + atoms[i].segi + "\x1f" +
                           (atoms[i].resi_text.empty()
                                    ? std::to_string(atoms[i].resi)
                                    : atoms[i].resi_text));
            }
            for (std::size_t i = 0; i < atoms.size(); ++i) {
                const std::string key =
                        atoms[i].chain + "\x1f" + atoms[i].segi + "\x1f" +
                        (atoms[i].resi_text.empty()
                                 ? std::to_string(atoms[i].resi)
                                 : atoms[i].resi_text);
                out[i] = hit.count(key) ? 1 : 0;
            }
            break;
        }
    }
}

}  // namespace

// --------------------------------------------------------------------------
// SelectionExpression
// --------------------------------------------------------------------------

SelectionExpression::SelectionExpression()
    : dialect_(SELECTION_DIALECT_UNMARKED) {}

SelectionExpression::SelectionExpression(const std::string& expression,
                                         SelectionDialect dialect)
    : expression_(expression), dialect_(dialect) {
    if (dialect_ == SELECTION_DIALECT_AUTO) {
        dialect_ = selection_dialect(expression);
    }
    const std::vector<std::string> tokens = selection_tokens(expression);
    Parser parser(tokens, expression);
    root_ = parser.parse();
    if (root_ && dialect_ == SELECTION_DIALECT_SPACE_LIST) {
        // `index` numbers atoms from 0 in this spelling and from 1 in the
        // other, so the same expression names different atoms. Shift the
        // values once, here, rather than at every comparison.
        shift_index_base(root_);
    }
}

SelectionExpression::SelectionExpression(const SelectionExpression& other)
    : expression_(other.expression_), dialect_(other.dialect_),
      root_(other.root_) {}

SelectionExpression& SelectionExpression::operator=(
        const SelectionExpression& other) {
    expression_ = other.expression_;
    dialect_ = other.dialect_;
    root_ = other.root_;
    return *this;
}

SelectionExpression::~SelectionExpression() {}

bool SelectionExpression::get_is_empty() const { return !root_; }

namespace {
bool uses_distance(const std::shared_ptr<SelectionExpression::Node>& n) {
    if (!n) return false;
    if (n->type == SelectionExpression::Node::DISTANCE) return true;
    return uses_distance(n->a) || uses_distance(n->b);
}
}  // namespace

bool SelectionExpression::get_needs_coordinates() const {
    return uses_distance(root_);
}

std::vector<int> SelectionExpression::evaluate(
        const std::vector<SelectionAtom>& atoms) const {
    std::vector<char> mask;
    evaluate_node(root_, atoms, mask);
    return std::vector<int>(mask.begin(), mask.end());
}

bool SelectionExpression::matches(const SelectionAtom& atom) const {
    if (!root_) return false;
    if (needs_whole_structure(root_)) {
        IMP_THROW("selection `" << expression_ << "` asks a question about a"
                  " residue, a chain or a distance, which no single atom can"
                  " answer -- use evaluate() over the whole structure.",
                  ValueException);
    }
    std::vector<SelectionAtom> one(1, atom);
    std::vector<char> mask;
    evaluate_node(root_, one, mask);
    return mask[0] != 0;
}

// --------------------------------------------------------------------------
// Structures in
// --------------------------------------------------------------------------

std::vector<std::string> protein_backbone_atom_names() {
    return protein_backbone_order();
}

// --------------------------------------------------------------------------
// Compiling to IMP's own selection algebra
// --------------------------------------------------------------------------



IMPBFF_END_NAMESPACE
