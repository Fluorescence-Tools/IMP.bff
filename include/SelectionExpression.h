/**
 *  \file IMP/bff/SelectionExpression.h
 *  \brief An atom selection written as text, compiled to an IMP selection.
 *
 * `chain A and resi 10-20 and not name CA+CB` is the notation a structural
 * biologist already types, and the notation the `strip_mask` field of an
 * fps.json is authored in. This is that language: a tokenizer, a precedence
 * parser, an evaluator over a flat atom table, and a compiler to
 * #IMP::atom::Selection.
 *
 * The vocabulary is what a structure itself carries -- names, residues,
 * chains, segments, elements, alternate locations, indices. Terms that would
 * need coordinates, bonds, or a viewer's object model (`within`, `around`,
 * `expand`, `neighbor`, `model`, `rep`, `color`, `state`) are refused by name
 * with the reason, rather than silently matching nothing: a selection that
 * quietly selects nothing is how a strip mask comes to be decoration.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_SELECTIONEXPRESSION_H
#define IMPBFF_SELECTIONEXPRESSION_H

#include <IMP/bff/bff_config.h>

#include <IMP/bff/IMPCompatibility.h>

#include <memory>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Which spelling of the language an expression is written in.
/*! The two vocabularies overlap almost everywhere, so an expression can be
    read without knowing which it is -- except for `index`, where they
    disagree about the first atom's number. Detecting the spelling is what
    settles that. */
enum SelectionDialect {
    //! Decide from the expression itself. The default.
    SELECTION_DIALECT_AUTO = 0,
    //! Nothing in the expression distinguishes the two (`chain A and name CA`).
    SELECTION_DIALECT_UNMARKED = 1,
    //! `resi`, `resn`, `segi`, `byres`, `n.`/`i.`/`c.`, `name CA+CB`.
    //! `index` counts from **1**.
    SELECTION_DIALECT_PLUS_LIST = 2,
    //! `resname`, `segname`, `serial`, `same residue as`, `50 to 60`,
    //! `name CA CB`, `water`, `noh`, `noself`. `index` counts from **0**.
    SELECTION_DIALECT_SPACE_LIST = 3,
    //! Markers of both, which is common in files written by more than one
    //! tool. Read as \c SELECTION_DIALECT_PLUS_LIST.
    SELECTION_DIALECT_MIXED = 4
};

//! Which spelling \p expression is written in, from its own markers.
/*! Neither vocabulary is preferred: what is counted is the keywords and
    separators only one of them has. `resid` counts for neither -- both take
    it. */
IMPBFFEXPORT SelectionDialect selection_dialect(const std::string& expression);

//! One atom, as a selection expression sees it.
/*! The fields a selection expression can name. `resi_text` carries the insertion
    code when there is one (`52A`), which is what `resi 52A` matches; `resi`
    is the bare number, which is what a range `resi 50-60` compares. */
struct IMPBFFEXPORT SelectionAtom {
    std::string name;       //!< atom name -- `name`, `n.`
    std::string resn;       //!< residue name -- `resn`, `resname`, `r.`
    std::string chain;      //!< chain id -- `chain`, `c.`
    std::string segi;       //!< segment id -- `segi`, `segid`, `s.`
    std::string elem;       //!< element symbol -- `elem`, `element`, `symbol`
    std::string alt;        //!< alternate location -- `alt`, `altloc`
    std::string resi_text;  //!< residue number as written, insertion code and all
    int resi;               //!< residue number -- `resi`, `resid`, `i.`
    int index;              //!< 1-based position in the table -- `index`, `idx.`
    int id;                 //!< the file's own serial -- `id`, `serial`
    bool hetatm;            //!< from a HETATM record -- `hetatm`, `het`
    double x, y, z;         //!< coordinates, for `within` and `beyond`

    SelectionAtom()
        : resi(0), index(0), id(0), hetatm(false), x(0.0), y(0.0), z(0.0) {}

    IMP_SHOWABLE_INLINE(SelectionAtom,
                        out << "SelectionAtom(" << chain << "/" << resn << resi
                            << "/" << name << ")");
};
IMP_VALUES(SelectionAtom, SelectionAtoms);

//! A parsed selection expression.
/*! Parsing happens once, in the constructor, and raises on anything the
    language does not accept or this layer cannot answer. Evaluating is a walk
    over the atoms.

    **Two dialects, one parser.** The same expression language is spelled two
    ways in this field, and both are accepted, including mixed:

    | | one spelling | the other |
    |---|---|---|
    | residue number | `resi 5`, `i. 5` | `resid 5` |
    | residue name | `resn ALA`, `r. ALA` | `resname ALA` |
    | a list | `name CA+CB` | `name CA CB` |
    | a range | `resi 50-60` | `resid 50 to 60` |
    | segment | `segi`, `segid` | `segname` |
    | file serial | `id` | `serial` |
    | whole residues | `byres SEL` | `same residue as SEL` |
    | water | `solvent` | `water` |
    | near something | `SEL within 5 of SEL2` | `within 5 of SEL2` |

    **Operators**, tightest first: `not`/`!`, `and`/`&`, `-` (and-not),
    `or`/`|`/`+`, then `byres`/`br.`/`same residue as` -- which binds
    *loosest*, so `byres chain A and name CA` extends the whole conjunction to
    complete residues. Parentheses group.

    **Selectors** taking a value: `name`/`n.`, `resn`/`resname`/`r.`,
    `resi`/`resid`/`residue`/`i.`, `chain`/`c.`,
    `segi`/`segid`/`segname`/`segment`/`s.`, `elem`/`element`/`symbol`/`e.`,
    `alt`/`altloc`, `index`/`idx.`, `id`/`serial`.

    **Selectors taking none**: `all`/`*`, `none`, `hydro`/`hydrogens`/`h.`,
    `noh`, `hetatm`/`het`, `polymer`/`pol.`, `protein`/`pro.`,
    `nucleic`/`nuc.`, `solvent`/`sol.`/`water`, `organic`/`org.`,
    `inorganic`/`ino.`, `backbone`/`bb.`, `sidechain`/`sc.`, `guide`.

    **Distance**: `within R of SEL` selects what lies within \p R Angstrom of
    `SEL`; `beyond R of SEL` the complement. `noself` excludes `SEL`'s own
    atoms, and `SEL1 within R of SEL2` is `SEL1 and (within R of SEL2)`.
    `SEL around R` is `within R of SEL` with `noself`. A `pbc` qualifier is
    refused rather than ignored: there is no cell here to wrap in.

    **Values**: `+` and `,` separate a list (`name CA+CB`), `-` and `to` give
    an inclusive numeric range (`resi 50-60`, `resid 50 to 60`), `*` and `?`
    are wildcards (`name C*`), and a value may be quoted when it contains any
    of those (`resn 'Cl-'`). Matching is case-insensitive except for `chain`
    and `segi`, where `a` and `A` are different chains.

    \note `index` counts from **1**, in the order the structure lists its
    atoms, and `id`/`serial` is the number the file itself carries. A viewer
    that counts `index` from 0 will disagree by one; where that matters, say
    `serial`.

    \note A bare word following a selector's value extends that value list, so
    `resn HOH SOL WAT` reads as `resn HOH+SOL+WAT`. A viewer would read the
    trailing words as names of objects or named selections; this module has no
    object table for them to name, and every `strip_mask` shipped with an
    fps.json is written in the space-separated form. */
class IMPBFFEXPORT SelectionExpression {
public:
    SelectionExpression();
    //! Parse \p expression. Empty selects nothing.
    /*! \throw ValueException on a syntax error, an unknown keyword, or a
               keyword this layer cannot answer -- the message names the token
               and, for the last case, says what it would need. */
    explicit SelectionExpression(const std::string& expression,
                                 SelectionDialect dialect = SELECTION_DIALECT_AUTO);
    SelectionExpression(const SelectionExpression& other);
    SelectionExpression& operator=(const SelectionExpression& other);
    ~SelectionExpression();

    //! 1 where the expression selects the atom, 0 elsewhere.
    /*! `byres` needs the whole table -- it extends a hit to every atom sharing
        its (chain, segi, resi, resn) -- which is why this takes all of them
        and #matches does not. */
    std::vector<int> evaluate(const std::vector<SelectionAtom>& atoms) const;

    //! Whether the expression selects \p atom, judged on that atom alone.
    /*! \throw ValueException when the expression contains `byres`, which no
               single atom can answer. */
    bool matches(const SelectionAtom& atom) const;

    //! The expression as given.
    std::string get_expression() const { return expression_; }

    //! The spelling this was read as -- detected, or as the caller declared.
    SelectionDialect get_dialect() const { return dialect_; }

    //! Whether anything was parsed at all.
    bool get_is_empty() const;

    //! Whether the expression measures a distance (`within`, `beyond`,
    //! `around`), and so cannot be answered by a table without coordinates.
    bool get_needs_coordinates() const;

    IMP_SHOWABLE_INLINE(SelectionExpression,
                        out << "SelectionExpression(" << expression_ << ")");

#ifndef SWIG
    struct Node;
    //! The parse tree, for the compiler in `SelectionExpression.cpp`.
    std::shared_ptr<Node> get_root() const { return root_; }
#endif

private:
#ifdef SWIG
    struct Node;
#endif
    std::string expression_;
    SelectionDialect dialect_;
    std::shared_ptr<Node> root_;
};
IMP_VALUES(SelectionExpression, SelectionExpressions);

//! The #IMP::atom::Selection an expression denotes.
/*! This is the language's real target: an expression compiles to IMP's own
    selection algebra, so what comes out composes with everything that takes a
    #IMP::atom::Selection -- restraints, rigid bodies, `get_selected_particles`
    -- rather than to a mask only this module understands.

    The mapping is term by term. `chain` becomes `set_chain_ids`, `resi` a
    `set_residue_indexes` (a range expanded), `resn` a `set_residue_types`,
    `name` a `set_atom_types`, `elem` a `set_element`; `and`, `or` and `-`
    become `set_intersection`, `set_union` and `set_difference`, and `not`
    subtracts from the whole hierarchy.

    Four kinds of term have no IMP predicate to compile to -- a wildcard
    (`name C*`), the fields IMP's hierarchy does not carry (`segi`, `alt`,
    `id`), `index`, and the `byres` operator. Those are evaluated over the
    atom table and enter the algebra as a #IMP::atom::Selection over exactly
    the particles that matched, which is still an IMP selection and composes
    the same way.

    \throw ValueException on a syntax error or an unsupported keyword. */

//! The indices of the atoms of \p hierarchy that \p expression selects.
/*! Indices into #IMP::atom::get_leaves order, which is the order
    #selection_atoms reports. */

//! The atoms a protein backbone is made of: `N`, `CA`, `C`, `O`, `OXT`.
/*! What `backbone` selects, and what a labelling site keeps when its side
    chain is stripped. One list, because two would drift. */
IMPBFFEXPORT std::vector<std::string> protein_backbone_atom_names();

//! The tokens of \p expression.
/*! Quoted runs stay whole, `!&|()%` are words of their own, and `;` ends a
    word. Exposed because a caller that reports a bad mask wants to point at
    the token. */
IMPBFFEXPORT std::vector<std::string> selection_tokens(
        const std::string& expression);

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_SELECTIONEXPRESSION_H
