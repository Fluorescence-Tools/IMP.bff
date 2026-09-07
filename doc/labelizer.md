Choosing a labelling site {#labelizer}
======================================

[TOC]

Before a FRET experiment there is a question that has nothing to do with
photons: **where does the dye go, and which pair do you measure?** Getting it
wrong is expensive in the way that matters — the protein is expressed, labelled
and measured before anyone finds out that the dye was buried, that it stuck to
a hydrophobic patch, or that the two sites never move relative to each other.

`IMP.bff` answers both questions with a native port of the label-site score of
Gebhardt *et al.*, *Nat. Commun.* **16**, 3305 (2025). Everything is prefixed
`ll_` / `Ll`, so the ported model is distinguishable from this package's own
physics without consulting a design note.

Why it is here at all
---------------------

The reference implementation is a Python package, and it does not run. Its
secondary-structure term shells out to DSSP and raises
`RuntimeError("Unknown platform")` on macOS outright; its **default**
solvent-exposure term shells out to MSMS, whose bundled binaries are 32-bit
ppc/i386 Mach-O; and it pins `labellib`, which is not used in this package. A
scoring model whose default term cannot be evaluated is not a model anyone can
check.

Meanwhile every piece of physics underneath it was already here — accessible
volumes, Förster radii, κ², PET quenchers, Shrake–Rupley surfaces, side-chain
packing — and what was missing was only the scoring layer on top. Nothing in
`include/` iterated residues and emitted a score column.

So the two external programs are replaced by native kernels
(#IMP::bff::ll_dssp, #IMP::bff::ll_residue_depth) and the difference is
**measured** rather than asserted. See `okf/validation/labelizer_ab.md`: on the
reference's own published output the secondary-structure and
cysteine-resemblance terms agree exactly on every residue of two proteins, and
the native surface matches MSMS with no bias.

What a score means
------------------

**A parameter score is a likelihood ratio, not a probability.** Each fitted
table holds

\f[
    \frac{P(\mathrm{labelable} \mid s)}{P(\mathrm{labelable})}
\f]

binned over an observable. A value above 1 says the observation makes a site
*more* likely to be labelable than the base rate. Values run from 0 to about
3.9, and because the terms are multiplied the combined score is **unbounded
above** — a cysteine on an exposed loop passes 2 easily. The `0.5` threshold
the pair layer filters on is a threshold on that unbounded quantity, not a
probability cut.

The published model combines four terms as a weighted geometric mean:

| term | tag | what it measures | why it matters |
|---|---|---|---|
| conservation | `cs` | evolutionary rate at the position | a conserved residue is doing something; mutating it to cysteine may cost function |
| solvent exposure | `se` | depth below the molecular surface | a buried position cannot be reached by a labelling reagent |
| cysteine resemblance | `cr` | how much the residue already resembles a cysteine | a conservative substitution perturbs least |
| secondary structure | `ss` | the local DSSP state | loops and turns tolerate a dye; helices and strands less so |

Two further terms — tryptophan proximity (`tp`) and charge environment (`ce`) —
are implemented and carry **weight zero in the published model**. They are
switched off there, not merely down-weighted. A third, methionine exclusion
(`me`), vetoes a site near an exposed methionine and is not part of the
published four.

Getting a number out
--------------------

    #include <IMP/bff/Labelizer.h>

    const std::vector<LlScore> scores = ll_score_structure(
            "protein.pdb", ll_model_paper(), LlOptions(), "grades.txt");

Each #IMP::bff::LlScore is one row of `_mmfdb_label_score`: a position, a
`score_type`, a value and a `status`. **A position that was not scored carries
no value at all** — `status` says why. That is not fastidiousness: the
reference writes `-1` for "excluded" and `0` for "no contribution" into the
same column as real scores, so its output cannot be read back, and the shipped
example's conservation column is degenerate as a direct consequence.

Conservation is **imported, never computed** (#IMP::bff::ll_read_consurf). An
alignment and a rate estimate are a different program; this reads its result.

Which pair to measure
---------------------

    LlFretOptions options;
    options.forster_radius = forster_radius(find_probe("Alexa488"),
                                            find_probe("Alexa647"));
    const std::vector<LlPairScore> pairs = ll_pair_scores(
            "protein.pdb", ll_combined_by_key(scores), options);

A pair is worth measuring when both sites are labelable **and** the dye–dye
distance sits where FRET responds to it — near \f$R_0\f$. With two
conformations (#IMP::bff::ll_pair_scores_two_states) the criterion changes: what
matters is that the distance *moves*.

Placing the dye is a cost ladder, and the rungs are not interchangeable:

| model | cost | what it is |
|---|---|---|
| #IMP::bff::PROBE_MODEL_CBETA | nothing | the Cβ itself; no dye at all |
| #IMP::bff::PROBE_MODEL_ALPHA_CONE | one neighbour search | the reference's analytic estimate of the dye's mean position |
| #IMP::bff::PROBE_MODEL_ACCESSIBLE_VOLUME | a grid search | a real accessible volume, this package's own |

**Do not leave `n_refine` at zero for real work.** The analytic cone
over-reaches by 2–5 Å, and because the pair score peaks sharply at
\f$R = R_0\f$ that is a large error in the quantity being ranked: measured on
one structure, *none* of the top five pairs by cone score were still in the top
five after rebuilding the clouds. The cheap screen is for deciding which pairs
are worth building volumes for, not for choosing between them.

The dyes
--------

\f$R_0\f$ is **derived, not supplied** — from the donor's quantum yield and the
overlap of its emission with the acceptor's absorption
(#IMP::bff::forster_radius). Naming the dyes makes the number traceable;
asserting `52.0` does not.

Names resolve loosely (#IMP::bff::resolve_probe_name), because the library keys
are vendor spellings — `AlexaFluor488`, `LumiprobeCy3`, `ATTO647N` — and nobody
types those. `Alexa488`, `Cy3` and `Atto647N` all work. Ambiguity raises rather
than guessing.

The library is a **probe** library (#IMP::bff::ProbeLibrary.h), not a dye
library: #IMP::bff::ProbeType distinguishes a dye from a fluorescent protein
from an EPR spin label, and the shipped rotamer stores are already partitioned
that way. It matters here because a spin label has no spectra by construction,
so #IMP::bff::forster_radius refuses one **by naming the probe type** rather
than reporting a missing spectrum — a category error and a missing file should
not read the same.

One unit trap, stated because it nearly shipped: `forster_radius` and the
container's #IMP::bff::probe_pto_forster_radius both return **Ångström**, and so
does everything that consumes an \f$R_0\f$ here. Spectra are in nanometres. An
\f$R_0\f$ near 5 is nanometres and is wrong by a factor of ten; near 50 it is
right.

The output is one file
----------------------

    ll_write_pto("scored.mmfdb.pto", "protein.pdb", scores, pairs, settings);

The reference writes six CSVs, four PDBs carrying a dimensionless score in the
B-factor column, a heat-map JSON and a zip. This writes one container
(#IMP::bff::ll_write_pto): the structure verbatim and recoverable against its
SHA-256, the scores as tables whose columns are named by MMFDB dictionary
items, the complete settings, and provenance edges tying them together. A
`.pto` is an EBML document, so a walker in any language can read it without
this package.

Reproducing the paper, defects and all
--------------------------------------

The default is #IMP::bff::LL_MODEL_PUBLISHED, which reproduces the reference
**including three defects**, because the paper's numbers were computed with
them: the joined label score of a pair uses `prod ** 0.5` where the geometric
mean is `prod ** (1/N)`; a term of weight *zero* that happens to read exactly
`0.0` still zeroes the whole score; and the configurable excluded amino acid is
written to one attribute and read from another, so it is always methionine.

#IMP::bff::LL_MODEL_CORRECTED fixes all three. Which was used is recorded in
the container, so a file always says which numbers it holds.

Where to look next
------------------

- `examples/labels/plot_labelizer_score.py` — the per-residue score, end to end.
- `examples/labels/plot_labelizer_fret_pair.py` — the two-conformation pair
  score on maltose-binding protein, with \f$R_0\f$ from named dyes.
- `imp_bff_labelizer` — the same from a command line.
- `okf/validation/labelizer_ab.md` — what agrees with the reference, and by how
  much.
- `okf/labelizer-correspondence.md` — which reference function became which
  symbol, and what was deliberately not ported.
