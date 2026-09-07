---
type: validation
title: "Retire `bff::Table`: it is a strict subset of `tttrlib::DataStore`, in the wrong repository, with one caller and that caller is its own test"
description: Feature-by-feature comparison of bff::Table against tttrlib::DataStore, and a search for callers across imp.bff, chisurf, imp-tricks and tttrlib. Table has zero consumers outside test/table/test_table.py; ndxplorer -- the migration target the design note wrote it for -- already calls DataStore.select_expression. Every Table capability except add_derived_column exists in DataStore, usually in a stronger form, and the layering rule puts tabular measurement data in tttrlib on both halves of the test. Decision 2026-08-31 - retire Table and delete Expression::compute_pointers with it; add one derived-column primitive to DataStore to close the only real gap.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, tttrlib, expression, table, datastore, layering, duplication, decision]
timestamp: '2026-08-31T00:00:00Z'
---

# `bff::Table` should be retired

**Decision: retire.** Not narrow. Delete `include/Table.h`,
`src/standalone/Table.cpp`, `test/table/test_table.py`, the nine lines of
SWIG that wrap them, and — with them — `Expression::compute_pointers`, whose
only caller anywhere is `Table.cpp`. `Expression`, `ChiSquared` and `Sampler`
are untouched: this ends Table's **storage** role, not the expression engine.

Three things force it, in descending order of weight.

1. **Table has no consumers.** A grep of `imp.bff`, `chisurf` (including the
   `modules/ndxplorer` submodule), `imp-tricks` and `tttrlib` for
   `bff.Table` / `bff::Table` returns exactly one non-build, non-note hit:
   `/Users/tpeulen/dev/imp.bff/test/table/test_table.py`. The class is
   exercised solely by the test written for it.
2. **The migration it was built to enable already happened, to the other
   class.** `tttrlib/okf/design-expression-selection.md` names
   `ndxplorer.core.DataSource.query_mask` and `_bff_query_table` as what
   `Table` would replace. `_bff_query_table` no longer exists, and
   `query_mask` at
   `/Users/tpeulen/dev/chisurf/modules/ndxplorer/ndxplorer/core/data_source.py:891`
   now reads `store.select_expression(query)` on a `tttrlib.DataStore`.
   Table's intended consumer went to DataStore without it.
3. **The layering rule points one way on both halves.** See
   [The layering rule, applied](#the-layering-rule-applied).

## The overlap, method by method

Read from `/Users/tpeulen/dev/imp.bff/include/Table.h` (135 lines),
`/Users/tpeulen/dev/imp.bff/src/standalone/Table.cpp` (330 lines),
`/Users/tpeulen/dev/tttrlib/modules/core/include/DataStore.h` (2435 lines) and
`/Users/tpeulen/dev/tttrlib/modules/core/src/DataStore.cpp` (690 lines).

### What Table does, and where DataStore already does it

| `bff::Table` | `tttrlib::DataStore` equivalent | verdict |
|---|---|---|
| `add_column(name, vector<double>)` | `add_column(name, ColumnType)` + `Column::set_f64`; Python `store.add(name, values)` | superset — 10 numeric types, `Bool` bitset, dictionary-encoded `String` |
| `add_column_array(name, double*, int)` | same, through the numpy typemap | equal |
| `add_column_float32_view(name, float*, int)` | none — DataStore always **owns** | the one structural difference; see below |
| `is_float32_column(name)` | `Column::type() == ColumnType::Float32` | superset |
| `has_column`, `get_column_names`, `get_number_of_columns`, `get_number_of_rows` | `find`, `column_names`, `n_columns`, `n_rows` | equal |
| `get_column(name, double**, int*)` — copies to a malloc'd double array | `Column::get_f64_view` / `get_f32_view` / … — **zero-copy** typed view, with the root-store keep-alive wired in `DataStore.i` | DataStore is better *and* cheaper |
| `select(query, double**, int*)` — one **double per row** | `select_expression(expr, Combine)` writing a bit-packed `BitMask`; `expression_mask(expr)` for the mask alone; `selection()` in Python | DataStore: 1 bit/row vs 64, and it composes |
| `count(query)` | `count_expression(expr)` | equal |
| `select_column(name, query, double**, int*)` | `select_expression` then `take(rows)` / `take_into` / a numpy view indexed by `selection()` | equal-to-better (DataStore takes **every** column, typed) |
| `add_derived_column(name, query)` | **nothing** | the only genuine gap |
| `describe()` | `__repr__`, `memory_report()`, `nbytes()`, `DataStoreRegistry::list()` | superset |

### What DataStore has that Table has no answer to

Verified present in `DataStore.h`; none of it exists in `Table`:

| capability | DataStore |
|---|---|
| typed columns | `ColumnType::{Float64,Float32,Int64,Int32,Int16,Int8,UInt64,UInt32,UInt16,UInt8,Bool,String}` — Table is `double` and a borrowed `float` |
| missing values | `Column::mask()`, `valid(i)`, `has_missing()`, `na_ranges()` with a *reason* string, `mask_non_finite()`; `mark_selected` refuses a row invalid in any column it read |
| selection algebra | `Combine::{Replace,And,Or,AndNot}`, `invert_selection`, `select_all/none`, `select_all_recursive` |
| selectors beyond an expression | `select_range`, `select_equal`, `select_interval`, `select_finite`, `select_rectangle`, `select_ellipse`, `select_quadratic`, `select_polygon`, `select_mask_image`, `select_mask_rows`, `where` |
| bit-packed masks | `BitMask` with `and_with`/`or_with`/`andnot_with`/`invert`/`count`, word-at-a-time |
| joins | `append_rows(other, Join::{Outer,Inner})`, `append_columns(..., OnDuplicate::{Refuse,KeepFirst})`, columns matched by **name** |
| row reshaping | `take_into`, `compact_into`, `copy_into` |
| hierarchy | groups: `add_group`, `ensure_group`, `group(path)`, `group_paths`, `walk`, `glob`/`rglob`, `tree` |
| column metadata | `units()`, `metadata()` (JSON, msgpack on disk), `attribute`/`set_attribute` |
| aggregation | `histogram`, `profile`, `fill_with`/`fill_sample_with` |
| accounting | `DataStoreRegistry` — every live store, its bytes, its label |
| bindings | Python, R, Java, JavaScript from one `.i`; Table is Python-only |
| integer correctness | `check_exactly_representable` refuses an `Int64`/`UInt64` beyond 2^53 rather than answer an `==` wrongly (`DataStore.cpp:448`) |

### The float32 view, honestly

`add_column_float32_view` is the only thing Table can do that DataStore
cannot: **borrow** a `float*` that lives elsewhere. That capability exists
only because Table is a query engine bolted onto data someone else owns —
its own header says so ("For data already held in C++ — a tttrlib store's
float32 columns"). Once the store *is* the owner, borrowing has nothing left
to borrow from; DataStore instead hands **views out** (`get_f32_view`), which
is the same saving in the direction that actually occurs. It is not a reason
to keep the class.

And the borrowed path is narrower than it reads:

- `bind()` (`Table.cpp:171`) **throws** on any query mixing a float32 view
  with an owned float64 column: *"mixing it with float64 columns in one query
  is not supported"*. DataStore widens the referenced columns into a scratch
  buffer held with the compiled program and answers the query
  (`DataStore.cpp:592-620`).
- `get_column` (`Table.cpp:151`) searches `columns_` only, so a float32 view
  column **cannot be read back at all**.
- `select_column` (`Table.cpp:288`) and `add_derived_column` (`Table.cpp:268`)
  likewise go through `columns_`/`bind`, so both are unusable on a float32
  view — `select_column` throws `no column`, `add_derived_column` throws the
  mixing error. Only `select` and `count` handle the float32 path.

Two smaller defects found while reading, recorded because they are the
signature of a class nobody uses: `count()` re-tests `all_float32` at
`Table.cpp:245` **after** `bind()` has already thrown for exactly that case —
ten lines that cannot execute; and `Float32ViewTests` in
`test/table/test_table.py` is defined *after* the file's
`if __name__ == "__main__": unittest.main()`, so it runs under pytest and
never under a direct `python test_table.py`.

### The one place Table is ahead: which engine evaluates

This is the only substantive point in Table's favour and it is worth stating
precisely, because it is easy to lose.

| path | parser | evaluator |
|---|---|---|
| `Table::select` / `count` / `add_derived_column`, float64 columns | bff `Expression` (ExprTk parse) | **bff's block-vectorised RPN engine**, via `compute_pointers` → `evaluate_vector_program` |
| `Table::select_float32` (`Table.cpp:105`) | ExprTk `parser<float>` | **ExprTk**, bound to the borrowed buffers |
| `DataStore::expression_mask` (`DataStore.cpp:537`), all types | ExprTk `parser<double>`/`parser<float>` | **ExprTk**, into a `std::vector`, then `mark_selected` packs bits |

So Table's float64 path is the only in-tree consumer of the fast engine over
table columns, and DataStore is still on ExprTk — which
`okf/handover-expression-engine.md` records as the losing evaluator.

That does not save Table, for two reasons. First, **it is already slower in
the measurement that exists**: at 200k rows and three float32 columns the
design note records Table at 0.41–0.43 ms and the handover records
`DataStore::select_expression` at **0.32 ms** on the same shape — DataStore's
ExprTk-plus-`BitMask` beats Table's ExprTk-plus-doubles, because writing one
bit per row instead of eight bytes is worth more than the evaluator
difference at that size. (Both figures are quoted from those notes; I did not
re-measure — the build lock is held elsewhere. Both were taken on a loaded
machine, so treat them as the same order rather than a precise 1.3x.)
Second, Table cannot use the fast engine on the float32 columns anyway — it
falls back to ExprTk there, which is the case a photon-derived burst table
actually is.

The real conclusion from that row is a **separate, larger ticket**: the
block-vector engine belongs where the columns are. `imp.bff` has no C++
dependency on `tttrlib` (verified: `dependencies.py` lists only IMP modules,
no `tttrlib` in `CMakeLists.txt`, and every `tttrlib` string in `include/` and
`src/` is a comment or a spec citation), and the layering is
tttrlib → imp.bff, so `DataStore` **cannot** call `bff::Expression`. Making
DataStore fast means moving the evaluator down into tttrlib, not keeping a
second table up in bff to reach it. Retiring Table does not foreclose that;
it removes the thing that made it look optional.

## The layering rule, applied

`AGENTS.md:22`: *photons/curves → tttrlib, coordinates → imp.bff, neither →
chisurf*, and *when input and consumer disagree, the consumer wins*.

- **Input.** Table's columns are burst-wise and photon-derived parameters —
  green/red/blue intensities, ratios, lifetimes. Its own header opens with
  "Burst and photon data already live in C++." Photon-derived → **tttrlib**.
  Not one column Table has ever held is a coordinate.
- **Consumer.** ndxplorer gating and chisurf burst selection. Already calling
  `DataStore` (`data_source.py:891`, `:925`; `core/tttrlib_selection.py`).
  → **tttrlib**.

Input and consumer **agree**, so the tie-break never fires and there is no
room for taste. Table is in the wrong repository on both halves of the test.
It is worth noticing that the same rule *keeps* `Expression` in bff, by the
tie-break rather than by the input test: a model curve is a curve
(→ tttrlib), but its consumers are `ChiSquared` and `Sampler`, which are bff
nodes, and the consumer wins. The rule discriminates cleanly between the two
classes; it is not being bent to reach the answer.

## Who calls `Expression::compute_pointers` afterwards: nobody. Delete it.

`compute_pointers` is `%ignore`'d in `pyext/swig.i-in:427` ("Raw column
pointers are Table's business, not Python's"), so no Python caller is
possible; the four in-tree call sites are `Table.cpp:224, 259, 282, 302`.
With Table gone it is unreachable code that Python cannot even see — the
worst kind to leave behind.

It is also nearly redundant. Reading `Expression.cpp:1341-1466`,
`compute_columns` and `compute_pointers` differ only in binding shape:
`compute_columns` takes one contiguous `n_vars × n_rows` block and mallocs
its output; `compute_pointers` takes one pointer per column and writes into
the caller's buffer. Both resolve names onto compiled slots, both dispatch to
`evaluate_vector_program` when a vector program exists, both fall back to
ExprTk bound in place otherwise. The genuinely distinct property —
**separate, non-contiguous column buffers with no staging copy** — has no
owner left in bff once Table goes, and by the section above it should not
acquire one: a C++ owner of tabular columns in bff is the thing being
retired. Thirty lines, trivially reinstated from `compute_columns` if a case
ever appears.

`compute_columns` (numpy 2-D in, array out) and `compute_mask` (numpy 2-D in,
byte mask out) stay: they are the wrapped surface, they are what the 55
expression tests exercise, and they cover every array-shaped caller.

## The one thing that dies with Table: derived columns

`add_derived_column(name, query)` — evaluate an expression over the columns
and keep the result as a new column — has **no DataStore equivalent**
(grepped `DataStore.h`, `DataStore.cpp`, `ext/python/DataStore.py` for
`derived`, `add_expression_column`, `eval_expression`: no hits). ndxplorer
currently does this in Python (`compute_values`, `data_source.py:980`, via
chisurf's `CompiledExpression`), so nothing regresses today — but the
capability is real and the machinery to provide it already exists.

Add to `DataStore`, in tttrlib:

```cpp
/// Evaluate `expr` over the columns and keep the result as a new column.
int add_expression_column(const std::string& name, const std::string& expr,
                          ColumnType type = ColumnType::Float64);
```

`expression_mask` already builds and caches an `ExpressionProgram` whose
`program.result` is the per-row double value before thresholding; this writes
that vector into a new column instead of into `mark_selected`. The refactor
is to split the existing `expression_mask` into *evaluate* and *threshold*
and call the first from both. Roughly 30 lines, no new dependency, and it
lands the capability at the layer that owns the data.

## Migration

Nothing downstream migrates — there is nothing downstream. The work is a
deletion in `imp.bff` plus one addition in `tttrlib`.

**imp.bff — delete**

| file | change |
|---|---|
| `include/Table.h` | delete (135 lines) |
| `src/standalone/Table.cpp` | delete (330 lines) |
| `src/Files.cmake` | drop `standalone/Table.cpp` from `cppfiles` |
| `src/CMakeLists.txt:135` | drop the `include/Table.h` entry |
| `pyext/swig.i-in:426-439` | drop the `%apply (float* IN_ARRAY1, int DIM1) {(float* in_floats, int n_floats)}`, `%apply (double* IN_ARRAY1, int DIM1) {(double* in_column, int n_column)}`, `%shared_ptr(IMP::bff::Table)`, `%include "IMP/bff/Table.h"` and the `%ignore …compute_pointers` line, with its comment |
| `test/table/test_table.py` | delete (194 lines, 21 tests) |
| `include/Expression.h:120` | remove the `compute_pointers` declaration and its "Not exposed to Python — `Table` is the caller" doc paragraph |
| `src/standalone/Expression.cpp:1406-1466` | remove the definition |
| `include/Expression.h:62` | `get_translated_expression()` keeps its doc but loses the sentence "Table compiles this in single precision…"; the method itself stays — `DataStore` needs the same Python→ExprTk rewrite and has its own `to_exprtk_syntax` |
| `okf/handover-expression-engine.md:26-27, 181-185` | replace "Probably should be retired" and the open question with a pointer to this note |

Also check whether the `%apply` typemaps removed above are used by another
`%apply` consumer in the same `.i` before deleting them — `%apply` is keyed
by parameter name and `in_column`/`in_floats` are Table-specific, so they
should come out cleanly, but the SWIG build is the check.

**tttrlib — add**

- `DataStore::add_expression_column` as above, with a test asserting
  agreement with numpy on the same expression.
- Optional, and the reason a follow-up ticket is worth opening: point
  `expression_mask` at a block-vectorised evaluator rather than ExprTk. That
  is the substance of the handover's "Is ExprTk worth keeping?" question and
  is where the engine's 1.3–2.2x mask advantage would finally reach the data.

## What breaks

| # | breakage | severity |
|---|---|---|
| 1 | 21 tests in `test/table/test_table.py` disappear from the suite (1569 → ~1548) | none — they test only the deleted class. Its pandas-parity queries are already covered for the engine by `test/expression/test_expression.py` and `test_expression_mask.py`, and for the store by tttrlib's own expression tests |
| 2 | `IMP.bff.Table` vanishes from the Python module | no in-tree or sibling-repo caller; an out-of-tree notebook using it would break, and the replacement is `tttrlib.DataStore` |
| 3 | `Expression::compute_pointers` vanishes | no Python caller possible (it was `%ignore`d); no C++ caller after (1) |
| 4 | The 200k-row `bff::Table` baseline in `tttrlib/okf/design-expression-selection.md` becomes unreproducible | it has already been beaten by `DataStore::select_expression` on the same shape; the note should record that and drop the table, not keep the class alive to serve it |
| 5 | `add_derived_column` has no home until `DataStore::add_expression_column` lands | ndxplorer computes derived columns in Python today, so nothing regresses; sequence the tttrlib addition first if you want zero window |
| 6 | The bff→tttrlib float32-view bridge (numpy view of a store column handed to `Table`) stops existing | it has no user; ndxplorer already gates inside the store |

Nothing in `doc/`, `examples/`, `benchmark/` or `bin/` references `Table`
(grepped). `benchmark/expression_curves.py` and `benchmark/expression_mask.py`
drive `Expression` directly, not `Table`, so both keep working.

## Verified vs inferred

**Verified by reading the source in this working tree** (all files read
whole, per the caution that git history is not current): `Table.h`,
`Table.cpp`, `Expression.h`, `Expression.cpp:1341-1466`, `swig.i-in:395-445`,
`src/Files.cmake`, `src/CMakeLists.txt`, `test/table/test_table.py`,
`DataStore.h`, `DataStore.cpp:400-690`, `DataStore.i`, `DataStore.py`,
`ndxplorer/core/data_source.py:880-945`, `AGENTS.md:1-48`. The caller search
covered `imp.bff`, `chisurf` (with `modules/ndxplorer`), `imp-tricks`,
`tttrlib`, `quest`, `ucfret`, `fpsimp`, excluding `.git` and build
directories.

**Quoted, not re-measured:** the 0.41–0.43 ms Table figures
(`tttrlib/okf/design-expression-selection.md`), the 0.32 ms
`select_expression` figure and the mask-vs-double ratios
(`okf/handover-expression-engine.md`). No build was run and no benchmark was
taken for this note — the build lock was held elsewhere.

**Inferred:** that `add_expression_column` is ~30 lines (from the shape of
`expression_mask`, not from writing it); that the `%apply` typemap removals
are clean (parameter names are Table-specific, but only a SWIG build proves
it); that no out-of-tree notebook uses `IMP.bff.Table` (unknowable from here,
though the class is recent and undocumented).
