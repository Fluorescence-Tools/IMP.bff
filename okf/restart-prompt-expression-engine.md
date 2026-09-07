# Restart prompt

Paste the block below into a fresh session. Read the handover first:
`~/dev/imp.bff/okf/handover-expression-engine.md`.

---

I am building a vectorised expression engine in C++ that three things use:
chisurf's parse models (equations from YAML, evaluated per fit iteration over
curves of 100-4096 points), tttrlib's table gating (boolean queries over
1e5-1e7 rows, answering with a bit-packed row mask), and bff fits (a whole fit
as a C++ node graph so a sampler move costs no interpreter).

Read `~/dev/imp.bff/okf/handover-expression-engine.md` before doing anything.
It has the current state, the measurements, and -- importantly -- four
approaches that were implemented, measured and reverted. Do not repeat them.

**Environment first, before any benchmark.** A previous session wasted several
rounds on phantom slowness that turned out to be a build with no
`CMAKE_BUILD_TYPE`, and took numbers on a machine at load average 19-33.
So: confirm `CMAKE_BUILD_TYPE=Release` in every build directory you measure
in, check `uptime`, and if the machine is loaded, say so and stop rather than
report numbers.

**The work, in order:**

1. **Fix `Expression::compute_mask()`.** It is a placeholder that allocates a
   full double buffer, runs the double evaluator and converts -- more work,
   not less. The block loop already carries a typed boolean stack; make the
   per-block final write emit bytes directly. Then re-measure; the current
   published number for it is measuring the stub.

2. **Decide ExprTk's fate.** It is now only a parser and a fallback; its
   evaluator lost decisively to the block engine. It also refuses valid input
   (>100 nested parens, `not x>2`). Establish what the three consumers
   actually need from a parser, and if our own tokenizer covers it, delete
   ExprTk -- it is 1.6 MB and ~20s of compile time. If something genuinely
   needs it, write down what, in the handover.

3. **Bit-packed output.** tttrlib's `DataStore` already answers selections
   with a `BitMask`, one bit per row -- 64x denser than numpy's bool array.
   That is where this engine can beat numpy rather than trail it, and it is
   the reason to put the engine in tttrlib at all.

4. **Common subexpression elimination.** FCS-shaped equations
   (`0.3+1/1.7*(1+x/1.2)**(-1)/sqrt(1+1/2.1**2*x/1.2)`) lose from ~512 points
   up. They compute `x/1.2` twice and repeat `1+...` structures. CSE in the
   RPN compiler should roughly halve the work without touching accuracy.

5. **Only then port to tttrlib**, replacing ExprTk in
   `DataStore::select_expression`, and re-measure against pandas. See
   `~/dev/tttrlib/okf/design-expression-selection.md` for the design and the
   seam (`scan_column`/`scan_typed`/`BitMask`/`Combine`).

**Standing requirements.** The engine must never crash the interpreter -- a
malformed equation from a YAML file or a query box has to raise, not
segfault. There is a robustness suite
(`test/expression/test_expression_robustness.py`) and a fuzz harness; keep
both green and extend them. All 86 shipped equations must keep matching numpy
to 1e-10. Do not report a speedup without saying what the machine load was.

**Also open, unrelated to speed:**
- `bff::Table` duplicates `tttrlib::DataStore` (which has typed columns, bit
  masks, groups, joins, missing values). It should probably lose its storage
  role -- see the design note.
- chisurf's `ParseModel` is still on Python `eval` and should stay there until
  the engine wins on FCS-shaped equations at realistic curve lengths.
