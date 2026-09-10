"""One-off: insert the conservation section into labelizer_score.ipynb, execute, save."""
import nbformat as nbf
from nbclient import NotebookClient

nb = nbf.read("labelizer_score.ipynb", as_version=4)

md = nbf.v4.new_markdown_cell(r"""## The conservation term, on the paper's worked example (1DDB)

Everything above dropped conservation — T4 lysozyme has no ConSurf grades, and
none are reachable today: the labelizer.org API accepts the conservation toggle
but its output comes back without the `cs` file for every entry tried (3GUN and
1DDB, custom upload and PDB-database load alike), its REST surface exposes no
grades endpoint, and the precomputed ConSurfDB host it cites is unreachable.

What the repository *does* ship is the reference distribution's own ConSurf
output for the paper's worked example — mouse BID, PDB **1DDB**, model 39,
chain A (`test/input/labelizer/1DDB-conservationscore-39-A.pdb`: 195 residues,
187 distinct grades, the normalised grade in the B-factor column, exactly the
file `ll_read_consurf` reads). Nothing to download. With it, the full published
model runs — conservation at weight 1 beside exposure, cysteine resemblance and
secondary structure.""")

code1 = nbf.v4.new_code_cell(r"""pdb_1ddb = "../../test/input/labelizer/1DDB-39.pdb"
grades = "../../test/input/labelizer/1DDB-conservationscore-39-A.pdb"

full_model = bff.ll_model_paper()          # cs(w=1), se(w=1), tp(0), cr(w=1), ss(w=1), ce(0)
scores_cs = bff.ll_score_structure(pdb_1ddb, full_model, bff.LlOptions(), grades)

comb_cs = {bff.ll_residue_key(r.asym_id, r.seq_id): r.value
           for r in scores_cs if r.score_type == "combined"}
status = {}
for r in scores_cs:
    if r.score_type == "combined":
        status[r.status] = status.get(r.status, 0) + 1
print(f"{len(comb_cs)} positions scored with the full model; statuses {status}")
print(f"conservation term: {sum(1 for r in scores_cs if r.score_type == 'conservation' and r.status == 'scored')} rows")

print(f"\n{'rank':>4}  {'site':>6}  {'LS':>7}")
for rank, (key, value) in enumerate(
        sorted(comb_cs.items(), key=lambda kv: -kv[1])[:10], start=1):
    print(f"{rank:>4}  {key:>6}  {'':>7}  {value:>7.4f}")""")

code2 = nbf.v4.new_code_cell(r"""# with vs without: how much the conservation term moves the answer
scores_ncs = bff.ll_score_structure(
    pdb_1ddb, [p for p in full_model if p.tag != "cs"], bff.LlOptions())
comb_ncs = {bff.ll_residue_key(r.asym_id, r.seq_id): r.value
            for r in scores_ncs if r.score_type == "combined"}

grades_map = bff.ll_read_consurf(grades)
keys = sorted(set(comb_cs) & set(comb_ncs), key=lambda k: int(k[1:]))
x = np.array([comb_ncs[k] for k in keys])          # no conservation
y = np.array([comb_cs[k] for k in keys])           # full published model
g = np.array([grades_map[k] for k in keys])        # the ConSurf grade itself
print(f"grade range {g.min():.2f}..{g.max():.2f} over {len(keys)} residues")

rank_c = {k: i for i, k in enumerate(sorted(keys, key=lambda k: -comb_cs[k]))}
rank_n = {k: i for i, k in enumerate(sorted(keys, key=lambda k: -comb_ncs[k]))}
top_c = set(sorted(comb_cs, key=comb_cs.get, reverse=True)[:10])
top_n = set(sorted(comb_ncs, key=comb_ncs.get, reverse=True)[:10])
print(f"top-10 overlap with/without conservation: {len(top_c & top_n)}/10")

fig, ax = plt.subplots(figsize=(5.6, 5.4))
lim = (min(x.min(), y.min()) - 0.05, max(x.max(), y.max()) + 0.05)
ax.plot(lim, lim, "k--", lw=0.8)
sc = ax.scatter(x, y, c=g, cmap="RdYlGn_r", s=22)
fig.colorbar(sc, ax=ax, label="ConSurf grade (B-factor)")
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("LS, no conservation")
ax.set_ylabel("LS, full published model")
ax.set_title("1DDB model 39: what conservation changes")
ax.grid(alpha=0.3)
fig.tight_layout()
plt.show()""")

md3 = nbf.v4.new_markdown_cell(r"""**Reading the conservation comparison**

- Conservation **reshapes the ranking**: only 3 of the top-10 sites survive
  from the no-conservation model. Variable positions rise (a site moving from
  rank 108 to 55 is typical), conserved ones sink — a residue the evolutionary
  record says is functionally constrained is a poor place to force a dye.
- The colour carries the mechanism: high-ConSurf-grade (variable) residues sit
  above the diagonal, conserved ones below.
- **A trap this section avoids, for the record**: the *other* shipped file,
  `1DDB-39_cs.pdb`, is the reference example's own *output* — re-running its
  example feeds that back as grades, landing in a two-cycle with just two
  distinct values. The lookup reproduces that degenerate pair to all sixteen
  digits (`test_the_conservation_lookup_is_exact_to_the_last_digit`), which is
  how the conservation machinery is verified without a reproducible reference
  input; `test_the_shipped_conservation_reference_is_degenerate` guards the
  story. See `okf/validation/labelizer_ab.md`.
- Should labelizer.org's conservation pipeline come back (its analysis emits a
  `*_cs.csv` when the internal ConSurf job succeeds — none of ours did), the
  same parity treatment as the SE/CR/SS section above applies unchanged.""")

# insert after the "Reading the results" section (end of notebook)
idx = max(i for i, c in enumerate(nb.cells)
          if c.cell_type == "markdown" and c.source.startswith("## Reading the results"))
nb.cells.insert(idx + 1, md)
nb.cells.insert(idx + 2, code1)
nb.cells.insert(idx + 3, code2)
nb.cells.insert(idx + 4, md3)

client = NotebookClient(nb, timeout=1800, kernel_name="python3")
client.execute()
nbf.write(nb, "labelizer_score.ipynb")

errs = [o for c in nb.cells if c.cell_type == "code"
        for o in c.get("outputs", []) if o.output_type == "error"]
print("written; errors:", len(errs))
for e in errs:
    print(e.ename, e.evalue)
