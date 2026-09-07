---
type: reference
title: "FPS export formats — everything the C# toolkit writes after a docking or screening run"
description: Byte-level specification of the six FPS outputs (per-result PyMOL script, Overlay, OverlayStates, R table, chi2 table, SimulationResults.bin), the arithmetic behind each number, and the defects a C++ port must decide to reproduce or fix.
resource: /Users/tpeulen/dev/imp.bff
tags: [reference, fps, export-format, pymol, docking, screening, rmsd, kabsch, porting, chisurf]
timestamp: '2026-08-31T00:00:00Z'
---

# FPS export formats

Everything below is derived from the C# at `/Users/tpeulen/dev/chisurf/junk/fps/Fps/`.
Every non-obvious claim cites `File.cs:line`. Read this page instead of the C#.

FPS has two modes, `FPSMode.Dock` and `FPSMode.Filter` (screening). The export
dialog (`SaveForm`) shows a different set of checkboxes for each
(`SaveForm.cs:22-28`) and every write goes through `okbutton_Click`
(`SaveForm.cs:159-197`).

## 0. The dialog and the file names

`MainForm.saveButton_Click` (`MainForm.cs:565-605`) fills the form. In Dock mode
it hands over `srs` (the whole result array, **by reference**), the selected
rows, the `MoleculeList`, the best-fit flag, the RMSD reference and the global
conversion polynomial. In Filter mode it hands over `frs` and the selected
`FilteringResult`s, and **never sets `BestFit` or `ConversionFunction`** — they
keep their defaults (`false` / `null`).

`DataToSaveAll`'s setter sorts the array in place by `InternalNumber`
(`SaveForm.cs:34-45`). Because `MainForm` passes `srs` itself, **opening the save
dialog re-sorts `MainForm.srs`**. Harmless in practice (it is already in that
order) but a real side effect. The same is true of `FilteringDataToSaveAll`
(`SaveForm.cs:62-73`).

Checkbox defaults: `etablecheckBox` and `pymolcheckBox` checked, radio button
`Selected` selected (`SaveForm.Designer.cs:72,84,164`). Choosing `Selected`
force-disables the binary export (`SaveForm.cs:506-516`) — **the `.bin` can only
ever contain all results, never a selection**.

The user picks a directory, not a file. `savepath` is that directory plus the
platform separator (`SaveForm.cs:176`). `filenametextBox.Text` is the prefix,
defaulting to `structure` in Dock mode and `screening_` in Filter mode
(`SaveForm.cs:28`).

| Mode | Checkbox | Output path |
|---|---|---|
| Dock | Pymol scripts | `<dir>/<prefix><InternalNumber>.pml`, one per result (`SaveForm.cs:182-183`) |
| Dock | Pymol saves pdb | adds `save <dir>/<prefix><N>.pdb` inside the `.pml` (`SaveForm.cs:224-228`) |
| Dock | Overlay | `<dir>/Overlay.pml` (`SaveForm.cs:184`) — prefix **not** used |
| Dock | Overlay (states) | `<dir>/OverlayStates.pml` (`SaveForm.cs:185`) — prefix **not** used |
| Dock | Model distances | `<dir>/Rtable_Rmp.txt` (+ `<dir>/Rtable_<DataType>.txt`) (`SaveForm.cs:186,398,423`) |
| Dock | Simulation results (bin) | `<dir>/SimulationResults.bin` (`SaveForm.cs:187`) |
| Dock | Chi2 table | `<dir>/chi2table.txt` (`SaveForm.cs:188`) |
| Filter | Model distances | `<dir>/<prefix>Rtable_Rmp.txt` (+ `<prefix>Rtable_<DataType>.txt`) (`SaveForm.cs:192`) |
| Filter | Chi2 table | `<dir>/<prefix>chi2table.txt` (`SaveForm.cs:193`) |

Note the inconsistency: `Overlay`, `OverlayStates`, `SimulationResults.bin` and
the Dock-mode `chi2table.txt`/`Rtable` ignore the user's prefix; the per-result
`.pml` and everything in Filter mode use it.

`<DataType>` is the `ToString()` of `DistanceDataType` — literally `Rmp`,
`RDAMean` or `RDAMeanE` (`Distances.cs:44-47`).

### Number formatting and locale

`Program.Main` clones the current culture and forces
`NumberFormat.NumberDecimalSeparator = "."` (`Program.cs:19-22`), so the decimal
point is always `.`. The `F` format specifier never inserts group separators, so
thousands separators cannot appear either. A C++ port should write everything in
the C locale and it will match.

`Vector3.ToString()` is `"{0:F3}, {1:F3}, {2:F3}"` — **three decimals, separated
by comma+space** (`MatrixVector3.cs:62-65`). Every `[x, y, z]` in every PyMOL
script uses exactly this.

`Boolean.ToString()` yields `True` / `False`. `SimulationMethods` is a `[Flags]`
enum (`SimulationResult.cs:5-13`), so its `ToString()` yields `Docking`,
`Unknown`, or a comma+space-joined list such as `Docking, Refinement`. **That
field can contain `, ` but never a tab**, so it is safe in a TSV, but a naive
comma-splitting reader will break.

## 1. Per-result PyMOL script — `SaveSimulationResult`

**Purpose.** Reconstruct one docking solution in PyMOL from the original,
untransformed PDB files, and optionally write the assembled complex to a PDB.

**Source.** `SaveForm.cs:199-262`.

**Format.** Plain text, `\r\n` line endings (`StreamWriter` default on Windows),
in this exact order:

1. One comment line:
   `# Energy = <E:F8> (converged)` or `# Energy = <E:F8> (not converged)`
   (`SaveForm.cs:201,206`). Note the leading space before the parenthesis is part
   of the `converged` string.
2. `load <FullFileName>` for **every** molecule, `i = 0 .. N-1`
   (`SaveForm.cs:207-208`). The path is unquoted, so a path containing a space
   breaks the script. `Molecule.FullFileName` is the path given at load time
   (`Molecule.cs:169`); `Molecule.Name` is `Path.GetFileNameWithoutExtension`
   (`Molecule.cs:199`), which is also the object name PyMOL will assign.
3. For `i = 1 .. N-1` (**molecule 0 is skipped**, `SaveForm.cs:209`), two lines:
   ```
   rotate [<ux:F3>, <uy:F3>, <uz:F3>], <theta:F3>, <Molecule.Name>, origin=[<cmx:F3>, <cmy:F3>, <cmz:F3>]
   translate [<tx:F3>, <ty:F3>, <tz:F3>], <Molecule.Name>
   ```
   `theta` is `Matrix3.AngleAndAxis(sr.Rotation[i], out u) * 180/pi`
   (`SaveForm.cs:211`); `u` is the unit rotation axis. `origin` is
   `Molecule.CM`, the **mass-weighted** centre of mass computed at PDB load
   (`Molecule.cs:222-223`). `translate` uses `sr.Translation[i]`.
4. If the best-fit checkbox is on (`SaveForm.cs:216-223`):
   ```
   select all
   translate [<BestFitTranslation:F3 x3>], sele
   rotate [<u:F3 x3>], <theta:F3>, sele, origin=[0, 0, 0]
   ```
   i.e. **translate first, then rotate about the origin** — that ordering is the
   contract, and `SimulationResult.RMSD` assumes it (`SimulationResult.cs:56-57`).
5. If "Pymol saves pdb" (`SaveForm.cs:224-228`):
   ```
   select all
   save <fname>.pdb, sele
   ```
   `<fname>` is the same absolute path prefix as the `.pml`, unquoted. This
   happens **before** any pseudoatoms are created, so the PDB contains only the
   molecules, post-best-fit.
6. If "Add labeling positions" (`SaveForm.cs:229-249`), for each labelling
   position in file order, four or five lines:
   ```
   pseudoatom <l.Name>, pos=[<rx:F3>, <ry:F3>, <rz:F3>]
   label <l.Name>, "<l.Name>"
   show spheres, <l.Name>
   color green, <l.Name>       # only if l.Dye == Donor
   color red, <l.Name>         # only if l.Dye == Acceptor
   ```
   with
   `r = sr.Rotation[im] * (lp_ref - molecules[im].CM) + molecules[im].CM + sr.Translation[im]`
   (`SaveForm.cs:241`), `im = molecules.FindIndex(l.Molecule)`. `lp_ref` is the
   **refined** labelling position from `sr.RefinedLabelingPositions` when the
   result carries the `Refinement` flag, otherwise the input position
   (`SaveForm.cs:238-240`). No colour line is written for `DyeType.Unknown`
   (`LabelingPositions.cs:7`).
   These pseudoatoms are placed at the **untransformed** (no best-fit) position.
7. If best-fit **and** at least one labelling position was written
   (`SaveForm.cs:250-258`):
   ```
   deselect
   select <name1> + <name2> + ...
   translate [<BestFitTranslation:F3 x3>], sele
   rotate [<u:F3 x3>], <theta:F3>, sele, origin=[0, 0, 0]
   ```
   The selection string is built as `" + " + name` repeated and then
   `.Substring(3)` to strip the leading `" + "` (`SaveForm.cs:247,253`) — correct.
   This is the second half of the split: molecules got the best-fit in step 4,
   pseudoatoms get it here, so the two end up consistent.
8. `deselect` (`SaveForm.cs:259`).

**Answer: is the best-fit applied to the written coordinates?** Yes, for the
PyMOL script — it is emitted as `translate`/`rotate` commands, so the geometry
PyMOL displays (and the optional `.pdb`, which is written after step 4) is
superposed. It is *not* baked into any number FPS itself computes for the tables.

### Defects in the per-result script

* **Molecule 0 is assumed to be identity.** The loop starts at `i = 1`
  (`SaveForm.cs:209`), which is valid because `SpringEngine` normalises every
  result so that `Translation[0] == 0` and `Rotation[0] == E`
  (`SpringEngine.cs:302-304`, mirrored in `MetropolisSampler.cs:54-57`). But
  `SaveOverlay`/`SaveOverlayStates` start at `i = 0` (`SaveForm.cs:281,326`). If
  a result ever carries a non-identity pose for molecule 0 — e.g. a hand-edited
  or foreign `.bin` — the per-result script and the overlays disagree. A port
  should always write all `N` transforms.
* **`rotate`/`translate` default to camera space.** PyMOL's `rotate` and
  `translate` take a `camera` argument that defaults to `1`, meaning the axis and
  the vector are interpreted in *camera* coordinates. FPS writes *model*-frame
  vectors. This works only because the camera rotation is identity in a fresh
  session right after `load`. Run the script into a session where the user has
  rotated the view and the structure lands wrong. **A C++ port should emit
  `camera=0` explicitly, or use `cmd.transform_selection` with a 4x4 matrix.**
* **Object-name collisions.** `Molecule.Name` is the bare basename. Two input
  PDBs with the same basename from different directories become the same PyMOL
  object name for the first, and `name_1` etc. for the rest — every subsequent
  `rotate`/`translate` targets the wrong object silently.
* **No quoting anywhere.** Paths with spaces break `load` and `save`.
* **Labelling-position names are used as PyMOL object names** with no
  sanitisation; a name containing `-`, `+` or a space corrupts both the
  `pseudoatom` line and the `select` list of step 7.
* `sr.BestFitRotation` is read but **never computed by this method** — see §7.

## 2. Overlay — `SaveOverlay`

**Purpose.** One PyMOL session containing every selected docking solution as a
**separate object**, so they can be shown/hidden and coloured individually.

**Source.** `SaveForm.cs:263-306`.

**Format.**

1. Header, written with `Write` for all but the last so it is one line:
   `# Overlay of structures <n0>, <n1>, ..., <nk>` (`SaveForm.cs:269-272`).
   These are `InternalNumber`s. **Crashes if the selection is empty**
   (`_dataToSave[_dataToSave.Length - 1]` with length 0).
2. `load <FullFileName>` for every molecule (`SaveForm.cs:274-275`).
3. For each result `j`, and for **each** molecule `i = 0 .. N-1`
   (`SaveForm.cs:281-288`):
   ```
   copy _tmp<i>, <Molecule.Name>
   rotate [<u:F3 x3>], <theta:F3>, _tmp<i>, origin=[<CM:F3 x3>]
   translate [<T:F3 x3>], _tmp<i>
   ```
   then
   ```
   select object _tmp*
   create _tmpjoin, sele
   ```
   then, if best-fit,
   ```
   translate [<BestFitTranslation:F3 x3>], _tmpjoin
   rotate [<u:F3 x3>], <theta:F3>, _tmpjoin, origin=[0, 0, 0]
   ```
   then
   ```
   copy <prefix><InternalNumber>, _tmpjoin
   delete _tmp*
   ```
4. After the loop, `delete <Molecule.Name>` for every molecule
   (`SaveForm.cs:301-302`) — the pristine originals are removed.

**Result:** `k` objects named `<prefix><InternalNumber>`, each a single-state
merge of all molecules. The user's prefix *is* used here (`SaveForm.cs:298`),
even though the file name is not.

**Defects.** `select object _tmp*` — `object` is not a PyMOL selection operator
in the way `model`/`m.` is; `select _tmp*` or `select model _tmp*` is the
canonical spelling. Verify against your PyMOL version before porting verbatim.
`delete _tmp*` also matches `_tmpjoin`, which is intended. The camera-space
issue of §1 applies identically.

## 3. OverlayStates — `SaveOverlayStates`

**Purpose.** The same overlay, but as **one object with N states**, for morphing
or for state-by-state playback.

**Source.** `SaveForm.cs:308-351`. Header, `load` block and the per-molecule
`copy`/`rotate`/`translate` block are **byte-identical** to `SaveOverlay`
(compare `SaveForm.cs:314-332` with `263-287`).

The difference is what happens after the transforms (`SaveForm.cs:334-344`):

```
select object _tmp*
translate [<BestFitTranslation:F3 x3>], sele     # only if best-fit
rotate [<u:F3 x3>], <theta:F3>, sele, origin=[0, 0, 0]   # only if best-fit
save _tmp.pdb, sele
load _tmp.pdb, <prefix>
delete _tmp*
```

There is no `create _tmpjoin`; the merge happens implicitly by saving the whole
selection to one PDB and loading it back into the object named `<prefix>`, which
appends a state each time.

**Differences from Overlay, in one line:** Overlay -> `k` single-state objects
named by `InternalNumber`; OverlayStates -> one object named `<prefix>` with `k`
states numbered `1..k` in `InternalNumber` order.

**Defects.**

* **`_tmp.pdb` is a relative path.** It is written into PyMOL's current working
  directory, not the export directory, and is **left behind** when the script
  finishes. If the cwd is not writable the whole script fails silently mid-way.
* **The `InternalNumber` of each state is lost.** States are 1..k; the only
  record of which state is which result is the `#` header line, and only if the
  reader knows the states follow the same order. A port should write a
  `set_title` / `alter` line per state.
* `delete _tmp*` would delete the accumulating object if the user set the prefix
  to something starting with `_tmp`.
* Round-tripping through PDB truncates coordinates to PDB's `%8.3f` and merges
  chains — the overlay object is not a faithful copy of the inputs.

## 4. The R table — `SaveDistances` (Dock) and `SaveFilterDistances` (Filter)

**Purpose.** For every structure, the model distance for every distance
restraint in the project — the table you regress against experiment.

### 4a. Dock mode — `SaveDistances` (`SaveForm.cs:392-447`)

Always writes `<base>_Rmp.txt`. Tab-separated, one header row:

```
Structure<TAB>Number<TAB><lp1>_<lp2><TAB><lp1>_<lp2><TAB>...
```

The column name is `l1.Name + '_' + l2.Name` for each entry of the project's
`DistanceList`, in file order (`SaveForm.cs:400-409`). **The names are joined
with `_` and are not escaped**, so a labelling position whose own name contains
`_` makes the header ambiguous.

One row per result (`SaveForm.cs:410-416`):

```
<prefix><InternalNumber><TAB><InternalNumber><TAB><R:F3><TAB><R:F3>...
```

The `Structure` and `Number` columns are redundant — both derived from
`InternalNumber`.

`R` is `SimulationResult.ModelDistance(l1, l2)` (`SimulationResult.cs:167-186`):

```
d = R[im1]*(lp1_ref - CM[im1]) + CM[im1] + T[im1]
  - R[im2]*(lp2_ref - CM[im2]) - CM[im2] - T[im2]
R = |d|
```

with `lp_ref` taken from `RefinedLabelingPositions` when the result carries the
`Refinement` flag and that list is non-null (`SimulationResult.cs:173-177`).
This is a **mean-position distance R_mp** and nothing else — in Dock mode a
labelling position is a single point, so there is no AV cloud to average over.

If the project's distance data type is not `Rmp`, a second file
`<base>_<DataType>.txt` is written with the same header and row layout
(`SaveForm.cs:421-446`), but each value is `Distance.Convert(R_mp, _globalcf)` —
a polynomial evaluated by Horner from the highest coefficient down
(`Distances.cs:36-41`) — **only when both labelling positions have a known dye
type**; otherwise the raw `R_mp` is written unchanged (`SaveForm.cs:438-439`).

**Defects.**

* The second file's values are `R_mp` pushed through a *global* R_mp->R_DA
  polynomial. It is **not** an AV-based <R_DA> or <R_DA>_E. A column headed
  `RDAMeanE` in a Dock-mode table is a polynomial estimate, not a simulated one.
  Filter mode's file of the same name *is* simulated (§4b) — **the same file name
  means two different quantities in the two modes.**
* If `dist.DataType != Rmp` and `_globalcf` is `null`, `Distance.Convert`
  dereferences a null array. `MainForm` sets `globalcf` to `null` in several
  reset paths (`MainForm.cs:349,375,505,1488`).
* The second loop re-looks-up `l1`/`l2` (`SaveForm.cs:428-429`) but only uses
  them for the header; the values still use the `l1index`/`l2index` arrays
  captured in the first loop. Correct, but dead code.

### 4b. Filter mode — `SaveFilterDistances` (`SaveForm.cs:448-494`)

Same two files, same column names, but the first two columns are

```
File<TAB>Number
```

and the row is `<fr.ShortFileName><TAB><j>` where **`j` is the zero-based array
index, not `fr.InternalNumber`** (`SaveForm.cs:463,486`).

> **This is the off-by-one.** `InternalNumber` is 1-based —
> `MainForm.saveButton_Click` indexes with `InternalNumber - 1`
> (`MainForm.cs:580,596`), and the Dock-mode R table and the chi2 table both
> print `InternalNumber`. The Filter R table prints `j`. So the `Number` column
> of `screening_Rtable_Rmp.txt` is offset by one from the `Number` column of any
> Dock table and from the identity used everywhere else in FPS. Worse, when the
> user exported a *selection*, `j` numbers the rows of the selection, not the
> screening run — it is not even a stable identifier. **Join on `File`, never on
> `Number`.**

The values:

* `<base>_Rmp.txt` writes `fr.RmpModel[i].R` — the true mean-position distance,
  `|rmp[l1] - rmp[l2]|` (`FilterEngine.cs:304`).
* `<base>_<DataType>.txt` writes `fr.RModel[i].R` — the distance of the
  project's data type: `R_mp` if either position has no AV
  (`FilterEngine.cs:305-307`), else a Monte-Carlo `RdaMeanEFromAv` or
  `RdaMeanFromAv` over the two AV clouds (`FilterEngine.cs:308-319`). This can be
  `Double.NaN` when an AV came out empty (`FilterEngine.cs:310,316`), which
  `ToString("F3")` renders as the literal `NaN`.

**Defect:** `FilterEngine` sets `fr.RmpModel.DataType = dist.DataType`
(`FilterEngine.cs:293`) — the *mean-position* list is tagged with the
*experimental* data type. It is only in-memory metadata and does not reach the
file, but any port that carries the tag forward will mislabel the column.

**Defect:** the two files are indexed positionally (`RmpModel[i]`, `RModel[i]`)
against `dist[i]`, while the header names come from `labelingpos.Find(...)`. That
is consistent only because `FilterEngine` builds both lists in `dist` order
(`FilterEngine.cs:294-322`). Preserve that invariant.

## 5. The chi2 table

### 5a. Dock mode — `SaveEnergyTable` (`SaveForm.cs:353-375`)

**Reference selection** (`SaveForm.cs:359`): `_rmsdref = _vsref ? _rmsdref : _dataToSave[0];`
`_vsref` is `rmsdvsref & (sr_rmsdreference.InternalNumber > 0)` (`MainForm.cs:587`).
When the user has not picked a reference in the GUI, the reference becomes **the
first exported result** — after the sort, the lowest `InternalNumber` in the
exported set, not necessarily the best-scoring one. The assignment mutates the
form's field.

**Header** (`SaveForm.cs:361-362`), tab-separated, 9 columns:

```
File	chi2	chi2_bond	chi2_clash	Convergence	Method	RMSD vs previous	RMSD vs <M>	RMSD vs <M> (selected)
```

where `<M>` is `_rmsdref.InternalNumber`. Columns 8 and 9 have **the same header
text** apart from the ` (selected)` suffix.

**First data row** (`SaveForm.cs:364-365`), format
`{0}{1}\t{2:F8}\t{3:F4}\t{4:F4}\t{5}\t{6}\t---\t{7:F6}\t{8:F6}`:

| Col | Content | Format |
|---|---|---|
| File | `<prefix><InternalNumber>` | — |
| chi2 | `sr.E` | `F8` |
| chi2_bond | `sr.Ebond` | `F4` |
| chi2_clash | `sr.Eclash` | `F4` |
| Convergence | `sr.Converged` | `True`/`False` |
| Method | `sr.SimulationMethod` | flags enum, e.g. `Docking, Refinement` |
| RMSD vs previous | literal `---` | — |
| RMSD vs M | `sr.RMSD(ref, bestfit)` | `F6` |
| RMSD vs M (selected) | `sr.RMSD(ref, bestfit, true)` | `F6` |

**Subsequent rows** (`SaveForm.cs:366-372`): identical, except column 7 is
`sr.RMSD(_dataToSave[i-1], bestfit)`. "Previous" means the previous row of the
exported, `InternalNumber`-sorted array — **not** the previous in simulation
time, and **not** the previous in the GUI's sort order.

**What the columns actually are.**

* **`E`** is the spring energy from the last `Iterate` (`SpringEngine.cs:259,281,292`),
  a sum over restraints of a quadratic penalty going linear past `MaxForce`
  (`SpringEngine.cs:432-444`). With `k = 2/sigma^2` and `dE = 0.5*k*dr^2`, `E`
  equals `sum (dr/sigma)^2` — a raw **chi-square, not a reduced one**
  (`SimulationResult.cs:19` says so).
  * **The header lies, and it disagrees with the GUI.** The grid divides by
    `dof = max(Ndist - 6*(Nmolecules-1), 1)` and shows `srs[i].E / dof`
    (`MainForm.cs:683,701,719`); `chi2table.txt` writes the **undivided** `E`.
    A port must decide which it means and say so; do not name the raw column
    `chi2`.
  * `E` **excludes clash energy** — `Eclash` is separate
    (`SpringEngine.cs:465,548`), added only in `SetState`
    (`SpringEngine.cs:209-211,325-327`), never into `sr.E`.
  * Under `OptimizeSelected == Selected`, `ForceAndTorque` skips unselected
    restraints (`SpringEngine.cs:428`), so **`E` is a partial sum** and is not
    comparable across projects. `SelectedThenAll`'s final pass covers everything.
* **`Ebond`** is the subset of the same `dE` terms whose `Distance.IsBond` is set
  (`SpringEngine.cs:445`, reset at `:410`). A "bond" is decided in
  `PrepareSimulation` (`SpringEngine.cs:111-116`): both labelling positions must
  have `Dye == Unknown`, `AVData.AtomID > 0` and `AVType == None` — a distance
  between two *atoms*, not two dyes. As a side effect both anchor atoms get
  `AtomData.vdWRNoClash` so they are excluded from clash detection
  (`SpringEngine.cs:118-134`). **`Ebond` is non-zero exactly when the project has
  an atom-atom restraint and it is violated**; in a pure-FRET project it is
  identically `0.0000`. It is a **subset** of `E`, not an additional term — do
  not add them.
* **`Eclash`** is the clash penalty (`SpringEngine.cs:465,548`) with
  `kclash = 2/ClashTolerance^2` (`:147`). Headed `chi2_clash`; **not** in `E`.
* **`Converged`** is `niter < MaxIterations` (`SpringEngine.cs:295`). Under
  `SelectedThenAll` the cap is `MaxIterations/2` (`:254-255`) and `niter` resets
  before the second pass (`:277`) — **so `Converged` is `True` unconditionally in
  that mode.** A real bug; the flag is meaningless there.
* **`RMSD`** — see §6. **Unweighted**, all atoms of all molecules, unless the
  `(selected)` variant restricts to `Molecule.Selected` molecules
  (`SimulationResult.cs:52`, `Molecule.cs:114`).

### 5b. Filter mode — `SaveFilterTable` (`SaveForm.cs:376-390`)

**Header** (`SaveForm.cs:381`), 7 tab-separated columns:

```
File	Chi2r	NaNs	RefRMSD	>1sigma	>2sigma	>3sigma
```

| Col | Content | Format |
|---|---|---|
| File | `fr.ShortFileName` (basename) | — |
| Chi2r | `fr.E` | `F3` |
| NaNs | `fr.InvalidR` | integer |
| RefRMSD | `fr.RefRMSD` | `F3` |
| >1sigma | `fr.Sigma1` | integer |
| >2sigma | `fr.Sigma2` | integer |
| >3sigma | `fr.Sigma3` | integer |

**Arithmetic** (`FilterEngine.cs:294-341`), with `dr = R_model - R_exp`:

```
E      += dr>0 ? dr^2/ErrPlus^2 : dr^2/ErrMinus^2   # only if !OptimizeSelected || IsSelected
activeR++
Sigma1 += (dr >  ErrPlus)  + (dr <  -ErrMinus)
Sigma2 += (dr > 2*ErrPlus) + (dr < -2*ErrMinus)
Sigma3 += (dr > 3*ErrPlus) + (dr < -3*ErrMinus)
E       = E / activeR
RefRMSD = nref_total == 0 ? 0 : sqrt(refrmsd_t / nref_total)
```

* `Chi2r` **is** reduced here — divided by `activeR` — but **not** by a
  degrees-of-freedom count; no parameters are subtracted. It is a mean squared
  normalised residual, and therefore *not* comparable with the Dock-mode `chi2`.
* If `activeR == 0`, `E = NaN` (`FilterEngine.cs:340`).
* `E` can be the sentinel `-999.9` (`MiscData.ENotCalculated`,
  `StaticData.cs:152`) when the molecule failed to load
  (`FilterEngine.cs:184`); the GUI hides it (`MainForm.cs:771-774`) but **the
  file writes `-999.900` as if it were data.**
* The sigma counters are **cumulative, not exclusive**: a 3.5-sigma outlier
  increments all three, so `Sigma1 >= Sigma2 >= Sigma3` and "exactly 1-2 sigma"
  is `Sigma1 - Sigma2`.
* NaN distances are counted in `InvalidR` and **excluded** from `E`, `activeR`
  and the sigma counters (`FilterEngine.cs:325-326`). A structure with many empty
  AVs therefore gets a flatteringly low `Chi2r` from few restraints. **Always
  read `NaNs` alongside `Chi2r`.**
* `RefRMSD` is **not** a structure-to-structure RMSD. It is the RMSD of the
  *reference-atom frames* used to place non-AV labelling positions onto the
  candidate (`FilterEngine.cs:225-278`) — a sanity check on the local geometry
  around the label, over all `AVSimlationType.None` positions. `0.0` when no
  position uses reference atoms.
* `FilteringResult.BestFitRotation` and the structure-vs-structure
  `FilteringResult.RMSD` (`SimulationResult.cs:221-242`) are computed for the GUI
  only (`MainForm.cs:792-793`) and **never written to any file**.

## 6. RMSD — what it actually measures

`SimulationResult.RMSD(other, bestfit, selected_only)` (`SimulationResult.cs:40-72`):

```
if (bestfit) this.CalculateBestFitRotation(sr_other, false);   // line 47
for each molecule i (skip if selected_only && !m.Selected):
    if bestfit:
        U = this.BestFitRotation * this.Rotation[i] - other.Rotation[i];
        t = this.BestFitRotation * (this.Translation[i] + m.CM + this.BestFitTranslation)
            - other.Translation[i] - m.CM - other.BestFitTranslation;
    else:
        U = this.Rotation[i] - other.Rotation[i];
        t = this.Translation[i] - other.Translation[i];
    for each atom j:
        r = (m.XLocal[j], m.YLocal[j], m.ZLocal[j]);   // CM-subtracted local coords
        sd += SquareNormDiff(U * r, t);                 // = |U*r - t|^2
    NAtomstotal += m.NAtoms;
return sqrt(sd / NAtomstotal);
```

**Reference:** the `other` argument — the previous exported row for column 7, the
`_rmsdref` structure for columns 8 and 9.

**Weighting:** none. Every atom counts 1. `CalculateBestFitRotation` is called
with `weighted = false` (`SimulationResult.cs:47`), and the mass-weighted branch
of `CalculateBestFitTranslation` (`:92-96`) is never exercised from the export
path — `SpringEngine.cs:307` calls `CalculateBestFitTranslation(false)`, the
unweighted centroid-to-origin branch (`:97-106`). The weighted code paths are
dead.

**Atom set:** all atoms of all molecules, or all atoms of `Selected` molecules.
Not CA-only, not by-chain.

> ### The sign bug
>
> The per-atom displacement between the two structures is
> `(R1*r + CM + T1) - (R2*r + CM + T2) = (R1-R2)*r + (T1-T2) = U*r + t`.
> The code accumulates `SquareNormDiff(U*r, t)`, and
> `SquareNormDiff(a, b) = |a - b|^2` (`MatrixVector3.cs:51-55`).
> **So FPS computes `|U*r - t|^2` where the correct quantity is `|U*r + t|^2`.**
> The translation difference enters with the wrong sign. The two agree only when
> `t == 0` or when `U*r` is perpendicular to `t` for every atom. The same
> inversion is in the best-fit branch (`SimulationResult.cs:55-57` vs line 67).
>
> Every RMSD FPS reports — the three RMSD columns of `chi2table.txt` and the
> GUI's RMSD column — is affected. A port that fixes it will not reproduce FPS's
> numbers; a port that reproduces them is propagating the bug. **Decide
> explicitly and record the decision.** `FilteringResult.RMSD`
> (`SimulationResult.cs:235-241`) does *not* have the bug — it compares
> coordinates directly.

`CalculateBestFitRotation` (`SimulationResult.cs:116-159`) is a standard Kabsch
via SVD of `R = Ry * Rxt`, with the determinant sign flip on the third singular
value (`:154`) and a final `Matrix3.RepairRotation` that re-orthogonalises
through axis-angle (`MatrixVector3.cs:147-154`). It short-circuits to the
identity if either `E` is `NaN` (`:118-122`).

`CalculateBestFitTranslation(false)` (`:97-106`) puts the **unweighted centroid
of all atoms of all molecules at the origin**. It is called once, at the end of
docking (`SpringEngine.cs:307`), so **`BestFitTranslation` is a property of the
structure alone, not of any pair.**

## 7. Where `BestFitRotation` comes from — the biggest trap

`SaveSimulationResult`, `SaveOverlay` and `SaveOverlayStates` all read
`sr.BestFitRotation` (`SaveForm.cs:220,294,338`) but **none of them computes it**.

It is computed as a side effect of `SimulationResult.RMSD`
(`SimulationResult.cs:47`), called from `MainForm.DisplayNewStructures` on the
live array elements (`MainForm.cs:694,708,718` — `srs[n].RMSD(...)` mutates
`srs[n]` because array element access on a value type yields a reference).

1. **The best-fit rotation baked into an exported `.pml` is whatever the GUI last
   computed.** With the default "RMSD vs previous" (`MainForm.cs:742-745`),
   `srs[i].BestFitRotation` is the fit of structure *i* onto structure *i-1* —
   **not** onto a common reference. The "Overlay" is then a chain of pairwise
   fits and the structures are *not* mutually superposed. Only after the user
   picks a reference (`MainForm.cs:1459-1460`) does every `BestFitRotation` refer
   to the same one.
2. **If the best-fit checkbox was off when the grid was last refreshed,
   `BestFitRotation` is the default-constructed `Matrix3` — all zeros, not the
   identity.** `AngleAndAxis` of the zero matrix gives `atan2(0,-1) = pi` about
   `(0,0,1)` (`MatrixVector3.cs:129-141`), so the script emits a spurious 180
   degree rotation.
3. **`SaveEnergyTable`'s RMSD calls do not fix this.** `tmp = _dataToSave[i]`
   copies the struct (`SaveForm.cs:363,368`), so the `BestFitRotation` computed
   inside `tmp.RMSD(...)` is discarded — and the `.pml` files are written first
   (`SaveForm.cs:180-188`) anyway.

**Answer to "is the best-fit applied or only reported?"**

* **PyMOL script / Overlay / OverlayStates:** applied, as emitted
  `translate`+`rotate` commands, using the stale `BestFitRotation` above and the
  structure's own `BestFitTranslation`. The optional `.pdb` is saved after the
  transform, so its coordinates are superposed.
* **R table (both modes):** **not applied.** `ModelDistance` uses only
  `Rotation`/`Translation`/`CM` (`SimulationResult.cs:183-184`); distances are
  rotation-invariant anyway.
* **chi2 table:** the best fit is used *inside* the RMSD computation only, and is
  recomputed fresh for each column's reference (`SimulationResult.cs:47`), so
  those three numbers are self-consistent even though the `.pml` files are not.
  Nothing about the best fit is written to the chi2 table itself.

## 8. The binary results file — `SaveSimulationResuntsBin`

Note the typo in the method name: **`SaveSimulationResunts`** (`SaveForm.cs:496`).

```csharp
BinaryFormatter bf = new BinaryFormatter();
using (FileStream fs = new FileStream(srfilename, FileMode.Create))
{ bf.Serialize(fs, _dataToSave); fs.Close(); }
```

**What it serializes:** the whole `SimulationResult[]`. `SimulationResult` is
`[Serializable]` (`SimulationResult.cs:16`), so every field is written **except
`Molecules`, which is `[NonSerialized]`** (`:31-32`). The file carries `E`,
`Eclash`, `Ebond`, `Translation[]`, `Rotation[]`, `BestFitTranslation`,
`BestFitRotation`, `Converged`, `SimulationMethod`, `InternalNumber`,
`ParentStructure` and `RefinedLabelingPositions` — and **not** the molecules, so
it is meaningless without the matching project. The reload path re-attaches the
currently-open `MoleculeList` **by position** (`MainForm.cs:467-468`) with **no
check that they are the same molecules in the same order**. Load a `.bin`
against a different project and you silently get garbage.

Because the `Selected` radio button disables the checkbox
(`SaveForm.cs:506-516`), this file always contains every result.

**Portability: none.** This is .NET `BinaryFormatter` — MS-NRBF, embedding
fully-qualified .NET type names (`Fps.SimulationResult`, `Fps.Vector3`,
`Fps.Matrix3`, `Fps.LabelingPosition`, `Fps.LabelingPositionList`) and the
assembly identity of `Fps.exe`, together with field names and an object graph.
Reading it from C++ would mean implementing MS-NRBF and hard-coding the field
layout of a long-dead assembly build. **A C++ reader is impractical and should
not be attempted.** `BinaryFormatter` is also disabled by default since .NET 5
and removed in .NET 9, so even a C# reader needs a legacy target.

**Recommended replacement in the port:** JSON or cereal, carrying exactly the
fields listed above plus a molecule manifest (file name, atom count, CM) so the
"re-attach by position" hazard becomes a checkable assertion.

## 9. What is lossy to regenerate from `fps.json` + rigid-body poses

**Not recoverable at all**

* **`Converged`** — depends on how many iterations the optimiser took versus
  `MaxIterations` at run time (`SpringEngine.cs:295`). Nothing in the geometry
  encodes it. (And it is always `True` in `SelectedThenAll` anyway — §5a.)
* **`SimulationMethod`** — the provenance flags `Docking | Refinement |
  ErrorEstimation | MetropolisSampling` (`SimulationResult.cs:5-13`). This is
  what tells you a row is an error-estimate replica rather than an independent
  solution; the GUI colours rows by it (`MainForm.cs:720-732`).
* **`ParentStructure`** — which structure a refinement or error-estimation run
  branched from (`SimulationResult.cs:29`, set at `SpringEngine.cs:298`). The
  tree structure of a run lives only here.
* **`InternalNumber`** — the stable 1-based identity used in every file name and
  every `Number`/`File` column.
* **`RefinedLabelingPositions`** — the refined mean positions from
  `Refinement.RedoAV` (`Refinement.cs:37`). `ModelDistance` and the `pseudoatom`
  lines both prefer them (`SimulationResult.cs:173-177`, `SaveForm.cs:238-240`).
  **For any result carrying the `Refinement` flag, the R table and the PyMOL
  labelling positions cannot be regenerated from `fps.json` alone.** The single
  most important lossy item in Dock mode.
* **Filter-mode `RModel` for `RDAMean`/`RDAMeanE`** — Monte-Carlo over AV clouds
  seeded with `rnd.Next()` (`FilterEngine.cs:311-318`); the seed is not recorded.
  `InvalidR`, `Sigma1..3`, `E` and `RefRMSD` all depend on them.
* **`FilteringResult.RefRMSD`** — needs the reference-atom lists and the
  candidate PDB, not just poses (`FilterEngine.cs:225-278`).

**Recoverable only with extra information**

* **`E`, `Eclash`, `Ebond`** — recomputable, but only against the same
  `SimulationParameters` (`MaxForce`, `ClashTolerance`), the same
  `OptimizeSelected` mode, the same `IsBond` classification *and* the same
  vdW-radius suppression of bond anchors (`SpringEngine.cs:118-134`). None of
  that is in the poses.
* **`BestFitTranslation`** — recomputable exactly; the unweighted centroid of all
  atoms (`SimulationResult.cs:97-106`).
* **`BestFitRotation`** — recomputable by Kabsch, but you must know *which*
  reference the stored value was fitted against, and §7 shows FPS does not record
  that. Treat the stored value as unreliable.

**Fully recoverable**: everything in the PyMOL scripts except the best-fit block
and the `# Energy` header; both R-table files in Dock mode for non-refined
results; the `File`/`Number` columns.
