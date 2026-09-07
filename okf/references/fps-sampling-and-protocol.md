---
okf_version: "0.2"
type: Reference
title: "FPS error estimation, Metropolis sampling and docking protocol — the C# source, specified for a C++ port"
description: "Line-by-line specification of the FPS toolkit's parametric-bootstrap error estimation, its Metropolis sampler, and every protocol knob that changes a docked answer (OptimizeSelected, the SetState shuffle, MaxForce, ClashTolerance, adaptive dt/viscosity, the 'bond' rule), with the published default parameters for all five modes and the defects found."
tags: [reference, fps, docking, error-estimation, metropolis, bootstrap, port, prd-121]
timestamp: '2026-08-31T00:00:00Z'
---

# FPS error estimation, Metropolis sampling and docking protocol

Source of truth: the FPS C# tree at `/Users/tpeulen/dev/chisurf/junk/fps/Fps/`.
All file:line citations refer to that tree. This page specifies the
**quantities** FPS computes, so an IMP-based C++ port can reproduce them even
though the optimiser underneath will differ (FPS uses a damped Verlet
rigid-body integrator; IMP uses its own movers and optimisers).

The one thing a port must preserve exactly is the **restraint algebra** and the
**protocol** — what is perturbed, what is randomised, what is held fixed. The
integrator is an implementation detail; the perturbation model is not.

## 0. The energy — what every mode is minimising

`SpringEngine.ForceAndTorque()` (`SpringEngine.cs:407-456`) is the whole
restraint term. For each distance `i` with target `R`, asymmetric errors
`ErrPlus`/`ErrMinus`, and model distance `absr`:

```
dr = absr - R
k+ = 2 / ErrPlus^2            (SpringEngine.cs:137)
k- = 2 / ErrMinus^2           (SpringEngine.cs:138)
```

In the harmonic region `E_i = 0.5*k*dr^2 = (dr/Err)^2` — **the energy is a
chi-square, not a half-chi-square.** Sign convention: `dr > 0` (model too long)
uses `ErrPlus`; `dr < 0` uses `ErrMinus` (`SpringEngine.cs:432-443`).

The UI divides by `dof = max(Ndistances - 6*(Nmolecules-1), 1)`
(`MainForm.cs:683`). `SimulationResult.E` itself is **not** reduced
(`SimulationResult.cs:18`). `E` as stored *includes* `Eclash` for the docking
path (`SpringEngine.cs:325-327`, becoming `sr.E` at `:292`).

## 1. Error estimation — a parametric bootstrap

**Class:** `ErrorEstimation : SpringEngine` (`ErrorEstimation.cs`), 62 lines.
Header comment at `:5`: `// Error estimation by bootstrapping`.

### 1.1 The algorithm

Per repetition (`SimulationJobManager.cs:96-104`):

1. **`SetState(sr)` — install the "truth".** `ErrorEstimation.SetState`
   (`:17-33`) **overwrites the target `R` of every distance** with the *model*
   distance of the parent docked structure:
   ```csharp
   rmodeltmp.R = sr.ModelDistance(LabelingPositions.Find(p1),
                                  LabelingPositions.Find(p2));   // :26-27
   ```
   The experimental distances are discarded for the run. The "truth" is **the
   docked model itself**, not the measurement. `ModelDistance`
   (`SimulationResult.cs:167-185`) uses `RefinedLabelingPositions` if the parent
   came out of refinement. The loop at `:23-29` runs over **all** distances,
   including deselected ones. It then calls `base.SetState(sr)`
   (`SpringEngine.cs:200-212`), placing every molecule at the parent pose; the
   returned energy is ~`Eclash` only, since every `dr` is 0 by construction.

2. **`Simulate()` — perturb once, then re-optimise.** (`:35-41`)
   `PerturbDistances(1.0)` then `base.Simulate()`, ORing
   `SimulationMethods.ErrorEstimation` into the flags (`:39`).
   **The structure is re-optimised from the parent pose, not re-docked from
   scratch.** There is no random shuffle in this mode.

3. **The perturbation** (`:48-60`):
   ```csharp
   if (!_distances[i].IsSelected) continue;           // :54
   rnorm = RandomNorm() * sigmafactor;                // :56, sigmafactor = 1.0
   dtmp.R += (rnorm > 0) ? rnorm * dtmp.ErrPlus
                         : rnorm * dtmp.ErrMinus;     // :57
   ```
   `RandomNorm()` (`SpringEngine.cs:607-641`) is a ratio-of-uniforms standard
   normal, seeded per engine instance (`SimulationBase.cs:129`). Only
   **selected** distances are perturbed.

4. **Resampling.** The job manager calls `SetState(sr)` **before every
   repetition** (`SimulationJobManager.cs:99`; `IsContinuous` is false here,
   unlike `MetropolisSampler.cs:11`), so targets reset to the unperturbed truth
   each time. **The perturbation is resampled per repetition and does not
   compound.**

5. **The ensemble.** `N` x (selected parents) results, tagged
   `ErrorEstimation`, carrying `ParentStructure` (`SpringEngine.cs:298`).

### 1.2 What it is, and what a user reads off it

A **parametric bootstrap** (Monte-Carlo error propagation), *not* a resampling
bootstrap and *not* a posterior sample: the data are replaced by the fitted
model's own distances (residuals zeroed), synthetic datasets are drawn from the
stated per-distance error model, each is re-fitted, and the scatter of the
re-fits is the reported uncertainty.

**The statistic a user reads off is the RMSD of each replica against the
parent.** `SaveForm.SaveEnergyTable` (`SaveForm.cs:359-371`) writes per-result
`RMSD vs <reference>` columns; the grid shows the same (`MainForm.cs:701,719`).
The spread of that column over the replicas is FPS's structural uncertainty. It
is a **precision** figure conditional on the restraint set and the error model,
with no model-misspecification contribution — because the residuals were zeroed.

### 1.3 Is the perturbation symmetric? No.

`:57` branches on the **sign of the standard normal draw**, not on any
density-matching rule. So `R'` is, with probability 1/2 each, a half-normal of
scale `ErrPlus` above `R` or of scale `ErrMinus` below it.

For `ErrPlus = 10`, `ErrMinus = 3`:

| quantity | value |
|---|---|
| P(R' > R) | 0.5 |
| median | exactly `R` |
| mean | `R + (10-3)/sqrt(2*pi)` = `R + 2.79` |
| sd | `sqrt((100+9)/2 - 2.79^2)` = `6.84` |
| density just above vs just below `R` | ratio `3/10` — **discontinuous** |

**Defect D1.** This is a *sign-split* normal, not the standard two-piece (split)
normal, which weights the halves `sigma_+/(sigma_+ + sigma_-)` and
`sigma_-/(sigma_+ + sigma_-)` so the density is continuous at the mode. FPS
gives each half mass 1/2, so the density jumps by `ErrMinus/ErrPlus` at `R` and
the perturbation carries a **systematic positive mean shift of
`(ErrPlus - ErrMinus)/sqrt(2*pi)`** whenever the errors are asymmetric. For
strongly asymmetric FRET errors this biases the ensemble outward along every
such restraint. Reproduce it for bit-comparability; diverge from it knowingly.

**Defect D2.** Step 1 rewrites `R` for *all* distances but step 3 perturbs only
the *selected* ones, and the error-estimation default is `OptimizeSelected = All`
(`ProjectData.cs:113`). Any **deselected** distance therefore becomes a
zero-noise restraint pinning the replica to the parent pose, biasing the
estimated uncertainty **downward**. Invisible in the common case
(`Distances.cs:138` selects everything on load), but it bites as soon as a user
unchecks a restraint.

## 2. Metropolis sampling

**Class:** `MetropolisSampler : SpringEngine`. `IsContinuous = true` (`:11,17`)
— the job manager calls `SetState(initial)` **once** and then `Simulate()`
repeatedly on the same engine (`SimulationJobManager.cs:107-111`), so successive
results are **snapshots along one Markov chain**, not independent runs.

### 2.1 Move set (`MetropolisSampler.Iterate`, `:64-111`)

1. **One rigid body, uniformly**: `i = _rnd.Next(translation.Length)` (`:67`).
   Molecule 0 is included, so the whole assembly can drift; the result is
   re-expressed in molecule 0's frame at save time (`:54-57`).
2. **Translation**: independent Gaussian per Cartesian component,
   **sigma = 0.05 A** (`:74-77`); RMS step `0.05*sqrt(3)` = 0.087 A.
   Hard-coded, *not* an `FPSParameters` field.
3. **Rotation axis**: uniform on the sphere via `z0 ~ U(-1,1)`,
   `phi ~ U(0,2*pi)` (`:84-87`). The comment at `:82` calls it "anisotropic",
   which is wrong — it is isotropic (**D18**).
4. **Rotation angle**: Gaussian, **sigma = 0.0017 rad** (`:80,:90`) = 0.0974 deg.
   Applied on the left, `rotation[i] = rot * rotation[i]` (`:91`).

Both step sizes are compile-time constants. There is no adaptive step-size
tuning and no acceptance-rate targeting.

### 2.2 Acceptance — and what `rkT` is

```csharp
double E = ForceAndTorque();                          // :94
CheckForClashes(true); E += Eclash;                   // :95-96
double recipr_kT = this.SimulationParameters.rkT;     // :99
bool accepted = E < Elast || _rnd.NextDouble() < Math.Exp((Elast - E) * recipr_kT);  // :100
```

Standard Metropolis on the *total* energy (restraints + clash); on rejection the
saved state is restored (`:103-107`).

**`rkT` is a reciprocal temperature, `beta = 1/kT`, despite the name.** The UI
label is "reciprocal kT" (`SimParametersForm.Designer.cs:311`), and the
multiplication `(Elast - E) * rkT` proves it.

**So a larger `rkT` means COLDER sampling.** With the default `rkT = 10`
(`ProjectData.cs:74,90,107,123,143`), an uphill move costing `dE = 0.1`
chi-square units is accepted with probability `exp(-1)` = 0.37; one costing
`dE = 1` with probability `4.5e-5`. **This is a very cold chain relative to the
chi-square scale of the energy** — it samples the immediate basin, not the
restraint-consistent manifold. Do not read `rkT = 10` as "kT = 10".

`rkT` appears in all five modes' defaults but is **read only by
`MetropolisSampler`**; it is dead everywhere else.

### 2.3 Termination

```csharp
do { nsuccessfull += this.Iterate() ? 1 : 0; ntotal++; }
while (nsuccessfull < this.SimulationParameters.MaxIterations);   // :32-37
```

**The loop counts ACCEPTED moves, not attempts.** `MaxIterations` (default 8000
for Sample, `ProjectData.cs:122`) is the number of accepted moves between
snapshots. **D3: there is no cap on `ntotal`** — the commented-out guard at
`:29,:37` shows it was removed. With a cold chain and a rejection-heavy region,
`Simulate()` can run unbounded. A port needs an attempt cap FPS does not have.

### 2.4 What is recorded (`:40-61`)

`E` (restraints + clash, recomputed at `:42-44`), `Eclash`, `Ebond`,
`ParentStructure`, `Translation[]`/`Rotation[]` re-referenced to molecule 0
(`:54-57`), `SimulationMethod |= MetropolisSampling` (`:49`), best-fit
translation (`:60`). Rotations re-orthonormalised at `:41` and again `:56`
(redundant).

**D4.** `sr.Converged = (ntotal < MaxIterations)` (`:48`). The loop exits only
when `nsuccessfull >= MaxIterations`, and `ntotal >= nsuccessfull` always, so
**`Converged` is always `false` for every Metropolis snapshot.**

## 3. Protocol details that change a docked answer

### 3.1 `OptimizeSelected` — and what "selected" means

Enum `Selected | All | SelectedThenAll` (`SimulationBase.cs:5-8`).

**"Selected" refers to *distances*, not molecules.** The gate is
`SpringEngine.cs:428`: `if (optimize_selected_now && !_distances[i].IsSelected) continue;`
`Distance.IsSelected` (`Distances.cs:14`) is the "ON" checkbox in the distances
grid, default `true` (`Distances.cs:138`), persisted as
`ProjectData.SelectedDistances`. `Molecule.Selected` (`Molecule.cs:114`) is a
**different, unrelated flag** — see 3.2.

* **`Selected`** — only checked distances contribute force, torque and energy.
* **`All`** — every distance contributes.
* **`SelectedThenAll`** — two phases, each capped at `MaxIterations / 2`
  (`:255`); `dt` and `viscosity` are re-derived between them (`:271-275`).

**D5.** `sr.Converged = (niter < MaxIterations)` (`:295`) but under
`SelectedThenAll` `niter` is capped at `MaxIterations/2`, so **`Converged` is
unconditionally `true`** in that mode.

Clash energy is *always* computed over all molecules (`Iterate` calls
`CheckForClashes(true)`, `:326`) regardless of `optimize_selected_now` — only
*distance* restraints are gated. A port must not gate clashes.

Screening reduces the enum to a Boolean (`MainForm.cs:1039`): `SelectedThenAll`
and `All` are both treated as "all".

### 3.2 The initial random shuffle — `SetState()` (`SpringEngine.cs:155-198`)

Used **only** by Dock (`SimulationJobManager.cs:83-89`); Refine, Error
estimation and Sample all enter through `SetState(SimulationResult)`.

```csharp
for (i = 0..N-1) { translation[i] = 0; rotation[i] = E; }   // :166-170
do {
  for (int i = 1; i < _molecules.Count; i++) {
    lastshaked = 0;                                          // :175  <-- see D6
    if (!_molecules[i].Selected) {
      tshake = (Vector3(U,U,U) + (-0.5)) * 10.0;             // :178
      translation[i] += tshake;                              // :179
      lastshaked = i;
      rotation[i] = rotshake * rotation[i];                  // :181-187
    } else {
      translation[i] = translation[lastshaked] + CM_i - CM_lastshaked;  // :191
      rotation[i]    = rotation[lastshaked];                 // :192
    }
  }
} while (CheckForClashes(false));                            // :196
return ForceAndTorque();
```

**Amplitude.** `Vector3 + double` is componentwise (`MatrixVector3.cs:26`), so
`tshake` is uniform in a **10 A cube centred on zero — each component in
[-5, +5] A**. Orientation is fully randomised: uniform axis, uniform angle in
`[0, 2*pi)` (`:185`). Hard-coded.

**Which molecules move.** Molecule 0 is never touched (loop starts at `i = 1`).
Molecules with **`Selected == false` are the ones shuffled** (`:176`) — the
*mobile* bodies. Molecules with `Selected == true` are pinned to the frame of
`lastshaked` (`:191-192`).

**D15.** `MainForm.cs:321` comments the handler as "change random-initial-position
property of a molecule", which reads as the opposite of what the code does:
**checking a molecule makes it NOT randomised.** `CheckForClashes`'s own comment
agrees with the code — `// not "all" means only initially mobile molecules`
(`:463`), "mobile" = `!Selected`.

**D6.** `lastshaked = 0` is assigned **inside** the loop (`:175`), so it resets
every iteration and is only ever `0` where it is read (`:191`); the
`lastshaked = i` at `:180` is dead. A `Selected` molecule is therefore *always*
glued to molecule 0's frame, never to the last-shuffled neighbour.

**Does the "no clashes" loop terminate?** Not provably. **D7: no iteration cap,
no timeout, no cancellation check.** In practice yes, and for a reason worth
stating: `translation[i] += tshake` at `:179` **accumulates** across passes (only
the `Selected` branch at `:191` assigns). Each failed pass adds another
independent `U[-5,5]^3` displacement, so the mobile bodies random-walk apart at
RMS `sqrt(n_passes) * 2.89 * sqrt(3)` A and clashes become vanishingly unlikely
after a few passes. The *rotation* is re-randomised rather than accumulated, so
it contributes no escape pressure. Keep the accumulating translation; add a cap.

**D8.** The clash acceptance test (`:466-469`) is
`clash |= CheckForClashes(im1, im2) && (!_molecules[im2].Selected | all);` —
only `im2`'s flag is tested and `im1 < im2` always, so a clash between a
`Selected` molecule at a lower index and a mobile one at a higher index is
counted while the reverse is ignored. Index-order dependent, almost certainly
unintended. (`|` where `||` was meant is harmless here.)

### 3.3 `MaxForce` — the force cap and the linear region

`PrepareSimulation` (`:141-142`): `drmaxplus = MaxForce/k+`,
`drmaxminus = -MaxForce/k-`, i.e. `drmax = MaxForce * Err^2 / 2`.

`ForceAndTorque` (`:432-443`) gives

```
E(dr) = (dr/Err)^2                      for |dr| <= drmax
E(dr) = MaxForce * (|dr| - drmax/2)     for |dr| >  drmax
```

a **Huber-like robustification**, joined C^1-continuously: outliers pull with
bounded force and cost only linearly, so one badly wrong restraint cannot
dominate. **A port using a pure harmonic restraint will not reproduce FPS's
docked poses when any restraint is grossly violated** — exactly the situation
during random-restart docking.

Mode dependence is large: `MaxForce = 400` for Dock/Sample/Screening against
`10000` for Refine/Error estimation. Note the cap engages at
`drmax = MaxForce*Err^2/2`, so it **robustifies the tight restraints and leaves
the loose ones harmonic** — with `Err = 0.1 A` and `MaxForce = 400`,
`drmax = 2 A`; with `Err = 3 A`, `drmax = 1800 A` and the cap is inert.

### 3.4 `ClashTolerance` — how it becomes a spring constant

`:147-149`: `kclash = 2/ClashTolerance^2`, `drmaxclash = -MaxForce/kclash`,
`dr4dt = min(ClashTolerance*0.5, minerror)`.

Same algebra as the distance restraints, so an overlap `dr` costs
`(dr/ClashTolerance)^2`. **`ClashTolerance` is the overlap in angstrom that costs
exactly one chi-square unit per atom pair.** Defaults 1.0 A for
Dock/Sample/Screening, 0.5 A for Refine/Error estimation — refinement penalises
clashes 4x harder per angstrom.

Evaluated by native code (`:540-543` -> `FpsNativeWrapper.CheckForClashes`) after
a two-level bounding-box + cluster-radius cull (`:494-536`).

**D11.** `drmaxclash` is computed at `:148` and **never used**: the native entry
point takes `kclash` only, with no force cap (`FpsNativeWrapper.cs:32-34`). The
capped form exists only in the commented-out managed reference
(`:591-592`). **In the shipping path the clash potential is purely quadratic and
uncapped**, while distance restraints are capped. Do not apply `MaxForce` to the
clash term if you want to match FPS.

**D16.** `InitRoutines` binds `CheckForClashes` only for `Win32`/`Win64`
(`FpsNativeWrapper.cs:218,225`); the Linux and Other branches (`:227-245`) leave
the delegate `null`. Docking cannot run outside Windows in this tree — the AV
routines have managed fallbacks, the clash routine does not.

### 3.5 Adaptive time step and viscosity

A damped position-Verlet on rigid bodies (`SpringEngine.Iterate`, `:318-356`),
each molecule carrying `Mass` and a scalar `SimpleI` = the **largest** principal
moment (`Molecule.cs:284`), used isotropically.

**Time step** (`:224`):
`dt = dr4dt / sqrt(2*max(E, Nmolecules)/minmass) * TimeStepFactor`.
`sqrt(2E/m)` is the velocity the *lightest* body would reach converting all
energy to kinetic; `dt` is the time to travel `dr4dt`. So **`dt` depends on the
current energy (floored at `Nmolecules`), the smallest molecular mass, the
smaller of `ClashTolerance/2` and the smallest error over all distances, and
`TimeStepFactor`.** It shrinks as energy rises — a stability guard, not physics.

**Viscosity** (`:231-234`):
`viscosity_per_dt = ViscosityFactor*2*dt*sqrt(ktot/Mtot)`,
`viscosity = viscosity_per_dt*dt`, clamped to `[1e-4, 0.2]` (`:9-10`), with
`ktot = sum over distances of max(k+, k-)`. `sqrt(ktot/Mtot)` is a system
frequency `omega`, so `viscosity = 2*ViscosityFactor*(omega*dt)^2` — damping
scaled to critical damping of the stiffest collective mode.

**`Cleanup`** (`:368-402`) fires when energy jumps (`Etot > Etot_last + 100`,
`maxdE` at `:11`), when `Eclash > 10`, or periodically (`:262`). It
re-orthonormalises rotations and, if the energy permits, **doubles `dt`** (`:396`).
`dt` is never decreased.

**D9.** The periodic trigger is `if ((niter << 28) == 0 || ...)` (`:262,284`). A
left shift by 28 in 32-bit `int` keeps only the low 4 bits, so the condition is
true **iff `niter % 16 == 0`**. It works; it reads as a no-op. Do not copy the
idiom.

**D10.** In the `SelectedThenAll` phase-2 re-derivation (`:273`) the viscosity is
computed as `ViscosityFactor*2*dt*sqrt(ktot/Mtot)` — **one factor of `dt` short**
of `:231-232` — and `viscosity_per_dt` is not updated, so the first `Cleanup`
recomputes `viscosity = viscosity_per_dt*dt` (`:398`) from the **stale phase-1**
value. Phase 2 briefly runs at the wrong damping then snaps to a stale one.

**D13.** `FPSParameters.ETolerance` (`SimulationBase.cs:20`) is UI-exposed
(`SimParametersForm.cs:41,63`) and stored in every mode's defaults, but is
**never read**. Convergence uses only `FTolerance`, `TTolerance`, `KTolerance`
(`:265-266`).

**D14.** `normF` and `normT` are documented as `Sum of F*F` / `Sum of T*T`
(`:314-315`) and accumulated as squares (`:338,340`), but compared **directly**
against `FTolerance`/`TTolerance` (`:265`). The tolerances are therefore on the
*squared* norms: `FTolerance = 0.001` means `sum|F|^2 < 1e-3`, i.e.
`|F| ~ 0.032`. `K` is compared against `KTolerance * Nmolecules` (`:266`), a
per-molecule kinetic-energy tolerance. A port comparing un-squared norms
converges at a much looser point.

## 4. What FPS calls a "bond"

`PrepareSimulation` (`SpringEngine.cs:111-116`):

```csharp
dist_i.IsBond =
    lp1.Dye == DyeType.Unknown && lp1.AVData.AtomID > 0 && lp1.AVData.AVType == AVSimlationType.None
 && lp2.Dye == DyeType.Unknown && lp2.AVData.AtomID > 0 && lp2.AVData.AVType == AVSimlationType.None;
```

**The test.** A distance is a bond iff **both** endpoints are plain atoms: no dye
(`Dye == Unknown`), a valid PDB atom id (`AtomID > 0`), and **no accessible
volume** (`AVType == None`). That is, an atom-to-atom restraint — a crosslink, an
EPR anchor, a symmetry or geometric constraint — not a dye-to-dye FRET restraint.

**The consequence — vdW radii are mutated** (`:119-134`). Both endpoint atoms get
`AtomData.vdWRNoClash = 0.4 A` (`StaticData.cs:33`; doc at `:26-28`: *"van der
Waals radius which replaces original vdwr for 'ATOM' type LPs"*). The write must
stay in sync across three places: `Molecule.vdWR[natom]` (`:125`),
`Molecule.ClusteredAtomvdwR[...]` (`:129`), and — via `Marshal.Copy` into the
16-byte-aligned float array the native kernel reads — the `w` lane of the packed
`XYZvdwR` vector (`:130-131`). Rationale: a bonded pair is *meant* to be pulled
to a covalent distance, so its atoms must not repel with a full vdW radius.

**D12.** `Molecule` is a reference type (`Molecule.cs:26`) and `SimulationBase`'s
setter copies the *list*, not the molecules (`SimulationBase.cs:45-49`). So the
radius mutation is applied to the **shared, session-global** `Molecule` objects —
it leaks into every clone, every other mode and every subsequent run, and is
**never reverted**. Deselect the bond, change modes, re-dock: the atom keeps its
0.4 A radius until the molecules are reloaded from disk. A port must copy the
radii per run or restore them.

**How it is reported.** `Ebond` is the sum of `dE` over bond distances only
(`:445`), stored on every result (`:294`, `MetropolisSampler.cs:47`;
`SimulationResult.cs:21`) and written as the **`chi2_bond`** column
(`SaveForm.cs:362,365,369`). **`Ebond` is a subset of `E`, not an addition to
it — do not sum them.**

## 5. Refinement — AVs in the docked context

**Yes, it recomputes them there.** `Refinement.RedoAV` (`Refinement.cs:62-87`):

1. `InitializeStructureForRefinement` (`:49-60`) builds `mt = new Molecule(sr)`.
   That constructor (`Molecule.cs:37-81`) **merges every subunit's atoms into one
   flat list, each transformed into the docked pose** (`Molecule.cs:66-68`).
2. `RedoAV` constructs `AVEngine av = new AVEngine(mt, ...)` (`Refinement.cs:67`)
   — on that **merged docked assembly**.
3. It re-runs `Calculate1R`/`Calculate3R` per position (`:77-80`), skipping
   `AVType == None` (`:74`).
4. The resulting `av.Rmp` is transformed **back** into the home subunit's local
   frame (`:81-82`) and stored as the new position (`:83`).

**So a partner subunit CAN occlude a dye, and refinement is precisely the step
that accounts for it.** Initial docking uses AVs on isolated subunits;
refinement replaces them with AVs computed inside the complex. This is a
*first-order self-consistent* step, not iterated to convergence: `Simulate`
(`:32-40`) does exactly one minimise -> re-AV -> re-score cycle per call.

`SetState` (`:23-30`) recomputes AVs on entry **unless** the input already
carries `Refinement` *and* non-null `RefinedLabelingPositions` (`:26-28`).

Downstream, `ModelDistance` (`SimulationResult.cs:172-176`) prefers
`RefinedLabelingPositions` when the flag is set — **so error estimation launched
on a refined parent takes its "truth" from the occlusion-corrected positions**
(§1.1 step 1).

Cost note: `AVGlobalParameters.ESamples = 200000` (`ProjectData.cs:162`), and
`RedoAV` runs a full AV grid search per position per refinement call.

## 6. The five modes and their default parameters

From `ProjectData()` (`ProjectData.cs:61-149`) — the values shipped in the tool,
and therefore the settings behind the published work unless a project overrides
them.

| Field | Dock | Refine | Error estimation | Sample | Screening |
|---|---|---|---|---|---|
| `ViscosityFactor` | 1.0 | **0.7** | **0.7** | 1.0 | 1.0 |
| `TimeStepFactor` | 1.0 | **0.5** | 1.0 | 1.0 | 1.0 |
| `MaxIterations` | 200 000 | **500 000** | 100 000 | **8 000** | 200 000 |
| `MaxForce` | 400.0 | **10 000.0** | **10 000.0** | 400.0 | 400.0 |
| `ClashTolerance` (A) | 1.0 | **0.5** | **0.5** | 1.0 | 1.0 |
| `rkT` (= 1/kT) | 10 | 10 | 10 | 10 (**only mode that reads it**) | 10 |
| `ETolerance` | 100 | 100 | 100 | 100 | 100 (**dead in all modes**) |
| `KTolerance` | 0.001 | **0.0005** | 0.001 | 0.001 | 0.001 |
| `FTolerance` | 0.001 | **0.0005** | 0.001 | 0.001 | 0.001 |
| `TTolerance` | 0.02 | **0.01** | 0.02 | 0.02 | 0.02 |
| `OptimizeSelected` | **Selected** | **All** | **All** | **Selected** | **Selected** |
| source lines | `:67-81` | `:84-98` | `:101-115` | `:118-132` | `:135-149` |

* `MaxIterations` is entered in the UI in **thousands** (`SimParametersForm.cs:37,59`;
  label "Max iterations (k)"). The table gives the stored integer.
* For **Sample**, `MaxIterations` counts **accepted** Metropolis moves per
  snapshot (§2.3) — 8 000 is not comparable to 200 000 *integrator steps*.
* `MaxForce` is 25x larger for Refine/Error estimation than for Dock — the
  outlier-robust linear region is effectively **off** for refinement and error
  estimation and **on** for docking (§3.3). The single largest behavioural
  difference between the modes.
* Screening uses a different engine (`FilterEngine`/`FilterJobManager`) and reads
  only `OptimizeSelected` and the conversion `R0` (`MainForm.cs:1036-1041`).

Non-`FPSParameters` defaults that matter:

| Field | Default | Line |
|---|---|---|
| `ConversionParameters.R0` | 52.0 A | `:152` |
| `ConversionParameters.PolynomOrder` | 3 | `:153` |
| `AVGlobalParameters.GridSize` (relative) | 0.2 | `:158` |
| `AVGlobalParameters.MinGridSize` (A) | 0.4 | `:159` |
| `AVGlobalParameters.LinkerInitialSphere` | 0.5 | `:160` |
| `AVGlobalParameters.LinkSearchNodes` | 3 | `:161` |
| `AVGlobalParameters.ESamples` | 200 000 | `:162` |
| repetitions spinner | 10 (max 100 000) | `MainForm.Designer.cs:870-888` |

## 7. What "repetitions" means, per mode

`SimulationJobManager.DoJob` (`:33-122`) with `FinalStatesPerInitialState` set in
`MainForm.cs:1054-1082`:

| Mode | `InitialStates` | `FinalStatesPerInitialState` | One "repetition" is |
|---|---|---|---|
| **Dock** | `null` (`MainForm.cs:1060`) | UI spinner (`:1062`) | One **independent random restart**: `SetState()` shuffle + minimisation (`SimulationJobManager.cs:83-89`). |
| **Refine** | selected result rows (`:1049-1052`) | **hard-coded 1** (`:1068`) | One deterministic refinement per parent. **The repetitions spinner is ignored.** |
| **Error estimation** | selected result rows | UI spinner (`:1073`) | One **fresh perturbation draw** re-optimised from the parent pose; i.i.d. across repetitions. |
| **Sample** | selected result rows | UI spinner (`:1078`) | One **snapshot** along a single continuous chain per parent (`IsContinuous`); consecutive snapshots are **correlated**. |
| **Screening** | n/a | n/a | `FilterJobManager`; one job per input structure (`:1082-1088`). |

Total planned work is `FinalStatesPerInitialState * max(1, InitialStates.Length)`
(`SimulationJobManager.cs:36`). Parallelism is capped at `NThreads` and, for
continuous jobs, additionally at `InitialStates.Length` so each chain lives on
one worker (`:42-43`). Each worker holds its **own clone** (`:50`) with its own
RNG (`SimulationBase.cs:129`) and its own distance/position lists — but
**shared `Molecule` objects** (D12).

## 8. Port checklist — what must match to reproduce FPS numbers

Ordered by how much a mismatch moves the answer.

1. **The chi-square restraint with asymmetric `k = 2/Err^2` and the `MaxForce`
   linear tail** (§0, §3.3). Pure harmonic will not match under docking's
   `MaxForce = 400`.
2. **The clash term `kclash = 2/ClashTolerance^2`, uncapped** (§3.4), with bonded
   atoms at `vdWR = 0.4 A` (§4).
3. **`OptimizeSelected` gating on *distance* `IsSelected`; clashes always
   global** (§3.1).
4. **The random-restart placement**: uniform orientation, +/-5 A per axis,
   accumulate-until-clash-free, molecule 0 fixed, `Selected` molecules glued to
   molecule 0 (§3.2).
5. **Error estimation's sign-split perturbation and its zeroed-residual "truth"**
   (§1), including the mean shift for asymmetric errors.
6. **`rkT` as `1/kT` = 10, and `MaxIterations` counting *accepted* moves**
   (§2.2, §2.3).
7. **Refinement's docked-context AV recomputation** (§5), and that a refined
   parent feeds occlusion-corrected distances into error estimation.
8. Convergence tolerances are on **squared** force/torque norms (§3.5).

## 9. Defect index

| # | Where | What |
|---|---|---|
| D1 | `ErrorEstimation.cs:57` | Sign-split, not density-matched, split normal: density discontinuous at `R`, mean shifted by `(Err+ - Err-)/sqrt(2*pi)`. |
| D2 | `ErrorEstimation.cs:23-29` + `ProjectData.cs:113` | Deselected distances become zero-noise pins to the parent pose under the default `OptimizeSelected = All`, biasing the estimated error downward. |
| D3 | `MetropolisSampler.cs:32-37` | No cap on attempted moves; only accepted moves are counted. Unbounded loop. |
| D4 | `MetropolisSampler.cs:48` | `Converged` is always `false` for sampling snapshots. |
| D5 | `SpringEngine.cs:295` + `:255` | `Converged` is always `true` under `SelectedThenAll`. |
| D6 | `SpringEngine.cs:175/180/191` | `lastshaked` reset inside the loop; `lastshaked = i` is dead. `Selected` molecules are always glued to molecule 0. |
| D7 | `SpringEngine.cs:171-196` | No cap on the "until no clashes" retry loop; terminates only because `tshake` accumulates. |
| D8 | `SpringEngine.cs:468` | Clash-acceptance test reads only `_molecules[im2].Selected`; index-order dependent. Uses `|` where `||` was meant. |
| D9 | `SpringEngine.cs:262,284` | `(niter << 28) == 0` is an obfuscated `niter % 16 == 0`. |
| D10 | `SpringEngine.cs:273` vs `:231-232,398` | `SelectedThenAll` phase 2 computes viscosity with one factor of `dt` missing, then reverts to a stale `viscosity_per_dt`. |
| D11 | `SpringEngine.cs:148` + `FpsNativeWrapper.cs:32-34` | `drmaxclash` computed but unused; the shipping clash potential is uncapped while distance restraints are capped. |
| D12 | `SpringEngine.cs:119-134` + `SimulationBase.cs:45-49` | Bond-driven vdW mutation applied to shared `Molecule` objects; leaks across modes/clones/runs and is never reverted. |
| D13 | `SimulationBase.cs:20` | `ETolerance` is UI-exposed and stored but never read. |
| D14 | `SpringEngine.cs:265` vs `:314-315` | Force/torque tolerances applied to **squared** norms; the UI labels imply otherwise. |
| D15 | `MainForm.cs:321-325` vs `SpringEngine.cs:176` | The comment says the checkbox sets "random-initial-position"; the code makes a checked molecule the one **not** randomised. |
| D16 | `FpsNativeWrapper.cs:227-245` | `CheckForClashes` bound only on Windows; docking cannot run elsewhere. |
| D17 | `SimulationResult.cs:40-71` | `RMSD(..., selected_only: true)` returns `sqrt(0/0)` = `NaN` when no molecule is `Selected` — the common case. That is the `"RMSD vs ref (selected)"` column. |
| D18 | `MetropolisSampler.cs:82` | Comment calls the axis construction "anisotropic"; it is isotropic. |
