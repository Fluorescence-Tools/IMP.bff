# The probability calculus of a Bayesian fluorescence decay analysis

*Written 2026-09-08 for the `bayesian_decays` notebook series. It is the
reference the notebooks point at: every symbol they use is defined here, and
every approximation they make is named here together with what it costs.*

---

## 0. What is being asked

A donor and an acceptor dye sit on a molecule at a distance $R$. Energy
transfer shortens the donor's fluorescence lifetime by a factor that depends
on $R$, so the shape of the donor's decay carries the distance. A sample is
not one molecule but an ensemble, and the ensemble is not one distance but a
**distribution** $p(R)$. The question is not "what is $R$" but

$$
p\big(p(R)\,\big|\,\text{the photons}\big),
$$

a posterior over a function, together with honest uncertainty on anything one
computes from it — its mean, its width, the weight it puts beyond some
distance.

Everything below is the machinery that turns a set of photon counts into that
posterior **without a sampler**: the model, the priors, the marginalisation
of the one hyperparameter that cannot be maximised, the Gaussian
approximation at each of its grid nodes, and the mixture of those Gaussians.

Distances are reported in units of the Förster radius, $r = R/R_0$, because
that is the coordinate the data actually constrain.

---

## 1. The data and their likelihood

### 1.1 What a TCSPC histogram is

A pulsed laser fires every $T$ nanoseconds. A detector records, for each
photon, the delay $t$ between the last pulse and the detection. Collecting
many photons gives a histogram over $t \in [0, T)$ — the **decay**. Every
photon is an independent event, so with a fixed acquisition time the count in
bin $b$ of histogram $k$ is Poisson:

$$
y_{kb} \sim \mathrm{Poisson}(\lambda_{kb}),
\qquad
\log L(y \mid \theta) = \sum_{k}\sum_{b \in M_k}
\Big[ y_{kb}\log \lambda_{kb}(\theta) - \lambda_{kb}(\theta) - \log y_{kb}! \Big].
$$

$M_k$ is a mask: the bins of histogram $k$ that enter the likelihood. It
exists because under interleaved excitation one histogram holds two pulses,
and one sometimes wants only one of them.

**The Poisson likelihood is not a detail.** Weighting residuals by
$\sqrt{y}$ instead of $\sqrt{\lambda}$ — Neyman instead of Pearson — biases
every amplitude downwards, because a bin that fluctuated low is given more
weight than one that fluctuated high. Counts of a few per bin are common in
the tail of a decay, which is exactly where the long-distance information
lives.

### 1.2 The time axis

The fine axis has $n$ channels of width $\delta t$ covering the period
$T = n\,\delta t$. Decays change fastest immediately after a pulse and are
nearly flat before the next, so the fine axis is **rebinned** by a matrix
$\mathbf{R} \in \{0,1\}^{n_\text{bin} \times n}$ whose rows are contiguous
blocks: narrow right after each pulse, wide before the next. Rebinning is a
sum, so it commutes with everything linear that comes before it, and the
model may be built on either axis.

---

## 2. The forward model: from parameters to expected counts

### 2.1 The chain

$$
\theta
\;\xrightarrow{\ \text{physics}\ }\;
\alpha_k \in \mathbb{R}^{K}
\;\xrightarrow{\ \text{instrument}\ }\;
a_k = e^{s_{\mathrm{samp}(k)}}\alpha_k + c^{\mathrm{sc}}_k\,\|\cdot\|\,e_0 + c^{\mathrm{bg}}_k\,\|\cdot\|\,e_{K-1}
\;\xrightarrow{\ \mathbf{B}_{d(k)}\ }\;
\lambda_k .
$$

$\alpha_k$ are **abstract amplitudes**: the weight of each of $K$ basis
columns. $s_{\mathrm{samp}}$ is one log scale per sample (molecules times
acquisition time). $c^{\mathrm{sc}}$ and $c^{\mathrm{bg}}$ are the scatter and
background fractions of that channel's total, loaded onto the first and last
columns.

### 2.2 The basis

For detector $d$ with instrument response $\mathrm{IRF}_d$, the columns of
$\mathbf{B}_d$ are

* column $0$: $\mathrm{IRF}_d$ itself — scattered laser light has no lifetime;
* columns $1 \dots K-2$: the **periodic** steady-state response to each grid
  lifetime $\tau_i$,
  $$
  \big(\mathrm{IRF}_d * e^{-t/\tau_i}\big)_{\text{periodic}}(t)
  = \sum_{m \ge 0} \big(\mathrm{IRF}_d * e^{-\cdot/\tau_i}\big)(t + mT),
  $$
  because the sample is still glowing from the previous pulse when the next
  arrives;
* column $K-1$: flat — uncorrelated background.

Every column is normalised to unit sum, so total counts are the sum of the
amplitudes and the amplitudes are counts.

The response is either analytic — a Gaussian of width $w_d$ centred at
$0.1T + \sigma_d$, skewed by an error-function factor for the detector's
transit-time tail — or measured. When the "response" is really a reference
dye's decay of lifetime $\tau_r$, the exact correction is the
delta-function convolution identity

$$
\mathrm{IRF} * e^{-\cdot/\tau_i}
= \mathcal{R} + \Big(\tfrac{1}{\tau_r} - \tfrac{1}{\tau_i}\Big)\,
  \big(\mathcal{R} * e^{-\cdot/\tau_i}\big),
\qquad \mathcal{R} = \mathrm{IRF} * e^{-\cdot/\tau_r},
$$

which needs no deconvolution (Zuker, Szabo, Bramall, Krajcarski & Selinger,
*Rev. Sci. Instrum.* **56**, 14, 1985).

### 2.3 The distance enters through a quenching rate

Under the **homogeneous approximation**, every component of a
multi-exponential donor is quenched by the same rate,

$$
k_{\mathrm{FRET}}(r) = \frac{1}{\tau_0}\, r^{-6},
\qquad r = R/R_0,
$$

with $\tau_0$ a single reference lifetime. Then the FRET factor separates
from the donor's own decay law: a donor with lifetime spectrum
$\{(x_i, \tau_i)\}$ at distance $r$ decays as
$\sum_i x_i \exp\!\big[-t(1/\tau_i + k_{\mathrm{FRET}}(r))\big]$, and the
whole distance dependence is one number per grid point (Peulen, Opanasyuk &
Seidel, *J. Phys. Chem. B* **121**, 8211, 2017). The alternative —
quenching each component at its own rate — is a different physical claim and
a different map; it is a switch in the code, not an approximation made
silently.

The transfer efficiency is

$$
E(r) = \frac{1}{1 + r^{6}} .
$$

Its derivative $|\mathrm{d}E/\mathrm{d}r| = 6r^5/(1+r^6)^2$ vanishes at both
ends, which is the whole reason a distance distribution has a **window** and
not a range: see §7.3.

### 2.4 Species, and why the amplitudes are linear in $p$

Write $u_j = (1-x_{D0})\,p_j$ for the number of doubly labelled molecules at
grid distance $r_j$, and $x_{D0}$ for the fraction carrying no acceptor. Then
for the labelled sample under the green pulse

$$
\alpha^{\mathrm{don}}
= Q_D \Big( x_{D0}\,\mathbf{base} + \sum_j u_j\, \mathbf{S}^{R}_j \Big) \mathbf{c}_D,
\qquad
\alpha^{\mathrm{acc}}
= Q_A \Big( \textstyle\sum_{j,l} u_j w^A_l\, \mathbf{S}^{Ag}_{jl} \mathbf{c}_D
  + \varepsilon_{AG} \big(\textstyle\sum_j u_j\big)\, \mathbf{a}^{\mathrm{dir}} \Big),
$$

with $\mathbf{c}_D$ the donor's lifetime spectrum, $w^A$ the acceptor's
lifetime weights, $\varepsilon_{AG}$ the probability that the green laser
excites the acceptor directly, and $\mathbf{S}^{R}, \mathbf{S}^{Ag}$
precomputed maps from (distance, lifetime) to basis columns.

**Both linearities matter.** $\alpha$ is linear in $u$ given everything else,
and linear in $\mathbf{c}_D$ given $u$. The first makes the Jacobian with
respect to the distance distribution a matrix rather than one forward pass
per coordinate; the second makes the starting values a non-negative least
squares problem. And putting all the weight on one grid point,
$p = e_j$, gives **one molecule's decay** — which is what lets the same
forward model simulate single molecules and fit ensembles.

### 2.5 Polarisation

Excitation by polarised light leaves the emission polarised until the dye
rotates. With fundamental anisotropy $r_0$ and rotational correlation time
$\rho$, the anisotropy decays as $r(t) = r_0 e^{-t/\rho}$ and

$$
I_\parallel \propto \tfrac{1}{3}\,I(t)\,[1 + 2r(t)],
\qquad
I_\perp \propto \tfrac{1}{3}\,I(t)\,[1 - r(t)] .
$$

Real detectors are imperfect: a factor $g$ for the relative sensitivity of the
perpendicular channel, and mixing coefficients $\ell_1, \ell_2$ for
depolarisation in the optics, so the model uses
$2 - 3\ell_1$ and $-1 + 3\ell_2$ in place of $2$ and $-1$. The rotational
time is not one number but a distribution on a grid, $w^\rho$, because a dye
on a linker has more than one motion.

### 2.6 Colour mixing

Emission does not stay in its own channel:

$$
\begin{pmatrix}\text{green detector}\\ \text{red detector}\end{pmatrix}
=
\begin{pmatrix} G_g C_{gD} & G_g C_{gA} \\ G_r C_{rD} & G_r C_{rA}\end{pmatrix}
\begin{pmatrix}\text{donor light}\\ \text{acceptor light}\end{pmatrix},
$$

with $G$ the detection efficiencies and $C$ the spectral crosstalks. All six
are parameters with priors, not constants.

### 2.7 Interleaved excitation: six scopes

With two lasers per period, "which sample" and "which pulse" are different
questions, and the pair decides the physics:

| scope | sample | pulse | what it holds |
|---|---|---|---|
| D0 | donor only | green | the unquenched donor |
| DR | donor only | red | the donor excited by the red laser ($\varepsilon_{DR}$, tiny) |
| DA | labelled | green | donor, sensitised acceptor, acceptor by direct excitation |
| DAr | labelled | red | the acceptors of the $1-x_{D0}$ molecules that have one |
| AG | acceptor only | green | the acceptor excited by the green laser |
| A0 | acceptor only | red | the acceptor under its own laser |

The two windows of one sample share that sample's concentration and
acquisition time, so **their ratio measures an excitation crosstalk**. That is
what makes $\varepsilon_{AG}$ identifiable, and it is why measuring the
reference samples under both lasers removes the donor-only fraction's bias
almost entirely (§7.2).

A histogram that holds both pulses is modelled as the sum of the two scopes'
expected counts, both periodic, so the cross-window tails are in the data and
in the model.

---

## 3. The unknowns

$\theta$ collects, with the prior each carries:

| group | symbol | prior |
|---|---|---|
| distance distribution | $c \in \mathbb{R}^{n_c-1}$, sum-to-zero | P-spline, §4 |
| roughness weight | $\log_{10}\lambda$ | uniform on $[-2, 5]$ |
| donor-only fraction | $x_{D0}$ | uniform on $[0, 0.5]$ |
| donor lifetime spectrum | $\epsilon \in \mathbb{R}^{K_\text{int}}$, softmax | $\mathcal{N}(0, 3^2)$ per point, plus continuity |
| acceptor lifetime weights | $w^A$ | logistic normal, free-form |
| rotational times | $w^\rho, w^{\rho_A}$ | logistic normal around a bump |
| anisotropies | $r_0^D, r_0^A$ | uniform on $[0.15,0.45]$, $[0,0.45]$ |
| calibration | $g,\ \ell_1,\ \ell_2,\ C_{gD},C_{gA},C_{rD},C_{rA},\ G_g,G_r,\ Q_D,Q_A,\ \varepsilon_{AG}$ | log normal or normal, widths from what a calibration measurement achieves |
| instrument | $\sigma_d$ (shift), $w_d$ (width), skew | normal / log normal per detector |
| per channel | $s_{\mathrm{samp}}$, $c^{\mathrm{sc}}_k$, $c^{\mathrm{bg}}_k$ | normal on the log; half a decade for the fractions |

Every parameter is optimised in an **unconstrained coordinate** $z$ with a
transform $\theta = T(z)$ — log, logit, additive log-ratio for simplex
variables, sum-to-zero for the spline gauge — and the prior carries the
Jacobian $\log|\partial T/\partial z|$. This is what lets an unconstrained
Newton method work on a constrained problem.

A note on one of these that was a bug: the sum-to-zero basis must be a
**fixed** orthonormal complement of the constant vector. Taking it from an
SVD of a degenerate subspace makes it machine-dependent, and a posterior
stored on one machine restores wrongly on another. Helmert contrasts are
canonical and platform-independent.

---

## 4. The distance distribution and its prior

### 4.1 The parameterisation

On the grid $r_1 \dots r_{n_R}$ with cubic B-spline basis
$\Phi \in \mathbb{R}^{n_R \times n_c}$,

$$
p(r_j) = \operatorname{softmax}(\Phi c)_j = \frac{e^{(\Phi c)_j}}{\sum_{j'} e^{(\Phi c)_{j'}}} .
$$

Positivity and normalisation hold by construction, and the softmax is shift
invariant, which is the gauge freedom removed by the sum-to-zero constraint on
$c$.

### 4.2 The roughness prior

A P-spline (Eilers & Marx, *Statist. Sci.* **11**, 89, 1996): a difference
penalty on the coefficients rather than a knot-placement decision,

$$
\pi(c \mid \lambda)
\;\propto\;
\lambda^{\,\mathrm{rank}(\mathbf{D}_m\Phi)/2}
\exp\!\Big(-\tfrac{\lambda}{2}\,\|\mathbf{D}_m c\|^2\Big),
$$

with $\mathbf{D}_m$ the $m$-th difference operator. It is **improper**: the
null space of $\mathbf{D}_m$ is unpenalised, which is why the normaliser
carries the *rank* and not the full dimension, and why a weak proper prior is
needed on the null directions.

**The order is a modelling statement about the null hypothesis.** With
$m = 2$ the unpenalised functions are linear in $\log p$, so
$\lambda \to \infty$ gives $p(r) \propto e^{\beta r}$ — an exponential, not a
flat distribution. That is where the shoulder at the edge of the grid in an
empty-sample fit comes from: it is the prior's null hypothesis showing
through, not a numerical artefact. With $m = 3$ the null space is quadratic in
$\log p$, that is, **one Gaussian population is free at every $\lambda$**.
Order 3 gives narrower bands and the same bias, and it was tried and not made
the default because the coverage rule declined it.

Penalising $\log p$ rather than $p$ or $\sqrt{p}$ is also a choice: in log
space the penalty is scale free, and a narrow population costs the same
wherever it sits. The alternatives were tried and withdrawn.

---

## 5. The posterior, and the one thing that must not be maximised

### 5.1 The joint

$$
\pi(\theta, \lambda \mid y)
\;\propto\;
L(y \mid \theta)\;\pi(c \mid \lambda)\;\pi(\theta_{\setminus c})\;\pi(\lambda).
$$

### 5.2 Why the joint mode is wrong

Maximising over $(c, \lambda)$ together fails, and not marginally. The
normaliser $\lambda^{r/2}$ *rewards* a large $\lambda$ whenever $c$ is smooth,
so the joint maximum runs away to a distribution with no populations at all.
Measured on a three-population truth: deviance per degree of freedom 1.16
against a floor of 0.99, and the three populations flattened into one broad
hump. This is the classical failure of joint MAP over a variance component,
and it is why the hyperparameter is treated as a hyperparameter.

### 5.3 The marginal, on a grid

Integrate $\lambda$ out. For each node $\lambda_\ell$ of a grid, do a Laplace
approximation of everything else:

$$
\hat\theta_\ell = \arg\max_\theta \; \log \pi(\theta, \lambda_\ell \mid y),
\qquad
\mathbf{H}_\ell = -\nabla^2 \log \pi(\theta,\lambda_\ell\mid y)\big|_{\hat\theta_\ell},
$$

$$
\log Z(\lambda_\ell)
\;\approx\;
\log \pi(\hat\theta_\ell, \lambda_\ell \mid y) + \frac{d}{2}\log 2\pi - \frac{1}{2}\log|\mathbf{H}_\ell| ,
$$

and the posterior of everything else is the **evidence-weighted mixture of
those Gaussians**

$$
\pi(\theta \mid y) \;\approx\; \sum_\ell w_\ell\, \mathcal{N}\!\big(\theta;\, \hat\theta_\ell,\, \mathbf{H}_\ell^{-1}\big),
\qquad
w_\ell = \frac{Z(\lambda_\ell)\,\pi(\lambda_\ell)}{\sum_{\ell'} Z(\lambda_{\ell'})\,\pi(\lambda_{\ell'})} .
$$

This is INLA's treatment of a hyperparameter (Rue, Martino & Chopin,
*J. R. Statist. Soc. B* **71**, 319, 2009), and marginal-likelihood choice of
a smoothing parameter (Wood, *JRSS-B* **73**, 3, 2011) with the choice
replaced by an average. The grid is refined around the evidence maximum, as
INLA explores around the mode of the hyperparameter.

There is no sampler anywhere: the weights are evidences and the draws, when
draws are wanted, come from Gaussians.

### 5.4 Finding each mode

Fisher scoring with Levenberg–Marquardt damping. The expected information of
Poisson counts with mean $\lambda = B a$ is

$$
\mathcal{I}(\theta) = \mathbf{J}^{\!\top} \operatorname{diag}(1/\lambda)\, \mathbf{J},
\qquad
\mathbf{J} = \frac{\partial \lambda}{\partial \theta},
$$

positive semi-definite by construction, which an observed Hessian is not.
Scoring converges only linearly where expected and observed information
differ, so an exact Newton step is taken periodically, with the curvature
reflected where it is negative.

The reported covariance is $\mathbf{H}^{-1}$ with $\mathbf{H}$ the expected
information plus the prior curvature.

---

## 6. Reporting: from a Gaussian in $c$ to a statement about $p(R)$

### 6.1 Why not just draw

Draws from $\mathcal{N}(\hat\theta, \Sigma)$ pushed through the softmax
overshoot: the Gaussian is a local approximation and the softmax is a wall,
so a draw that is unremarkable in $c$ can be a spike in $p$. Summaries are
therefore computed by the **delta method** rather than by drawing.

### 6.2 Summaries

For a functional $T(\theta)$ — the mean $\sum_j r_j p_j$, the standard
deviation, the mass beyond a distance, the donor-only fraction —

$$
\mathbb{E}[T] \approx T(\hat\theta),
\qquad
\operatorname{Var}[T] \approx \nabla T(\hat\theta)^{\!\top}\, \Sigma\, \nabla T(\hat\theta),
$$

per node, and then over the mixture

$$
\mu = \sum_\ell w_\ell \mu_\ell,
\qquad
\sigma^2 = \sum_\ell w_\ell\big[\sigma_\ell^2 + (\mu_\ell - \mu)^2\big].
$$

The band on $p(r_j)$ itself is the same construction at every grid point, in
$p$ and **clipped at zero** — not in $\log p$. Where $p$ is empty the Laplace
standard deviation of $\log p$ is enormous (a flat direction), and
$\exp(m + 2\sigma)$ reaches absurd values, which makes the figure unreadable
and the claim false.

Intervals from this construction are Gaussian and can cross a boundary a
quantity cannot: a negative lower edge on a mass means the data do not exclude
zero, not that negative mass is possible.

### 6.3 Ranks

For calibration one needs $\Pr(T < T_{\text{true}} \mid y)$ under the mixture,
which is available in closed form as
$\sum_\ell w_\ell\, \Phi\!\big((T_{\text{true}} - \mu_\ell)/\sigma_\ell\big)$.

---

## 7. Whether to believe it

### 7.1 Rule 0: the fit statistic and the residuals

For Poisson counts the deviance is twice the gap to the saturated model:

$$
D = 2\sum_{k,b}\Big[ y_{kb}\log\frac{y_{kb}}{\lambda_{kb}} - (y_{kb} - \lambda_{kb})\Big],
\qquad y\log y \equiv 0 \text{ at } y=0 .
$$

Two things must be said about $D$ and usually are not.

**Its reference is not the number of degrees of freedom.** For counts of the
size a decay actually has, $\mathbb{E}[D]$ of a *correct* model is close to
the number of bins, while $\mathrm{dof} = n_{\text{bins}} - \dim\theta$
subtracts the full parameter dimension even though a penalised spline uses far
fewer of them. The honest reference is measured: draw Poisson realisations at
the fitted means and compute their deviance against those means. A fit is good
when $D/\mathrm{dof}$ sits inside that scatter.

**A statistic without residuals is not a check.** The weighted (Pearson)
residual

$$
w_{kb} = \frac{y_{kb} - \lambda_{kb}}{\sqrt{\lambda_{kb}}}
$$

is plotted for every histogram and its signs are put through a runs test: a
model with the right total misfit and the wrong shape has structured
residuals, and only the runs test sees it.

### 7.2 Simulation-based calibration

A posterior that fits is not a posterior that is calibrated. Draw a truth from
the prior, simulate, fit, and record the rank of the truth in the posterior.
Over many truths those ranks must be **uniform** (Talts, Betancourt, Simpson,
Vehtari & Gelman, arXiv:1804.06788, 2018). Non-uniformity says which way the
posterior is wrong: a peak in the middle means over-wide intervals, peaks at
the ends mean over-narrow.

The measured state of this model, in short: the calibration constants pass;
the distribution's width and its far-tail mass come out about one standard
deviation high and the donor-only fraction low, and those three move together
because information the nuisances do not pin down is absorbed by $p(R)$ as
width. Measuring the reference samples under **both** lasers removes the
donor-only fraction's bias almost entirely and two thirds of the width bias,
lifting the width's interval coverage from 0.52 to 0.83.

### 7.3 The window, and the resolution

Since $E(r) = 1/(1+r^6)$ flattens at both ends, the distance is unbounded
wherever $E$ is within $\delta$ of $0$ or $1$:

$$
r_{\min} = \Big(\frac{\delta}{1-\delta}\Big)^{1/6},
\qquad
r_{\max} = \Big(\frac{1-\delta}{\delta}\Big)^{1/6}
\;\;\Rightarrow\;\;
r \in [0.46,\ 2.15] \ \text{ at } \delta = 0.01 .
$$

Outside that window a fitted curve is the prior, and every figure greys it.

Inside it, resolution was measured rather than asserted, by refitting
noiseless two-population truths: two populations of width $0.06\,R/R_0$
separate reliably at a centre separation of $0.30$, sometimes at $0.20$, never
at $0.15$ — and never at all above about $1.4\,R/R_0$, where the truth is no
longer even representable on the grid.

---

## 8. Single molecules

### 8.1 Bursts

In a confocal single-molecule experiment the photons arrive in bursts, one
molecule at a time. Each photon carries a macro time (which pulse) and a micro
time (delay after it). A burst search on the macro times selects the bright
stretches; the micro times of the selected photons make the decays.

The selection is not free of consequences. It keeps bright crossings, and
brightness depends on $r$, so the distribution behind the selected bursts is
not the one in the cuvette. The decays measure the former.

### 8.2 Classify, then pool

Each burst carries its own low-precision estimate of the same distance:

$$
E^* = \frac{n_{\text{red}|\text{green pulse}}}{n_{\text{green}|\text{green pulse}} + n_{\text{red}|\text{green pulse}}},
\qquad
S = \frac{n_{|\text{green pulse}}}{n_{\text{total}}},
\qquad
\langle t \rangle_D = \frac{1}{n}\sum_{\text{donor photons}} t_i .
$$

$E^*$ is intensity, $\langle t\rangle_D$ is the clock, and they are
independent measurements of the same quantity; $S$ separates molecules whose
acceptor does not work at all. With $N$ photons the per-burst precision is
$\sigma_{E^*} \approx \sqrt{E(1-E)/N}$ and
$\sigma_{\langle t\rangle} \approx \tau/\sqrt{N}$, so a burst of a few hundred
photons localises $E$ to a few percent — enough to *group* molecules even
though it is far too coarse to place one.

**Per-burst trust, and a quantity not to use for it.** A clustering can report
a membership strength per point, and it is tempting to weight or threshold
bursts by it. It is a rank *within* a cluster: the denominator is that
cluster's own death density, so every cluster attains 1 somewhere however
diffuse it is, and equal strengths in two different clusters mean different
things. The comparable quantity is the density at which a point leaves its
cluster — the condensed tree's $\lambda = 1/d$ on the row where that point is
the child — which is in the units of the standardised observables and can be
compared between a tight group and a diffuse one.

Grouping the bursts and pooling each group factorises the problem:

$$
\pi\big(\{p_g\}, \theta_{\text{nuis}} \mid y\big)
\;\propto\;
\prod_g L\big(y_g \mid p_g, \theta_{\text{nuis}}\big)\;
\prod_g \pi(p_g \mid \lambda_g)\;\pi(\theta_{\text{nuis}}),
$$

with the nuisances shared. Each $p_g$ is a much narrower distribution than the
pooled $p = \sum_g w_g p_g$, and a narrow distribution is what the roughness
prior and the resolution limit of §7.3 handle best. **The gain is real but the
step is an approximation**: the assignment is hard, so a misassigned burst
puts its photons in the wrong group, and the groups' distributions are
correspondingly contaminated. The honest version of this is a mixture model
with soft responsibilities; hard assignment is its zero-temperature limit and
is what is done here.

---

## 9. Every approximation in one place

| approximation | where | what it costs |
|---|---|---|
| homogeneous quenching | §2.3 | a different physical model if the donor's components are quenched differently |
| Laplace at each node | §5.3 | wrong in proportion to the posterior's skew; the softmax wall is the worst offender |
| finite grid of $\lambda$ | §5.3 | negligible once refined around the maximum |
| expected information for the curvature | §5.4 | equals the observed Hessian at the mode; differs on the way there |
| delta method for summaries | §6.2 | first order; a second-order correction exists and moves little |
| Gaussian intervals for bounded quantities | §6.2 | intervals can cross a boundary |
| grid discretisation of $p(R)$ | §4.1 | sets the resolution floor: nothing narrower than a grid cell |
| hard burst assignment | §8.2 | contaminates groups by the misassignment rate |

---

## 10. References

* Campello, Moulavi & Sander, *Density-Based Clustering Based on Hierarchical Density Estimates*, PAKDD 2013.
* Eilers & Marx, *Flexible smoothing with B-splines and penalties*, Statist. Sci. **11**, 89, 1996.
* Kennedy & O'Hagan, *Bayesian calibration of computer models*, JRSS-B **63**, 425, 2001.
* McInnes, Healy & Astels, *hdbscan*, J. Open Source Software **2**, 205, 2017.
* Peulen, Opanasyuk & Seidel, *Combining graphical and analytical methods with molecular simulations*, J. Phys. Chem. B **121**, 8211, 2017.
* Rue, Martino & Chopin, *Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations*, JRSS-B **71**, 319, 2009.
* Talts, Betancourt, Simpson, Vehtari & Gelman, *Validating Bayesian inference algorithms with simulation-based calibration*, arXiv:1804.06788, 2018.
* Tierney & Kadane, *Accurate approximations for posterior moments and marginal densities*, JASA **81**, 82, 1986.
* Wood, *Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models*, JRSS-B **73**, 3, 2011.
* Zuker, Szabo, Bramall, Krajcarski & Selinger, *Delta function convolution method*, Rev. Sci. Instrum. **56**, 14, 1985.
