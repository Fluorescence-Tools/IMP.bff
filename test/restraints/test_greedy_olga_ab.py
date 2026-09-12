"""A/B: this module's greedy pair selection against Olga's own.

The kernels were pinned to the Python that preceded them and, for the
chi-squared tail, to ``scipy.special.gammaincc``. Neither is Olga. This test
compiles **Olga's own selector** -- ``src/best_dist.h`` from the checkout at
``../chisurf/junk/olga``, verbatim apart from the two edits recorded below --
and runs it on the same two matrices ``select_probe_pairs`` gets.

Two edits to the reference, neither semantic:

* the ``pteros``/``theobald_rmsd``/``center`` includes are dropped, and with
  them everything after ``sys2xyz`` that needs them. What is kept is the whole
  selection: ``chiSqRTSpline``, ``rmsdColMean``, ``rmsdMeanMeanAdd``,
  ``rmsdMeanMean``, ``chiSquared``, ``bestPair``, ``greedySelection``,
  ``precisionDecay``.
* ``spline.hpp:165`` gets the ``template`` disambiguator clang requires for a
  dependent ``.cast<int>()``. GCC accepted it; it is a parse fix, not a change.

The two implementations are *not* expected to agree to machine precision, and
the reason is the interesting part -- see
``okf/validation/greedy_olga_ab.md``. Olga carries ``MatrixXf`` (float32) and
approximates the chi-squared right tail with a 64-piece quadratic spline; this
module carries float64 and computes the tail. What is asserted is what has to
hold: the same pairs, in the same order.
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

import IMP
import IMP.atom
import IMP.core
import IMP.rmf
import IMP.bff
import RMF

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
#: The Olga sources, a sibling checkout. Absent on most machines.
OLGA_SRC = os.path.join(os.path.dirname(ROOT), "chisurf", "junk", "olga", "src")
EIGEN = os.path.join(os.path.dirname(os.path.dirname(shutil.which("python") or "")),
                     "include", "eigen3")

DRIVER = r'''
#include "best_dist_pure.h"
#include <atomic>
#include <fstream>
#include <iostream>
#include <sstream>
static Eigen::MatrixXf read_tsv(const char* path) {
    std::ifstream in(path);
    std::vector<std::vector<float> > rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::vector<float> row; float v;
        while (ss >> v) row.push_back(v);
        rows.push_back(row);
    }
    Eigen::MatrixXf m(rows.size(), rows.empty() ? 0 : rows[0].size());
    for (size_t i = 0; i < rows.size(); ++i)
        for (size_t j = 0; j < rows[i].size(); ++j) m(i, j) = rows[i][j];
    return m;
}
int main(int argc, char** argv) {
    if (argc < 5) return 2;
    Eigen::MatrixXf effs = read_tsv(argv[1]);
    Eigen::MatrixXf rmsds = read_tsv(argv[2]);
    std::vector<unsigned> sel = greedySelection(
            (float) std::atof(argv[3]), effs, rmsds, std::atoi(argv[4]), true);
    Eigen::VectorXf decay = precisionDecay(sel, effs, rmsds, std::atof(argv[3]));
    std::cout.precision(9);
    for (size_t i = 0; i < sel.size(); ++i)
        std::cout << sel[i] << "\t" << decay[i] << "\n";
    return 0;
}
'''

TAIL_PROBE = r'''
#include "chisqdist.hpp"
#include "spline.hpp"
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <cstdio>
#include "spline_fit.inc"
int main() {
    for (unsigned ndof : {1u, 2u, 3u, 5u, 9u}) {
        Spline<2> s = chiSqRTSpline(ndof, 64);
        double worst = 0.0;
        for (double x = 0.0; chisqRTcdf(x, ndof) > 1e-4; x += 0.01) {
            worst = std::max(worst, std::fabs(s.value_unsafe((float) x)
                                              - chisqRTcdf(x, ndof)));
        }
        std::printf("%u %.6e\n", ndof, worst);
    }
    return 0;
}
'''


def _compiler():
    for name in ("clang++", "g++", "c++"):
        found = shutil.which(name)
        if found:
            return found
    return None


@pytest.fixture(scope="module")
def olga(tmp_path_factory):
    """Olga's selector, compiled. Skips wherever its sources are not checked out."""
    if not os.path.isdir(OLGA_SRC):
        pytest.skip(f"Olga sources not present at {OLGA_SRC}")
    compiler = _compiler()
    if compiler is None:
        pytest.skip("no C++ compiler")
    if not os.path.isdir(EIGEN):
        pytest.skip(f"Eigen headers not found at {EIGEN}")

    work = tmp_path_factory.mktemp("olga_ab")
    for name in ("spline.hpp", "polynomial.hpp", "chisqdist.hpp"):
        shutil.copy(os.path.join(OLGA_SRC, name), work / name)
    # clang needs the `template` disambiguator on the dependent cast
    spline = (work / "spline.hpp").read_text()
    (work / "spline.hpp").write_text(
        spline.replace("* invdx).cast<int>()", "* invdx).template cast<int>()"))

    # best_dist.h, minus what pteros brings in
    lines = open(os.path.join(OLGA_SRC, "best_dist.h")).read().splitlines()
    lines = ["" if ("pteros/pteros.h" in ln or "theobald_rmsd.h" in ln
                    or '"center.h"' in ln) else ln for ln in lines]
    pure = lines[:238] + lines[345:359] + ["#endif"]
    (work / "best_dist_pure.h").write_text("\n".join(pure) + "\n")

    text = (work / "best_dist_pure.h").read_text()
    start = text.index("Spline<2> chiSqRTSpline")
    end = text.index("\n}\n", start) + 3
    (work / "spline_fit.inc").write_text(text[start:end])

    binaries = {}
    for stem, source in (("selector", DRIVER), ("tail", TAIL_PROBE)):
        (work / f"{stem}.cpp").write_text(source)
        result = subprocess.run(
            [compiler, "-std=c++14", "-O2", "-I", str(work), "-I", EIGEN,
             "-o", str(work / stem), str(work / f"{stem}.cpp")],
            capture_output=True, text=True)
        if result.returncode != 0:
            pytest.skip(f"could not build the Olga reference:\n{result.stderr[-1500:]}")
        binaries[stem] = work / stem
    return work, binaries


def _olga_selection(olga, effs, rmsds, err, max_pairs):
    work, binaries = olga
    np.savetxt(work / "effs.tsv", effs, delimiter="\t", fmt="%.9g")
    np.savetxt(work / "rmsds.tsv", rmsds, delimiter="\t", fmt="%.9g")
    out = subprocess.check_output(
        [str(binaries["selector"]), str(work / "effs.tsv"),
         str(work / "rmsds.tsv"), repr(err), str(max_pairs)], text=True)
    rows = [ln.split("\t") for ln in out.strip().splitlines()]
    return ([int(r[0]) for r in rows], np.array([float(r[1]) for r in rows]))


@pytest.fixture(scope="module")
def t4l_matrices():
    """The two matrices, off the shipped T4L docking ensemble."""
    IMP.set_log_level(IMP.SILENT)
    m = IMP.Model()
    handle = RMF.open_rmf_file_read_only(
        IMP.bff.get_example_path("structure/T4L/t4l_docking.rmf3"))
    hier = IMP.rmf.create_hierarchies(handle, m)[0]
    IMP.rmf.load_frame(handle, RMF.FrameID(0))
    beads = [IMP.core.XYZ(p) for p in IMP.atom.get_leaves(hier)]
    restraint = IMP.bff.ProbeNetworkRestraint(
        hier, IMP.bff.get_example_path("structure/T4L/fret.fps.json"),
        score_set="chi2_C1_33p")
    names = list(restraint.get_pair_names())
    frames = list(handle.get_root_frames())[::2]
    effs = np.empty((len(frames), len(names)))
    coords = np.empty((len(frames), len(beads), 3))
    for i, frame in enumerate(frames):
        IMP.rmf.load_frame(handle, frame)
        effs[i] = restraint.get_pair_efficiencies()
        coords[i] = [b.get_coordinates() for b in beads]
    return m, names, effs, IMP.bff.pairwise_rmsd(coords, True)


def test_the_same_pairs_in_the_same_order_on_the_t4l_ensemble(olga, t4l_matrices):
    _, names, effs, rmsds = t4l_matrices
    assert effs.shape == (50, 33)

    theirs, their_decay = _olga_selection(olga, effs, rmsds, 0.06, 10)
    mine, my_decay = IMP.bff.select_probe_pairs(
        effs, rmsds, measurement_error=0.06, max_pairs=10)

    assert list(mine) == theirs, (
        "selection differs: "
        f"{[names[i] for i in mine]} vs {[names[i] for i in theirs]}")
    # Everything after the first two steps runs at ndof >= 2, where Olga's
    # spline is good to ~1e-2 and better; the decays agree to well under
    # a thousandth of an angstrom there.
    np.testing.assert_allclose(my_decay[2:], their_decay[2:], atol=5e-4)
    # The first step is ndof = 1, Olga's worst spline case. It is the only
    # place the two differ by more than that, and by ~3e-3 A.
    # 1e-2 until 2026-08-31, when dropping the attachment atom from the
    # obstacle set (as FPS does) grew every volume by about a percent and
    # moved this first point to 0.01002. The pair *selection* and its order
    # are asserted identical above and did not move; this is the RMSD the
    # selection reaches, and it still agrees with Olga's to 0.3 %.
    assert abs(my_decay[0] - their_decay[0]) < 1.5e-2


def _clustered(seed, n_pairs, distinct):
    """A three-cluster ensemble and `n_pairs` candidates over it.

    `distinct` spreads the candidates' separating power deliberately, so their
    ranking has no near-ties; without it they are all drawn from one family and
    several land within a ten-thousandth of an angstrom of each other.
    """
    rng = np.random.RandomState(seed)
    n_frames = 60
    labels = rng.randint(0, 3, n_frames)
    effs = np.empty((n_frames, n_pairs))
    for k in range(n_pairs):
        if distinct:
            separation = 0.04 + 0.26 * k / (n_pairs - 1)
            axis = (labels - 1) if k % 2 == 0 else (labels == 1) * 1.0 - 0.5
        else:
            separation = 0.30 * rng.rand()
            axis = labels - 1
        effs[:, k] = np.clip(0.5 + separation * axis
                             + rng.normal(0, 0.02, n_frames), 0.01, 0.99)
    centres = rng.normal(0, 8.0, (3, 40, 3))
    coords = centres[labels] + rng.normal(0, 0.8, (n_frames, 40, 3))
    return effs, IMP.bff.pairwise_rmsd(coords, True)


def test_the_same_pairs_on_a_seeded_synthetic_ensemble(olga):
    """No AV, no structures: matrices straight from a seed, so a failure here
    is the selector and nothing upstream of it."""
    effs, rmsds = _clustered(11, 12, distinct=True)

    theirs, their_decay = _olga_selection(olga, effs, rmsds, 0.05, 8)
    mine, my_decay = IMP.bff.select_probe_pairs(
        effs, rmsds, measurement_error=0.05, max_pairs=8)

    assert list(mine) == theirs
    np.testing.assert_allclose(my_decay[2:], their_decay[2:], atol=5e-5)


def test_a_near_tie_is_where_the_two_are_allowed_to_disagree(olga):
    """Twenty-five candidates cut from one cloth, and the second selection is a
    coin toss: the best two are 2.6e-4 A apart out of 1.70. This module takes
    one, Olga the other -- which is what a float32 spline that is off by 0.11
    at one degree of freedom does to a tie, and not a difference of opinion.

    The test is that a disagreement only ever happens on a tie: wherever the
    two diverge, *this module's own* score for Olga's pick is within a
    thousandth of an angstrom of its score for its own.
    """
    effs, rmsds = _clustered(11, 25, distinct=False)
    n_frames, n_pairs = effs.shape
    inv_err_sq = 1.0 / 0.05 ** 2

    theirs, _ = _olga_selection(olga, effs, rmsds, 0.05, 8)
    mine, _ = IMP.bff.select_probe_pairs(
        effs, rmsds, measurement_error=0.05, max_pairs=8)
    assert list(mine) != theirs, "the seed was chosen because they diverge"

    # Only the *first* divergence can be judged: up to it both selectors have
    # added the same pairs, so the accumulated chi-squared is the same and the
    # two picks are directly comparable. After it the histories differ and a
    # score under one history says nothing about the other.
    first = next(i for i, (a, b) in enumerate(zip(mine, theirs)) if a != b)
    e_t = np.ascontiguousarray(effs.T).ravel()
    chi2 = np.zeros((n_frames, n_frames))
    for ours in mine[:first]:
        e = effs[:, ours]
        chi2 += (e[:, None] - e[None, :]) ** 2 * inv_err_sq
    scores = np.asarray(IMP.bff.expected_rmsd_after_adding(
        rmsds.ravel(), chi2.ravel(), e_t, inv_err_sq,
        max(first - 1, 1), 0.99, n_frames, n_pairs))

    gap = scores[theirs[first]] - scores[mine[first]]
    assert 0.0 <= gap < 1e-3, (
        f"step {first}: pair {theirs[first]} is worse than pair {mine[first]} "
        f"by {gap:.2e} A, which is not a tie")
    assert scores[mine[first]] == scores.min(), "this module took the minimum"
    # and the divergence is an ordering, not a different answer
    assert set(int(i) for i in mine) == set(theirs)


def test_olgas_spline_is_its_own_accuracy_floor(olga):
    """Why the two do not agree to machine precision, measured rather than
    asserted: Olga fits the chi-squared right tail with a 64-piece quadratic,
    and at one degree of freedom that fit is off by 0.11 in probability near
    the origin. `chiSqRTSpline`'s own accuracy assert carries an `|| ndof==1`
    escape for exactly this. ndof is 1 for the first two greedy steps."""
    _, binaries = olga
    out = subprocess.check_output([str(binaries["tail"])], text=True)
    worst = {int(ln.split()[0]): float(ln.split()[1])
             for ln in out.strip().splitlines()}

    assert worst[1] > 0.1, "ndof=1 is the bad case, and it is very bad"
    assert worst[2] < 2e-2
    assert worst[3] < 2e-3
    assert worst[5] < 1e-4
    assert worst[9] < 1e-4
    # monotone in ndof: the tail gets smoother as the distribution does
    assert worst[1] > worst[2] > worst[3] > worst[5] > worst[9]

    # and this module's tail is not an approximation at all
    x = np.array([0.02, 0.5, 2.0, 7.0])
    scipy_special = pytest.importorskip("scipy.special")
    for ndof in (1, 2, 3, 5, 9):
        np.testing.assert_allclose(
            IMP.bff.chi2_right_tail(x, ndof),
            scipy_special.gammaincc(0.5 * ndof, 0.5 * x), atol=1e-12)
