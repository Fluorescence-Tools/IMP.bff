/**
 *  \file IMP/bff/ProbePotentialRestraints.h
 *  \brief The coarse-grained protein potentials, as IMP scores and restraints.
 *
 * **These are IMP scores, not a second scoring world.** A coarse-grained
 * potential needs coordinates, a residue lookup and a pair loop with a cutoff,
 * and IMP has all three -- particles in a #IMP::Model, an
 * #IMP::atom::Hierarchy over them, and #IMP::container::ClosePairContainer.
 * So:
 *
 *  - the typed contact potentials are #IMP::core::StatisticalPairScore
 *    subclasses, the same shape #IMP::atom::DopePairScore has, with a
 *    `add_*_score_data` that puts the residue type on the particles;
 *  - the terms IMP has no form for -- the four-channel hydrogen bond, the
 *    truncated-LJ Go model, generalized Born, the Ramachandran map -- are
 *    #IMP::Restraint subclasses;
 *  - the terms IMP *does* have a form for are built out of it rather than
 *    written again: the clash term is #IMP::core::SoftSpherePairScore, the
 *    C-alpha internal coordinates are #IMP::core::Harmonic restraints, and the
 *    radius of gyration is #IMP::atom::get_radius_of_gyration.
 *
 * The kernels are public as well, on flat arrays, because a caller holding
 * numpy should not have to build a model to score one -- and because they are
 * what the parity tests pin.
 *
 * **Distances are plain, everywhere.** Every C-alpha distance is a distance,
 * not a squared one, and each kernel squares what it needs. Mixing the two is
 * how a prefilter silently stops filtering -- the hydrogen bond's 8 A C-alpha
 * cutoff is the documented case, `okf/validation/hbond_ca_cutoff.md`.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PROBEPOTENTIALRESTRAINTS_H
#define IMPBFF_PROBEPOTENTIALRESTRAINTS_H

#include <IMP/bff/bff_config.h>

#include <IMP/Key.h>
#include <IMP/Model.h>
#include <IMP/PairScore.h>
#include <IMP/Restraint.h>
#include <IMP/atom/Atom.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Residue.h>
#include <IMP/core/StatisticalPairScore.h>
#include <IMP/core/Cosine.h>
#include <IMP/bff/RotamerScoring.h>
#include <IMP/bff/Base.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

// ---------------------------------------------------------------------------
// Typing the particles
// ---------------------------------------------------------------------------

//! The residue types the contact tables are indexed by -- their own key space.
/*! **Not #IMP::atom::ResidueType.** An IMP `Key` table grows at run time: read
    a structure with a residue IMP has not seen, or name a `ResidueType("DYE")`,
    and the next index is one further out. A statistical table is a fixed
    square, so a type index that grows past its side reads outside it -- which
    is a wrong number at best and a crash at worst, and it appears only once
    something else in the process has widened the key table.

    #IMP::atom::DopeType is a `Key` of its own for the same reason, and this is
    that key for the residue-typed tables. Its index space is closed: the names
    are the ones the table names, registered when the table is read, and a
    residue the table says nothing about is never typed. */
/*! \note The number is this key family's identity and has to be one nobody
    else uses: `783462` is #IMP::atom::ProteinLigandType, `6453462` is
    #IMP::atom::DopeType and `6453472` is `LoopStatisticalType`. Sharing one
    shares the name table, which is how this key first came back with two
    hundred names in it that were not residues. */
typedef IMP::Key<8064531> ResidueContactType;
IMP_VALUES(ResidueContactType, ResidueContactTypes);

//! The attribute the typed coarse-grained scores read a residue type from.
/*! The value is a #IMP::bff::ResidueContactType index, which is what
    #add_residue_type_score_data writes and #MiyazawaJerniganPairScore reads --
    the shape #IMP::atom::add_dope_score_data and #IMP::atom::DopePairScore
    have. */
IMPBFFEXPORT IMP::IntKey get_residue_type_key();

//! Register the residue names a table covers, so they have indices.
/*!
    The types are the table's: reading the table is what puts the names in the
    key space, and this does the same without building a score. Calling it
    twice is free.

    \param[in] table the table to take the names from, `mj` or `unres`
    \param[in] path the container; empty reads the shipped one
*/
IMPBFFEXPORT void load_residue_contact_types(std::string table = "mj",
                                             std::string path = "");

//! Put the residue type on one representative atom of every residue.
/*!
    The Miyazawa-Jernigan potential scores C-beta pairs and the UNRES centroid
    potential scores side-chain centroids; both need the particle to know which
    residue it belongs to. This is #IMP::atom::add_dope_score_data's job for
    DOPE, done for residues.

    A residue the table says nothing about -- a ligand, a modified base, a
    dye -- is **not** typed and takes no part, which is what keeps the index
    space closed. So is a residue without \p site: a glycine has no `CB`.

    \param[in] hierarchy the structure
    \param[in] site which atom stands for the residue
    \param[in] table which table's residue set to type against; empty types
               against whatever names are already registered
    \return the particles that were typed, ready for a container
    \throw ValueException when nothing could be typed, which means the table's
           names were never registered
*/
IMPBFFEXPORT IMP::ParticlesTemp add_residue_type_score_data(
        IMP::atom::Hierarchy hierarchy,
        IMP::atom::AtomType site = IMP::atom::AT_CB,
        std::string table = "mj");

// ---------------------------------------------------------------------------
// Typed contact potentials: the table is the potential
// ---------------------------------------------------------------------------

//! Miyazawa-Jernigan residue contact energy between two representative atoms.
/*!
    A contact potential: a pair inside the cutoff contributes the matrix
    element of its two residue types, and a pair outside it contributes
    nothing. That is a #IMP::core::StatisticalPairScore whose table has one
    value per type pair, so there is no loop here -- the pair loop is the
    caller's #IMP::container::ClosePairContainer and the arithmetic is IMP's.

    Score the particles #add_residue_type_score_data returned:

        ps = IMP.bff.add_residue_type_score_data(h)
        c = IMP.container.ClosePairContainer(
                IMP.container.ListSingletonContainer(m, ps), 6.5)
        r = IMP.container.PairsRestraint(IMP.bff.MiyazawaJerniganPairScore(), c)

    \see Miyazawa and Jernigan, J. Mol. Biol. 256, 623 (1996).
*/
class IMPBFFEXPORT MiyazawaJerniganPairScore
    : public IMP::core::StatisticalPairScore<IMP::bff::ResidueContactType,
                                             false, false> {
    typedef IMP::core::StatisticalPairScore<IMP::bff::ResidueContactType, false,
                                            false>
            P;

public:
    //! The shipped table.
    /*! \param[in] threshold the contact distance, A */
    MiyazawaJerniganPairScore(double threshold = 6.5);
    //! A table of the caller's, in IMP's PMF format.
    MiyazawaJerniganPairScore(double threshold, IMP::TextInput data_file);
};

//! UNRES side-chain centroid contact energy between two centroids.
/*!
    A distance-binned energy per residue-type pair -- the GBV side-chain
    potential -- which is #IMP::core::StatisticalPairScore exactly. The
    repulsion below `min_dist` and the clamp above `max_dist` are **in the
    table**: the bins below the minimum hold the penalty and the bins above the
    last measured distance repeat it, so the behaviour is the data's and not a
    branch in a loop.

    The shipped table is in kcal/mol; multiply the restraint's weight by 0.593
    for kT at 298 K, which is what chisurf's `CEPotential` did.
*/
class IMPBFFEXPORT UNRESCentroidPairScore
    : public IMP::core::StatisticalPairScore<IMP::bff::ResidueContactType,
                                             false, false> {
    typedef IMP::core::StatisticalPairScore<IMP::bff::ResidueContactType, false,
                                            false>
            P;

public:
    //! The shipped table.
    /*! \param[in] threshold the centroid distance beyond which a pair scores
        zero, A */
    UNRESCentroidPairScore(double threshold = 19.0);
    //! A table of the caller's, in IMP's PMF format.
    UNRESCentroidPairScore(double threshold, IMP::TextInput data_file);
};

//! Truncated Lennard-Jones between two beads with one equilibrium distance.
/*!
    \f$E = s(s - 2)\f$ with \f$s = (r_m^2/r^2)^3\f$: the (6,12) well in its
    \f$r_m\f$ form, with no \f$\epsilon\f$ -- the restraint's weight is the
    depth. What chisurf called the "LJ bead" potential over C-alphas.

    Unlike #IMP::atom::LennardJonesTypedPairScore this needs no types and no
    force field: one distance for every pair, which is what a one-bead-per-
    residue model has.
*/
class IMPBFFEXPORT LennardJonesBeadPairScore : public IMP::PairScore {
    double rm_;

public:
    //! \param[in] rm the equilibrium separation; the default is the C-alpha to
    //!            C-alpha distance of a trans peptide bond
    LennardJonesBeadPairScore(double rm = 3.8208650279);

    double get_rm() const { return rm_; }
    void set_rm(double rm) { rm_ = rm; }

    virtual double evaluate_index(
            IMP::Model* m, const IMP::ParticleIndexPair& p,
            IMP::DerivativeAccumulator* da) const override;
    virtual IMP::ModelObjectsTemp do_get_inputs(
            IMP::Model* m, const IMP::ParticleIndexes& pis) const override;
    IMP_OBJECT_METHODS(LennardJonesBeadPairScore);
};

// ---------------------------------------------------------------------------
// The terms IMP has no form for
// ---------------------------------------------------------------------------

//! The native contacts of a conformer: a well depth and a minimum, per pair.
struct IMPBFFEXPORT GoContacts {
    //! The two residues of each contact, parallel to #energies and #distances.
    std::vector<int> residue_i, residue_j;
    //! The well depth of each contact.
    std::vector<double> energies;
    //! The equilibrium distance of each contact, A.
    std::vector<double> distances;
    //! How many of them were inside the native cutoff.
    int n_native;

    GoContacts() : n_native(0) {}
    int get_n_contacts() const { return static_cast<int>(energies.size()); }

    void get_energies(double** out_view, int* n_out_view) const;
    void get_distances(double** out_view, int* n_out_view) const;

    IMP_SHOWABLE_INLINE(GoContacts, out << "GoContacts(" << energies.size()
                                        << " pairs, " << n_native
                                        << " native)");
};
IMP_VALUES(GoContacts, GoContactsList);

//! A Go model over a set of beads: truncated Lennard-Jones on native contacts.
/*!
    Every pair gets a well whose minimum is the distance it has **in the
    structure the restraint was built from** -- that is what makes it a Go
    model -- with the full depth inside \p cutoff and \p nn_e_factor of it
    outside. The energy is

    \f$E = \sum \epsilon_{ij}(-s + s^2/4 + 0.00818)\f$, \f$s = 2(r_m/r)^6\f$,

    inside \f$r_m < r < 2.5 r_m\f$ and \f$-\epsilon_{ij}\f$ outside it.

    \note #IMP::bff::create_go_restraints in `RotamerScoring.h` is a *harmonic* native
    contact over force-field sites. This is the truncated-LJ form over beads;
    they are different terms and both stay.

    \note The Python in chisurf carried its accumulator across pairs
    (`tmp -= sr` on the previous pair's value), so every energy after the first
    was wrong by whatever came before it. The imp-tricks kernel fixed that with
    a local, and this is that one.
*/
class IMPBFFEXPORT GoRestraint : public IMP::Restraint {
    IMP::ParticleIndexes pis_;
    double epsilon_, cutoff_, nn_e_factor_;
    GoContacts contacts_;

public:
    //! \param[in] m,pis the beads, one per residue
    /*! \param[in] epsilon the native well depth
        \param[in] cutoff the native contact distance, A
        \param[in] nn_e_factor the non-native fraction of \p epsilon; 0
                   switches non-native contacts off entirely */
    GoRestraint(IMP::Model* m, const IMP::ParticleIndexes& pis,
                double epsilon = 1.0, double cutoff = 6.5,
                double nn_e_factor = 0.7);

    //! Retake the native contacts from where the beads are now.
    void set_native_contacts();
    const GoContacts& get_contacts() const { return contacts_; }
    int get_n_native() const { return contacts_.n_native; }
    int get_n_non_native() const;

    virtual double unprotected_evaluate(
            IMP::DerivativeAccumulator* accum) const override;
    virtual IMP::ModelObjectsTemp do_get_inputs() const override;
    IMP_OBJECT_METHODS(GoRestraint);
};

//! Backbone hydrogen bonds, from a four-channel distance lookup.
/*!
    A bond is an amide H within \p cutoff_h of a carbonyl O; each one reads
    four distances -- O-H, O-N, C-H and C-N -- out of the four channels of the
    table. Both directions of a residue pair are tested, so an antiparallel
    pair can contribute twice. The channels can be switched off one at a time,
    which is what chisurf's `oh`/`on`/`cn`/`ch` flags did.

    The residues' N, C, O and H are resolved once, when the restraint is built.

    \note The C-alpha cutoff **is applied**, unlike in the kernel this was
    ported from; see `okf/validation/hbond_ca_cutoff.md`.
*/
class IMPBFFEXPORT HydrogenBondRestraint : public IMP::Restraint {
    // Four atoms per residue, -1 where the residue has none: N, C, O, H.
    std::vector<int> n_, c_, o_, h_;
    IMP::ParticleIndexes pis_;
    std::vector<double> table_;
    int n_bins_;
    double cutoff_ca_, cutoff_h_, bin_width_;
    bool ch_, on_, oh_, cn_;
    mutable int n_hbonds_;

public:
    //! \param[in] m,hierarchy the structure; its residues are walked once
    /*! \param[in] table the lookup, flat `4 * n_bins`, channels CH, ON, OH, CN
        \param[in] n_bins bins per channel; 0 takes it from the table's length
        \param[in] cutoff_ca residue pairs further apart than this are not
                   examined, A
        \param[in] cutoff_h the O to H distance that makes a bond, A
        \param[in] bin_width the lookup's bin width, A */
    HydrogenBondRestraint(IMP::Model* m, IMP::atom::Hierarchy hierarchy,
                          const std::vector<double>& table, int n_bins = 0,
                          double cutoff_ca = 8.0, double cutoff_h = 3.0,
                          double bin_width = 0.01);

    //! Which of the four channels count.
    void set_channels(bool ch, bool on, bool oh, bool cn);
    //! How many bonds the last evaluation found.
    int get_n_hbonds() const { return n_hbonds_; }

    virtual double unprotected_evaluate(
            IMP::DerivativeAccumulator* accum) const override;
    virtual IMP::ModelObjectsTemp do_get_inputs() const override;
    IMP_OBJECT_METHODS(HydrogenBondRestraint);
};

//! Generalized-Born implicit solvation over a set of atoms.
/*!
    \f$E = \frac{1}{8\pi}(\frac{1}{\epsilon_0} - \frac{1}{\epsilon})
    \sum_{i \le j} q_i q_j (\frac{1}{f_{GB}} - \frac{1}{f_c})\f$ with
    \f$f_{GB} = \sqrt{r^2 + a_i a_j e^{-r^2/4a_ia_j}}\f$. The sum includes each
    atom with itself, which is where the Born self-energy comes from -- so this
    is a restraint over a set and not a #IMP::PairScore over pairs.

    Charges come from #IMP::atom::Charged and radii from #IMP::core::XYZR,
    which is where #IMP::atom::CoulombPairScore reads them too.

    \note The shift that smooths the cutoff uses `cutoff` where the squared
    distance would be dimensionally right (`sqrt(cutoff + ...)`, not
    `sqrt(cutoff^2 + ...)`). That is what the reference implementations do,
    and changing it would move every energy this has ever produced; it is kept
    deliberately, and said out loud here.
*/
class IMPBFFEXPORT GeneralizedBornRestraint : public IMP::Restraint {
    IMP::ParticleIndexes pis_;
    double epsilon_, epsilon0_, cutoff_;

public:
    //! \param[in] m,pis the charged atoms
    /*! \param[in] epsilon the solute dielectric constant
        \param[in] epsilon0 the solvent dielectric constant
        \param[in] cutoff pairs beyond this are dropped, A */
    GeneralizedBornRestraint(IMP::Model* m, const IMP::ParticleIndexes& pis,
                             double epsilon = 4.0, double epsilon0 = 80.1,
                             double cutoff = 12.0);

    virtual double unprotected_evaluate(
            IMP::DerivativeAccumulator* accum) const override;
    virtual IMP::ModelObjectsTemp do_get_inputs() const override;
    IMP_OBJECT_METHODS(GeneralizedBornRestraint);
};

//! The Ramachandran pseudo-energy of a chain's backbone dihedrals.
/*!
    \f$\phi_i\f$ from `C(i-1), N(i), CA(i), C(i)` and \f$\psi_i\f$ from
    `N(i), CA(i), C(i), N(i+1)`; the first residue has no \f$\phi\f$ and the
    last no \f$\psi\f$, and both are skipped. Each residue's
    \f$(\phi, \psi)\f$ is read out of a grid as \f$-\log(P/P_{max})\f$, so a
    uniform map scores zero. A map whose values never go positive is taken to
    hold pseudo-energies already and is returned negated.

    chisurf's class of this name returned 0.0 and logged a warning -- the C
    function it called had been lost. This is the imp-tricks implementation.
*/
class IMPBFFEXPORT RamachandranRestraint : public IMP::Restraint {
    // Four atoms per residue: C(i-1), N, CA, C, and N(i+1); -1 where absent.
    std::vector<int> c_prev_, n_, ca_, c_, n_next_;
    std::vector<int> channel_;
    IMP::ParticleIndexes pis_;
    std::vector<double> grids_;
    int n_channels_, n_bins_;
    double empty_penalty_;

public:
    //! \param[in] m,hierarchy the structure
    /*! \param[in] grids flat `n_channels * n_bins * n_bins`, row-major in phi
        \param[in] n_channels one grid per residue class; the shipped map has
                   three, general / proline / glycine
        \param[in] n_bins the side of one grid; 0 derives it from the length
        \param[in] empty_penalty what a zero cell costs */
    RamachandranRestraint(IMP::Model* m, IMP::atom::Hierarchy hierarchy,
                          const std::vector<double>& grids, int n_channels = 3,
                          int n_bins = 0, double empty_penalty = 10.0);

    //! The backbone dihedrals as they stand, radians; NaN at the chain ends.
    void get_phi(double** out_view, int* n_out_view) const;
    void get_psi(double** out_view, int* n_out_view) const;

    virtual double unprotected_evaluate(
            IMP::DerivativeAccumulator* accum) const override;
    virtual IMP::ModelObjectsTemp do_get_inputs() const override;
    IMP_OBJECT_METHODS(RamachandranRestraint);
};

// ---------------------------------------------------------------------------
// Built out of IMP's own scores
// ---------------------------------------------------------------------------

//! The soft-sphere clash restraint over a hierarchy's atoms.
/*!
    #IMP::core::SoftSpherePairScore is \f$\frac{1}{2}k(\sigma - r)^2\f$ on
    overlap, and the potential this ports is \f$((\sigma - r)/t)^2\f$, so
    \f$k = 2/t^2\f$ -- the same term, and IMP's, evaluated over an
    #IMP::container::ClosePairContainer rather than every pair.

    The Python skipped pairs closer than a `covalent_radius`, meaning "these
    are bonded, do not count them". Bonded pairs are excluded here by
    #IMP::atom::StereochemistryPairFilter when the hierarchy carries bonds,
    which is what that parameter was reaching for and is exact rather than a
    distance guess.

    \param[in] hierarchy the atoms
    \param[in] clash_tolerance the \f$t\f$ above: larger is softer
    \param[in] slack the container's slack, A
    \return one restraint, named `clash`
*/
IMPBFFEXPORT IMP::Restraint* create_clash_restraint(
        IMP::atom::Hierarchy hierarchy, double clash_tolerance = 2.0,
        double slack = 2.0);

//! Harmonic restraints holding a C-alpha trace near a reference conformer.
/*!
    The bond lengths, angles and dihedrals of consecutive C-alphas, each with
    an #IMP::core::Harmonic at the value the reference has. chisurf's
    `internal_potential` summed \f$k(x - x_0)^2\f$ over the three; an IMP
    harmonic scores \f$\frac{1}{2}k(x-x_0)^2\f$, so the constants here are
    twice the reference's and the energy is the same.

    \note #IMP::atom::CAAngleRestraint and #IMP::atom::CADihedralRestraint are
    IMP's *statistical* C-alpha terms -- a score per angle bin. These are the
    harmonic-to-a-reference kind, which is a different question.

    \param[in] m the model
    \param[in] cas the C-alpha particles, in sequence
    \param[in] reference the conformer to hold them near, same length
    \param[in] k_bond,k_angle,k_dihedral the force constants; a zero omits its
               family entirely
    \return the restraints, bonds then angles then dihedrals
*/
IMPBFFEXPORT IMP::Restraints create_ca_internal_restraints(
        IMP::Model* m, const IMP::ParticleIndexes& cas,
        const IMP::ParticleIndexes& reference, double k_bond = 1.0,
        double k_angle = 0.2, double k_dihedral = 0.1);

//! Residue-level solvent-accessible surface area of a hierarchy, A^2.
/*!
    One representative sphere per residue -- the Shrake-Rupley construction
    with a single site standing for the whole residue, which is what a
    coarse-grained model has. Composed from #IMP::bff::sphere_points and
    #IMP::bff::solvent_accessible_surface_area rather than written again.

    \param[in] hierarchy the structure
    \param[in] site which atom stands for the residue
    \param[in] n_sphere samples per sphere
    \param[in] probe the probe radius, A
    \param[in] radius the radius every representative is given, A

    \note The kernel this was ported from read its representative out of column
    6 of a residue lookup table, which is `CG` in the layout it was fed -- so it
    reported the area of the residues that have a `CG` and nothing else.
*/
IMPBFFEXPORT double residue_solvent_accessible_surface(
        IMP::atom::Hierarchy hierarchy,
        IMP::atom::AtomType site = IMP::atom::AT_CA, int n_sphere = 590,
        double probe = 1.0, double radius = 2.5);

// ---------------------------------------------------------------------------
// The kernels
// ---------------------------------------------------------------------------

//! Soft-sphere overlap energy over every atom pair, on flat coordinates.
/*!
    \f$E = \sum_{i<j} ((\sigma_{ij} - r_{ij})/t)^2\f$ over the pairs with
    \f$r_c < r_{ij} < \sigma_{ij}\f$. What #create_clash_restraint scores
    through IMP, for a caller holding arrays.

    \param[in] xyz flat coordinates, three per atom
    \param[in] vdw one radius per atom
    \param[in] clash_tolerance the denominator \f$t\f$
    \param[in] covalent_radius pairs closer than this are not counted
*/
IMPBFFEXPORT double clash_energy(const std::vector<double>& xyz,
                                 const std::vector<double>& vdw,
                                 double clash_tolerance = 2.0,
                                 double covalent_radius = 1.5);

//! The native contacts of a set of beads, from their distances.
IMPBFFEXPORT GoContacts go_native_contacts(const std::vector<double>& xyz,
                                           double epsilon,
                                           double nn_e_factor = 0.1,
                                           double cutoff = 6.5);

//! The truncated Lennard-Jones energy of a set of beads over its contacts.
IMPBFFEXPORT double go_energy(const std::vector<double>& xyz,
                              const GoContacts& contacts);

//! Truncated Lennard-Jones over beads, on flat coordinates.
IMPBFFEXPORT double lennard_jones_bead_energy(
        const std::vector<double>& xyz, double rm = 3.8208650279);

//! Generalized-Born solvation energy, on flat coordinates.
IMPBFFEXPORT double generalized_born_energy(const std::vector<double>& xyz,
                                            const std::vector<double>& radii,
                                            const std::vector<double>& charges,
                                            double epsilon = 4.0,
                                            double epsilon0 = 80.1,
                                            double cutoff = 12.0);

//! Residue-level solvent-accessible surface area, on flat coordinates.
/*! \param[in] xyz the representative atoms' coordinates, three per residue */
IMPBFFEXPORT double residue_asa(const std::vector<double>& xyz,
                                int n_sphere = 590, double probe = 1.0,
                                double radius = 2.5);

//! The Ramachandran pseudo-energy of one \f$(\phi, \psi)\f$, from a grid.
/*!
    \param[in] phi,psi radians; NaN (a chain end) scores zero
    \param[in] grid one channel, flat `n_bins * n_bins`, row-major in phi
    \param[in] n_bins the side of it
    \param[in] empty_penalty what a zero cell costs
*/
IMPBFFEXPORT double ramachandran_energy(double phi, double psi,
                                        const std::vector<double>& grid,
                                        int n_bins,
                                        double empty_penalty = 10.0);

// -------- from RotamerScoring.h (the restraint factories over a typed dye system) --------
// These are the IMP side of RotamerScoring.h: everything above in that header is
// arithmetic over arrays and typed systems and stays in the core; these five
// return IMP restraints and scores over IMP particles, so they live here in
// the connection layer (PRD-137 step 5).

//! The \c IMP::core::Cosine a torsion type spells, converted from CHARMM.
/*!
    cgprobe stores torsions in the CHARMM convention
    \f$V = k(1 + \cos(n\phi - \delta))\f$, and \c IMP::core::Cosine scores
    \f$k(1 - \cos(n\phi - \delta'))\f$ -- so the phase needs
    \f$\delta' = \delta + \pi\f$. Getting that wrong put every conjugated
    torsion's minimum at 90 degrees (PRD-108); it is one conversion, in one
    place, for that reason.
*/
IMPBFFEXPORT IMP::core::Cosine* torsion_cosine(const FFTorsionType& type);

//! Bonded and element-aware nonbonded restraints for a typed dye system.
/*!
    Bonds and angles become harmonics, torsions the CHARMM cosine above,
    impropers a harmonic about the *current* geometry's dihedral, and every
    non-excluded site pair a lower-bound harmonic at its Lennard-Jones
    \f$R_{min}\f$. A site with no particle is skipped rather than an error:
    a system may describe more than the caller decorated.

    \param[in] model the model the restraints score in
    \param[in] system the typed system (bonds, angles, torsions, impropers)
    \param[in] site_ids,particles the decorated sites, parallel; a site id
               that appears twice takes its last particle
    \param[in] nonbonded add the repulsion (#create_steric_restraint). False
               leaves it to the caller, which is what a run wants when it
               scores rigid moves against the repulsion *alone*.
    \throw ValueException when the two sequences disagree in length

    \note A bond whose equilibrium length is not set (zero, which no real bond
    has) is restrained about the geometry as it stands, and an angle likewise.
    That is what a system built without a template carries, and refusing it
    would refuse the systems this module builds itself.
*/
IMPBFFEXPORT IMP::Restraints create_probe_restraints(
        IMP::Model* model, const ProbeForceFieldSystem& system,
        const std::vector<std::string>& site_ids,
        const IMP::ParticleIndexes& particles, bool nonbonded = true);

//! The steric term alone: one soft-sphere restraint over the whole system.
/*!
    Every site pair the system does not exclude, in a single
    #IMP::container::PairsRestraint over a #IMP::core::SoftSpherePairScore.

    **This is the repulsion, for molecular dynamics and Monte Carlo alike.**
    It is differentiable, so a dynamics run can use it; it is one restraint
    over a container rather than one per pair, so a long run can afford it;
    and it depends only on where the spheres are, so a rigid move can be
    scored against it and nothing else -- a rigid move cannot change a bond,
    an angle or a torsion, and scoring those during one is work whose answer
    never changes. One repulsion serves both: soft spheres for the Monte-Carlo
    step *and* for dynamics, rather than these plus per-pair Lennard-Jones
    lower bounds elsewhere, which would be two implementations of one piece of
    physics.

    The Lennard-Jones parameters remain what
    #IMP::bff::IntramolecularEnergy *evaluates*: an energy of a
    conformation is a different question from keeping two atoms apart, and a
    sphere overlap answers the second.

    The exclusions are the system's own (#ProbeForceFieldSystem::get_exclusions),
    so this and the bonded terms cannot disagree about which pairs are 1-2,
    1-3 or 1-4.

    \param[in] model,system,site_ids,particles as for #create_probe_restraints
    \param[in] k the soft-sphere force constant; negative takes the system's
               own `nonbonded.k`
    \return the restraint; null when the system's non-bonded term is off or
            it has no non-excluded pair
*/
IMPBFFEXPORT IMP::Restraint* create_steric_restraint(
        IMP::Model* model, const ProbeForceFieldSystem& system,
        const std::vector<std::string>& site_ids,
        const IMP::ParticleIndexes& particles, double k = -1.0);

//! Place a guest rigidly around a host and keep the best-scoring pose.
/*!
    A random search, not an optimisation: each trial rotates the guest about a
    random axis by a random angle and puts its centre a fixed distance from
    the host's, in a random direction; the pose that scores lowest is the one
    left in the model. That is enough to start a simulation somewhere
    plausible rather than wherever the input files happened to put the two
    components -- which for a dye and a protein read from separate files is
    usually on top of each other.

    \param[in] scoring_function what a trial is judged by
    \param[in] model the model both sets of particles belong to
    \param[in] host,guest the two particle sets; only \p guest moves
    \param[in] distance the separation of the two centres, A
    \param[in] n_trials how many poses to try
    \param[in] seed the random seed, so a run repeats
    \return the best score, or NaN when there is nothing to place

    The guest keeps its own shape throughout: every trial applies one rotation
    and one translation to the same input coordinates, so this cannot distort
    what it is placing.
*/
IMPBFFEXPORT double place_guest_by_score(
        IMP::ScoringFunction* scoring_function, IMP::Model* model,
        const IMP::ParticleIndexes& host, const IMP::ParticleIndexes& guest,
        double distance = 62.5, int n_trials = 200, int seed = 42);

//! Native-contact (Go) restraints that hold a component in its own shape.
/*!
    Every heavy-atom pair of a component that is closer than \p cutoff and is
    not already within two bonds gets a harmonic at the distance it currently
    has. That is what keeps a component's fold while its parts are free to
    move: the bonded terms fix the local geometry, these fix the tertiary one.

    Hydrogens are left out (a name beginning with `H`): they add pairs without
    adding shape, and they are the majority of the atoms.

    \param[in] model,system the model and the typed system
    \param[in] site_ids,particles the decorated sites, parallel
    \param[in] site_atom_names the atom name of each site id, for the hydrogen
               test; a site missing from it counts as heavy
    \param[in] component which component to restrain
    \param[in] only_sites when non-empty, a pair is kept only if at least one
               of its two sites is in this set -- how a mostly-rigid component
               releases a few atoms and restrains only what they touch
    \param[in] k the force constant
    \param[in] cutoff the contact distance, A
    \return one restraint per contact, named `go_<component>_<a>_<b>`

    \note The equilibrium is `max(d, 1 A)`: two atoms that a structure has
    placed on top of each other would otherwise get a harmonic at zero, which
    is a singularity a minimiser walks straight into.
*/
IMPBFFEXPORT IMP::Restraints create_go_restraints(
        IMP::Model* model, const ProbeForceFieldSystem& system,
        const std::vector<std::string>& site_ids,
        const IMP::ParticleIndexes& particles,
        const std::map<std::string, std::string>& site_atom_names,
        const std::string& component,
        const std::vector<std::string>& only_sites = std::vector<std::string>(),
        double k = 3.0, double cutoff = 6.0);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_PROBEPOTENTIALRESTRAINTS_H
