#ifndef MCNUX_ESCAPE_HXX
#define MCNUX_ESCAPE_HXX

// Per-species escape tallies at the outer domain boundary
// ([MCNX-GPU-05] of specs/particle-container-and-gpu.md: "Packets whose
// positions have left the outer domain boundary are removed and tallied as
// escapes (per the `README.md` boundary policy) rather than redistributed").
//
// Pure arithmetic only: the escape predicate, the tally slot map, and the
// cumulative accumulator type. The runtime detection/removal lives at the
// two live silent-drop sites — the geodesic push (mcnux_geodesic.cxx) and
// the episode driver (mcnux_interactions.cxx) — which tally into a device
// buffer, invalidate the escaped packet (`make_invalid`, the absorption
// idiom; the id VALUE stays retired per [MCNX-GPU-04]), and fold into the
// run accumulator behind escape_run_tally() (defined in mcnux_escape.cxx).
// The emission loop's Redistribute() is PROVABLY never an escape site
// (creation positions land in the walked cell's own half-open interval —
// see the rationale at the mcnux_emission.cxx Redistribute call), so it
// carries no escape code.
//
// Storage-surface decision (recorded per the escape-tallies task): escapes
// are surfaced through the gated diagnostic grid array
// MCNuX::mcnux_escape_diag (SIZE = escape_num_slots, the mcnux_ledger_diag /
// mcnux_stats_diag DISTRIB=CONSTANT precedent), fed from the cumulative
// accumulator here. They are deliberately NOT a LedgerAudit slot: an escape
// deposits nothing into rad_force/lepton_source
// (specs/hydro-coupling-source-terms.md:115 — "none (escape tallies, not
// source terms)"), so it can never join the grid-vs-event cancellation of
// ledger_closure; folding escapes into the emission/interaction step audits
// or MCNuX_LedgerClosure would break closure by construction.
//
// Tally semantics, per species s (species-major slot layout below):
//   count  += 1        one discrete escaped packet;
//   number += N        the packet weight (PIdx::w);
//   energy += N * p^t  with p^t the UPPER-index null-closure energy
//                      (p_t_closure of mcnux_tetrad.hxx) — the project-wide
//                      ledger energy convention
//                      (specs/hydro-coupling-source-terms.md:106-116,
//                      mcnux_srcterms.hxx), NOT the covariant pk_pt
//                      diagnostic column.
// Energy anchor convention: an escaped packet's position is outside the
// domain, where the metric gather may not be evaluated (the
// VertexMetricGather bounds assert), so p^t is evaluated with the FINAL
// momentum and the metric gathered at the LAST IN-DOMAIN position — the
// pre-step position in the push kernel, the last in-domain episode-start
// position in the episode driver. Exact on Minkowski (p^t = |p|
// everywhere, the escape-freestream benchmark).
//
// Escape predicate: half-open per axis — escaped iff any d has
// x[d] < prob_lo[d] || x[d] >= prob_hi[d] — matching the half-open
// cell-anchor convention of the fluid/metric gathers (a packet bitwise ON
// the upper domain face has left the last cell). Geometry::
// outsideRoundoffDomain is host-only and cannot be used in kernels, hence
// the explicit device-safe form here.
//
// Known deferred limitation (recorded per the task decision): the run
// accumulator is CUMULATIVE across the run (per-step +=, never reset,
// unlike the recompute-fresh diag precedents), but mcnux_escape_diag is
// checkpoint="no" like every mcnux_*_diag table and the accumulator itself
// is not checkpointed either — a restarted run's tally restarts from zero.
// No restart benchmark exists yet; restart persistence extends this header
// explicitly when one does.
//
// Deliberately free of cctk.h, AMReX, and CarpetX includes: plain portable
// constexpr C++ (the mcnux_srcterms.hxx pattern), so host code, device
// kernels, and the self-test battery all include it.

#include "mcnux_units.hxx"

namespace MCNuX {

// Slot layout of the tally and of MCNuX::mcnux_escape_diag: species-major,
// NUM_SPECIES x {count, number, energy}.
inline constexpr int escape_num_channels = 3; // count, number, energy
inline constexpr int escape_channel_count = 0;
inline constexpr int escape_channel_number = 1;
inline constexpr int escape_channel_energy = 2;
inline constexpr int escape_num_slots = NUM_SPECIES * escape_num_channels;

// The slot index map. Species indices come from species_index/NUM_SPECIES
// of mcnux_units.hxx ([MCNX-CNV-06]) — never re-enumerated here.
constexpr int escape_slot(int species_idx, int channel) noexcept {
  return species_idx * escape_num_channels + channel;
}

constexpr int escape_slot(Species s, int channel) noexcept {
  return escape_slot(species_index(s), channel);
}

// The half-open escape predicate documented above. Device-safe (plain
// arithmetic on by-value-captured bounds).
constexpr bool outside_domain(const double x[3], const double prob_lo[3],
                              const double prob_hi[3]) noexcept {
  for (int d = 0; d < 3; ++d)
    if (x[d] < prob_lo[d] || x[d] >= prob_hi[d])
      return true;
  return false;
}

// The cumulative per-species escape tally: escape_num_slots doubles in the
// species-major slot layout above.
struct EscapeTally {
  double v[escape_num_slots] = {};

  // One escaped packet of species index `species_idx`, weight N, and
  // upper-index null-closure energy pt (the semantics documented above).
  constexpr void accumulate(int species_idx, double N, double pt) noexcept {
    v[escape_slot(species_idx, escape_channel_count)] += 1.0;
    v[escape_slot(species_idx, escape_channel_number)] += N;
    v[escape_slot(species_idx, escape_channel_energy)] += N * pt;
  }

  // Fold a per-step slot buffer (the sites' device accumulators, copied to
  // the host) into the running totals.
  constexpr void add_slots(const double slots[escape_num_slots]) noexcept {
    for (int i = 0; i < escape_num_slots; ++i)
      v[i] += slots[i];
  }
};

// The run-cumulative accumulator behind MCNuX_EscapeDiag; defined in
// mcnux_escape.cxx (TU-local storage, the interaction_step_audit()
// mechanism — but NOT part of the mcnux_srcterms.hxx ledger family, per the
// storage-surface decision above).
EscapeTally &escape_run_tally();

// ---------------------------------------------------------------------------
// Compile-time verification (exact binary-fraction fixtures throughout)
// ---------------------------------------------------------------------------

static_assert(escape_num_slots == 9,
              "the escape tally is 3 species x {count, number, energy}");
static_assert(escape_slot(Species::NuE, escape_channel_count) == 0 &&
                  escape_slot(Species::NuE, escape_channel_energy) == 2 &&
                  escape_slot(Species::NuEBar, escape_channel_count) == 3 &&
                  escape_slot(Species::NuX, escape_channel_number) == 7 &&
                  escape_slot(Species::NuX, escape_channel_energy) == 8,
              "species-major escape slot map: slot = 3 * species + channel");

namespace detail {

// Predicate fixtures on the [-1, 1]^3 benchmark domain: interior inside;
// exactly on the LOWER face inside; exactly on the UPPER face escaped (the
// half-open convention); beyond a single axis escaped.
constexpr bool escape_predicate_holds() noexcept {
  constexpr double lo[3] = {-1.0, -1.0, -1.0};
  constexpr double hi[3] = {1.0, 1.0, 1.0};
  constexpr double interior[3] = {0.125, -0.25, 0.5};
  constexpr double on_lo[3] = {-1.0, 0.0, 0.0};
  constexpr double on_hi[3] = {1.0, 0.0, 0.0};
  constexpr double one_axis[3] = {0.0, 1.25, 0.0};
  return !outside_domain(interior, lo, hi) && !outside_domain(on_lo, lo, hi) &&
         outside_domain(on_hi, lo, hi) && outside_domain(one_axis, lo, hi);
}

// Tally arithmetic on a two-species hand fixture: two nu_e escapes
// (N, p^t) = (1.5, 2.0) and (0.5, 4.0), one nu_x escape (2.0, 0.25).
constexpr EscapeTally escape_tally_fixture() noexcept {
  EscapeTally t{};
  t.accumulate(species_index(Species::NuE), 1.5, 2.0);
  t.accumulate(species_index(Species::NuE), 0.5, 4.0);
  t.accumulate(species_index(Species::NuX), 2.0, 0.25);
  return t;
}

constexpr bool escape_tally_holds() noexcept {
  constexpr EscapeTally t = escape_tally_fixture();
  return t.v[0] == 2.0 && t.v[1] == 2.0 && t.v[2] == 5.0 && // nu_e
         t.v[3] == 0.0 && t.v[4] == 0.0 && t.v[5] == 0.0 && // nu_e_bar
         t.v[6] == 1.0 && t.v[7] == 2.0 && t.v[8] == 0.5;   // nu_x
}

} // namespace detail

static_assert(detail::escape_predicate_holds(),
              "the escape predicate must be half-open: lower face inside, "
              "upper face escaped ([MCNX-GPU-05])");
static_assert(detail::escape_tally_holds(),
              "EscapeTally::accumulate must add {1, N, N p^t} into exactly "
              "the species' three slots");

} // namespace MCNuX

#endif // MCNUX_ESCAPE_HXX
