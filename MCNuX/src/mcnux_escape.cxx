// Runtime side of the [MCNX-GPU-05] escape tallies
// (specs/particle-container-and-gpu.md; the escape clause: packets that
// leave the outer domain boundary are removed and tallied, never silently
// dropped): the run-cumulative accumulator behind escape_run_tally() and
// the gated diagnostic writer MCNuX_EscapeDiag.
//
// The detection/tally/removal itself lives at the two live escape sites —
// the geodesic push (mcnux_geodesic.cxx) and the episode driver
// (mcnux_interactions.cxx) — which fold their per-step device buffers into
// the accumulator here after their Redistribute(). Detection is ALWAYS ON
// (binding correctness); only the diagnostic table below is gated
// (MCNuX::test_escape_tallies, the `escape-freestream` benchmark).
//
// The accumulator is CUMULATIVE across the run (per-step +=, never reset —
// unlike the recompute-fresh diag precedents); its restart persistence is a
// recorded deferred limitation (see mcnux_escape.hxx).

#include "mcnux_escape.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Groups.h>
#include <cctk_Parameters.h>

namespace MCNuX {

namespace {

// The run-cumulative tally (TU-local storage; the interaction_step_audit()
// mechanism, deliberately outside the mcnux_srcterms.hxx ledger family —
// escapes never join the closure audit, per the mcnux_escape.hxx decision).
EscapeTally g_escape_tally{};

} // namespace

EscapeTally &escape_run_tally() { return g_escape_tally; }

// Mirror the running totals into MCNuX::mcnux_escape_diag (species-major,
// 3 species x {count, number Sum N, energy Sum N p^t}), cross-checking the
// declared SIZE against escape_num_slots (the SIZE-agreement abort idiom of
// mcnux_geodesic.cxx / mcnux_selftest.cxx).
extern "C" void MCNuX_EscapeDiag(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_EscapeDiag;
  DECLARE_CCTK_PARAMETERS;

  const int gi = CCTK_GroupIndex("MCNuX::mcnux_escape_diag");
  if (gi < 0)
    CCTK_VERROR("MCNuX::mcnux_escape_diag is not a known group");
  const CCTK_INT *const *const sizes = CCTK_GroupSizesI(gi);
  if (!sizes || CCTK_GroupDimI(gi) != 1)
    CCTK_VERROR("MCNuX::mcnux_escape_diag must be a 1-dimensional grid array");
  const int nrows = int(*sizes[0]);
  if (nrows != escape_num_slots)
    CCTK_VERROR("MCNuX::mcnux_escape_diag has SIZE=%d but the escape tally "
                "carries %d slots (mcnux_escape.hxx); the two must agree "
                "(MCNuX/interface.ccl)",
                nrows, escape_num_slots);

  const EscapeTally &t = escape_run_tally();
  for (int i = 0; i < escape_num_slots; ++i)
    esc_tally[i] = t.v[i];

  for (int s = 0; s < NUM_SPECIES; ++s)
    CCTK_VINFO("MCNuX escape tally species %d: count = %.17g, "
               "number = %.17g, energy = %.17g",
               s, t.v[escape_slot(s, escape_channel_count)],
               t.v[escape_slot(s, escape_channel_number)],
               t.v[escape_slot(s, escape_channel_energy)]);
}

} // namespace MCNuX
