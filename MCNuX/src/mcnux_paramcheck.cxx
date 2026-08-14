#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace MCNuX {

// Single-rate pin of specs/carpetx-thorn-integration.md [MCNX-CTX-06].
//
// MCNuX transport is operator-split and advances every packet on every
// refinement level exactly once per coarsest-level Delta t ([MCNX-CTX-05]).
// With CarpetX::use_subcycling = yes the refined levels advance in per-level
// batches with different Delta t, and neither the transport cadence nor the
// packet redistribution point is specified for that regime — it is an open
// question of the spec, and running unspecified is not permitted. The guard
// below keeps that region unreachable: it aborts at CCTK_PARAMCHECK, before
// any evolution can start.
//
// use_subcycling enters scope unqualified via DECLARE_CCTK_PARAMETERS from
// the SHARES: Driver / USES BOOLEAN block of MCNuX/param.ccl.

extern "C" void MCNuX_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_ParamCheck;
  DECLARE_CCTK_PARAMETERS;

  if (use_subcycling)
    CCTK_VERROR(
        "MCNuX requires CarpetX::use_subcycling = no, but this run sets it to "
        "yes. MCNuX transport is single-rate: it advances all packets on all "
        "refinement levels exactly once per coarsest-level dt, which is "
        "well-defined only when every level shares that timestep. Transport "
        "under subcycling is unspecified (see the open question of "
        "specs/carpetx-thorn-integration.md [MCNX-CTX-06]). Set "
        "CarpetX::use_subcycling = no, or deactivate the MCNuX thorn.");
}

} // namespace MCNuX
