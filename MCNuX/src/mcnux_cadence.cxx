// Transport cadence of specs/carpetx-thorn-integration.md [MCNX-CTX-05], and
// the driver-native grid-access surface of [MCNX-CTX-01].
//
// This translation unit holds the two routines of the MCNuX transport-step
// group declared in MCNuX/schedule.ccl: the counter initializer and the
// cadence observer. The observer is deliberately minimal — it advances an
// MCNuX-owned grid scalar and reports it — but it already exercises the three
// obligations that travel with the driver-native access idiom:
//
//   1. [MCNX-CTX-02] the hard `REQUIRES CarpetX` of MCNuX/configuration.ccl,
//      which is what makes the driver-internal include below well-defined;
//   2. [MCNX-CTX-03] global-mode scheduling, so the routine runs exactly once
//      per CCTK_EVOL traversal rather than once per (patch, level, component);
//   3. [MCNX-CTX-04] complete READS:/WRITES: clauses for every grid variable
//      touched, enforced at runtime by the driver's poison/checksum/presync
//      machinery.

// TODO: Don't include files from other thorns; create a proper interface
//
// CarpetX exposes no public capability header for `ghext` (its
// `PROVIDES CarpetX { }` block is empty), so reaching the driver's data model
// uses the stack's own idiom: the TODO-flagged relative include of the
// driver's `driver.hxx` (CarpetX/ODESolvers/src/solve.hxx). MCNuX lives in its
// own arrangement, i.e. at arrangements/MCNuX/MCNuX/src/, one directory deeper
// than ODESolvers relative to the CarpetX arrangement — hence the extra `../`.
// This is an implementation choice, not a contract: [MCNX-CTX-01] binds the
// objects accessed, not the include mechanism.
#include "../../../CarpetX/CarpetX/src/driver.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstddef>

namespace MCNuX {

namespace {

// Structural summary of the driver's data model, read through the singleton
// `CarpetX::ghext` (driver.hxx): the per-patch AmrCore and per-level group
// data of the [MCNX-CTX-01] access table hang off these containers, and the
// packet container and gather/deposit kernels of the later specs are built on
// them. Nothing here dereferences grid data, so no READS:/WRITES: clause
// beyond the counter arises from this traversal.
struct DriverExtent {
  int num_patches;
  int max_num_levels;
};

DriverExtent driver_extent() {
  if (!CarpetX::ghext)
    CCTK_VERROR("MCNuX reached the transport step, but the CarpetX driver "
                "singleton `ghext` is null. MCNuX accesses grid data through "
                "the driver's native AMReX objects "
                "(specs/carpetx-thorn-integration.md [MCNX-CTX-01]) and cannot "
                "run without them.");

  DriverExtent extent{int(CarpetX::ghext->patchdata.size()), 0};
  for (const auto &patchdata : CarpetX::ghext->patchdata)
    if (int(patchdata.leveldata.size()) > extent.max_num_levels)
      extent.max_num_levels = int(patchdata.leveldata.size());
  return extent;
}

} // namespace

extern "C" void MCNuX_TransportInit(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_TransportInit;

  *transport_step_count = 0;
}

extern "C" void MCNuX_TransportObserver(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MCNuX_TransportObserver;

  const DriverExtent extent = driver_extent();

  *transport_step_count += 1;

  // The cadence observable of [MCNX-CTX-05]: exactly one such line per
  // coarsest-level Delta t, emitted after ODESolvers_Solve has completed the
  // fluid step. `cctk_delta_time` in global mode is that coarsest-level step.
  CCTK_VINFO("MCNuX transport step %d: iteration %d, t = %.17g, dt = %.17g, "
             "patches = %d, max levels = %d",
             int(*transport_step_count), cctk_iteration, double(cctk_time),
             double(cctk_delta_time), extent.num_patches,
             extent.max_num_levels);
}

} // namespace MCNuX
