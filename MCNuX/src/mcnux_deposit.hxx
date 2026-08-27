#ifndef MCNUX_DEPOSIT_HXX
#define MCNUX_DEPOSIT_HXX

// The shared atomic source-term deposition helper (the C2 helper): the ONE
// place a packet kernel turns a pre-negation LedgerDelta into the five
// atomic += adds on the MCNuX-owned source-term grid variables
// {rf_t, rf_x, rf_y, rf_z, lep_src} of MCNuX/interface.ccl.
//
// Governing specs:
//   * [MCNX-HYD-02] (specs/hydro-coupling-source-terms.md) — the ledger
//     arithmetic itself is source_from_delta of mcnux_srcterms.hxx, applied
//     verbatim here (code units, per unit coordinate volume per unit
//     coordinate time, undensitized).
//   * [MCNX-GPU-02]/[MCNX-RNG-08] (specs/particle-container-and-gpu.md:
//     156-171) — many packets of one kernel may deposit into the same cell,
//     so every add is an amrex::Gpu::Atomic::AddNoRet on the cell's value.
//     This helper is the packet-deposition idiom; the grid-point-loop `+=`
//     of MCNuX_SyntheticDeposit (one contribution per cell inside a grid
//     loop) is NOT it and must not be copied for packet events.
//
// Consumers: the emission loop (mcnux_emission.cxx) now; the episode driver
// and absorption/removal operators later — all reuse this ONE helper.

#include "mcnux_srcterms.hxx"

#include <AMReX_Array4.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cctk.h>

#include <type_traits>

namespace MCNuX {

// The two MCNuX source-term groups (group names per MCNuX/interface.ccl:
// MCNuX::rad_force = {rf_t, rf_x, rf_y, rf_z}, MCNuX::lepton_source =
// {lep_src}). Mirrors metric_groups() of mcnux_gather.hxx.
struct SourceGroups {
  int rad_force, lepton_source;
};

inline SourceGroups source_groups() {
  SourceGroups g{CCTK_GroupIndex("MCNuX::rad_force"),
                 CCTK_GroupIndex("MCNuX::lepton_source")};
  if (g.rad_force < 0 || g.lepton_source < 0)
    CCTK_VERROR("MCNuX deposition needs the MCNuX::rad_force and "
                "MCNuX::lepton_source groups "
                "(specs/hydro-coupling-source-terms.md [MCNX-HYD-01])");
  return g;
}

// The by-value view pack a deposit kernel captures: one single-component
// Array4 per source variable, laid out to match {rf_t, rf_x, rf_y, rf_z,
// lep_src} (the SourceContribution order of mcnux_srcterms.hxx). Built per
// box from the group MultiFabs via array(box, comp) — rad_force components
// 0..3, lepton_source component 0.
struct SourceViews {
  amrex::Array4<amrex::Real> rf_t, rf_x, rf_y, rf_z, lep_src;
};

static_assert(std::is_trivially_copyable_v<SourceViews>,
              "[MCNX-GPU-02] SourceViews is captured by value into packet "
              "kernels and must be trivially copyable");

// Deposit one event's LedgerDelta into cell (i, j, k): the [MCNX-HYD-02]
// normalization source = -delta/(dV dt) via source_from_delta (never
// re-derived here), then five atomic adds (the [MCNX-GPU-02] deposition
// idiom — deterministic up to floating-point addition order; bitwise
// reproducibility of deposits is governed by [MCNX-RNG-08], a later task's
// concern, not this helper's). dV and dt are the COORDINATE (code-unit)
// cell volume and transport step; the caller's event-side ledger audit must
// use the identical pair ([MCNX-HYD-05] pairing).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
deposit_delta(const SourceViews &v, int i, int j, int k, const LedgerDelta &d,
              double dV, double dt) noexcept {
  const SourceContribution s = source_from_delta(d, dV, dt);
  amrex::Gpu::Atomic::AddNoRet(&v.rf_t(i, j, k), amrex::Real(s.Gt));
  amrex::Gpu::Atomic::AddNoRet(&v.rf_x(i, j, k), amrex::Real(s.Gx));
  amrex::Gpu::Atomic::AddNoRet(&v.rf_y(i, j, k), amrex::Real(s.Gy));
  amrex::Gpu::Atomic::AddNoRet(&v.rf_z(i, j, k), amrex::Real(s.Gz));
  amrex::Gpu::Atomic::AddNoRet(&v.lep_src(i, j, k), amrex::Real(s.Sl));
}

} // namespace MCNuX

#endif // MCNUX_DEPOSIT_HXX
