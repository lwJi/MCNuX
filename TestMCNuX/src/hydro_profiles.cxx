// Analytic hydro test profiles
// (specs/verification-suite-design.md [MCNX-VER-04]):
//
//   uniform sphere:      (rho, T, Ye)(x) = (rho0, T0, Ye0)          if |x - x0| <= R_s
//                                          (rho_amb, T_amb, Ye_amb) otherwise
//   equilibration box:   (rho, T, Ye)(x) = (rho0, T0, Ye0)          everywhere
//
// Both profiles set the fluid at rest (vel^i = 0) and zero every other
// HydroBaseX variable the writer touches (eps, press, entropy, Bvec, Avec*).
// Units pass through verbatim in MCNuX's microphysics convention: rho in
// g/cm^3, temperature in MeV ([MCNX-HYD-04]), Ye dimensionless; no unit
// conversion happens anywhere in this writer. Discontinuity convention at the
// sphere surface: the cell-centered value is the profile evaluated at the
// cell center -- no volume-fraction smoothing.
//
// Writer idiom follows CarpetX/HydroBaseX/src/hydrobase.cxx: the
// DECLARE_CCTK_ARGUMENTSX accessor style, one grid.loop_all_device<1,1,1>
// for the cell-centered (ccc) fields, then three separate staggered loops
// for the vector potentials (Avecx cvv, Avecy vcv, Avecz vvc -- a single
// <1,1,1> loop cannot write them).

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace TestMCNuX {
using namespace Loop;

extern "C" void TestMCNuX_HydroProfiles(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestMCNuX_HydroProfiles;
  DECLARE_CCTK_PARAMETERS;

  // Hoisted host-side branch: the equilibration box is the sphere's
  // everywhere-inside branch.
  const bool use_sphere = CCTK_EQUALS(initial_hydro, "uniform sphere");

  const CCTK_REAL rho_in = rho0;
  const CCTK_REAL T_in = T0;
  const CCTK_REAL Ye_in = Ye0;
  const CCTK_REAL rho_out = rho_amb;
  const CCTK_REAL T_out = T_amb;
  const CCTK_REAL Ye_out = Ye_amb;
  const CCTK_REAL x0 = sphere_x0;
  const CCTK_REAL y0 = sphere_y0;
  const CCTK_REAL z0 = sphere_z0;
  const CCTK_REAL radius = sphere_radius;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using std::sqrt;
        const CCTK_REAL dx = p.x - x0;
        const CCTK_REAL dy = p.y - y0;
        const CCTK_REAL dz = p.z - z0;
        const CCTK_REAL dist = sqrt(dx * dx + dy * dy + dz * dz);
        const bool inside = !use_sphere || dist <= radius;

        rho(p.I) = inside ? rho_in : rho_out;
        temperature(p.I) = inside ? T_in : T_out;
        Ye(p.I) = inside ? Ye_in : Ye_out;

        velx(p.I) = 0;
        vely(p.I) = 0;
        velz(p.I) = 0;
        eps(p.I) = 0;
        press(p.I) = 0;
        entropy(p.I) = 0;
        Bvecx(p.I) = 0;
        Bvecy(p.I) = 0;
        Bvecz(p.I) = 0;
      });

  grid.loop_all_device<1, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avecx(p.I) = 0; });
  grid.loop_all_device<0, 1, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avecy(p.I) = 0; });
  grid.loop_all_device<0, 0, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avecz(p.I) = 0; });
}

} // namespace TestMCNuX
