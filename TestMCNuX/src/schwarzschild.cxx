// Isotropic-coordinate Schwarzschild initial data
// (specs/verification-suite-design.md [MCNX-VER-02/03]):
//
//   r = sqrt(x^2 + y^2 + z^2),   psi = 1 + M/(2 r)
//   gamma_ij = psi^4 delta_ij,   K_ij = 0
//   alpha = (1 - M/(2 r)) / (1 + M/(2 r))
//
// The shift is not written here: ADMBaseX's own "zero" initializer is exact
// for this static beta^i = 0 solution. The formulas are evaluated as written
// at every grid point; behavior at r <= M/2 is unspecified by the spec and
// unreachable in the pinned benchmark domains (all placed in r > M/2).
//
// Writer idiom follows CarpetX/ADMBaseX/src/admbase.cxx: vertex-centered
// GF3D2 wrappers and a grid.loop_all_device<0,0,0> device loop.

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>
#include <cmath>

namespace TestMCNuX {
using namespace Loop;

extern "C" void TestMCNuX_Schwarzschild_InitialData(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestMCNuX_Schwarzschild_InitialData;
  DECLARE_CCTK_PARAMETERS;

  const std::array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> gxx_(layout, gxx);
  const GF3D2<CCTK_REAL> gxy_(layout, gxy);
  const GF3D2<CCTK_REAL> gxz_(layout, gxz);
  const GF3D2<CCTK_REAL> gyy_(layout, gyy);
  const GF3D2<CCTK_REAL> gyz_(layout, gyz);
  const GF3D2<CCTK_REAL> gzz_(layout, gzz);

  const GF3D2<CCTK_REAL> kxx_(layout, kxx);
  const GF3D2<CCTK_REAL> kxy_(layout, kxy);
  const GF3D2<CCTK_REAL> kxz_(layout, kxz);
  const GF3D2<CCTK_REAL> kyy_(layout, kyy);
  const GF3D2<CCTK_REAL> kyz_(layout, kyz);
  const GF3D2<CCTK_REAL> kzz_(layout, kzz);

  const CCTK_REAL mass = M;

  const GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using std::sqrt;
        const CCTK_REAL r = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        const CCTK_REAL psi = 1 + mass / (2 * r);
        const CCTK_REAL psi2 = psi * psi;
        const CCTK_REAL psi4 = psi2 * psi2;

        gxx_(p.I) = psi4;
        gxy_(p.I) = 0;
        gxz_(p.I) = 0;
        gyy_(p.I) = psi4;
        gyz_(p.I) = 0;
        gzz_(p.I) = psi4;

        kxx_(p.I) = 0;
        kxy_(p.I) = 0;
        kxz_(p.I) = 0;
        kyy_(p.I) = 0;
        kyz_(p.I) = 0;
        kzz_(p.I) = 0;
      });
}

extern "C" void TestMCNuX_Schwarzschild_InitialLapse(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestMCNuX_Schwarzschild_InitialLapse;
  DECLARE_CCTK_PARAMETERS;

  const std::array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> alp_(layout, alp);

  const CCTK_REAL mass = M;

  const GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using std::sqrt;
        const CCTK_REAL r = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        const CCTK_REAL half_m_over_r = mass / (2 * r);
        alp_(p.I) = (1 - half_m_over_r) / (1 + half_m_over_r);
      });
}

} // namespace TestMCNuX
