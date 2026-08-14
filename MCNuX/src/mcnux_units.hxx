#ifndef MCNUX_UNITS_HXX
#define MCNUX_UNITS_HXX

// Constants, unit conversions, and the neutrino species enumeration.
//
// Governing spec: specs/conventions-and-units.md, requirements
// [MCNX-CNV-02..06]. This header is the *single* place in MCNuX where a
// physical constant literal may be written (consistency sweep,
// conventions-and-units.md:132): every conversion elsewhere in the thorn must
// be expressed in terms of the factors defined below.
//
// [MCNX-CNV-02]: the two unit systems are fixed at compile time. There is no
// runtime parameter — and there must never be a param.ccl entry — that selects
// units or constants; golden data must not be configuration-dependent.
//
// Deliberately free of cctk.h, AMReX, and CarpetX includes: plain portable
// C++ so that host code, device code, and standalone unit tests can all
// include it.

namespace MCNuX {

namespace detail {

// Constexpr |x|; <cmath>'s std::abs is not usable in a constant expression on
// every supported toolchain.
constexpr double cabs(double x) noexcept { return x < 0.0 ? -x : x; }

// Mixed relative/absolute closeness predicate of
// conventions-and-units.md:129, |a - b| <= rtol*|b| + 1e-30. The absolute
// floor keeps comparisons near zero well defined.
constexpr bool approx_eq(double a, double b, double rtol) noexcept {
  return cabs(a - b) <= rtol * cabs(b) + 1e-30;
}

// Tolerance tier for agreement with the spec's *pinned* 10-significant-digit
// table values. It cannot be the machine tier: a value truncated to 10
// significant digits differs from the full-precision one by up to ~5e-10
// relative — the spec itself exhibits this: its MeV -> K derivation
// (conventions-and-units.md:102) carries 17 digits and then states that the
// table value is that number "pinned to 10 digits". [MCNX-CNV-04]
// (conventions-and-units.md:123,141) explicitly permits carrying more digits
// internally, so the pinned values are the *agreement anchor* at their own
// precision: one unit in the 10th significant digit.
constexpr double rtol_pinned = 1e-9;

// The true machine tier of conventions-and-units.md:123,129-130, applied to
// identities and round trips among full-precision computed values.
constexpr double rtol_machine = 1e-14;

// cgs -> code -> cgs (and back) is the identity to the machine tier for every
// conversion factor (conventions-and-units.md:130). Probed over several
// decades so that a factor that "round trips" only near 1 cannot pass.
constexpr bool round_trips(double factor) noexcept {
  const double xs[] = {1.0, 1.2345678901234567e-7, 9.876543210987654e13,
                       3.14159265358979e5, 2.718281828459045e-3};
  for (const double x : xs) {
    if (!approx_eq(x / factor * factor, x, rtol_machine))
      return false;
    if (!approx_eq(x * factor / factor, x, rtol_machine))
      return false;
  }
  return true;
}

} // namespace detail

// ---------------------------------------------------------------------------
// Pinned defining constants  [MCNX-CNV-03]
// ---------------------------------------------------------------------------
// Exactly these five values are the inputs of every derivation below
// (conventions-and-units.md:68-74). No other constant set carried by the
// reference stack (HydroBaseX's README, CarpetX's openPMD metadata) is used
// anywhere in MCNuX; the discrepancy is recorded and resolved by this pinning
// (conventions-and-units.md:146).

// Speed of light [cm/s], SI exact.
inline constexpr double c_cgs = 2.99792458e10;

// Boltzmann constant [erg/K], SI-2019 exact.
inline constexpr double k_B_cgs = 1.380649e-16;

// Elementary charge e = 1.602176634e-19 C (SI-2019 exact), used throughout in
// its equivalent energy form: 1 MeV = 1.602176634e-6 erg, exact.
inline constexpr double MeV_to_erg = 1.602176634e-6;

// Gravitational constant [cm^3 g^-1 s^-2], CODATA 2018.
inline constexpr double G_cgs = 6.67430e-8;

// Nominal solar mass parameter [cm^3 s^-2], IAU 2015 B3, exact by convention.
inline constexpr double GM_sun_cgs = 1.3271244e26;

// ---------------------------------------------------------------------------
// Derived conversion factors  [MCNX-CNV-04]
// ---------------------------------------------------------------------------
// Geometrized units G = c = M_sun = 1 for spacetime and transport. Each factor
// below is *computed* from the five defining constants, never transcribed from
// the spec's table; the table values appear only as static_assert anchors at
// the bottom of this header. Names are directional: `x_code_to_cgs` multiplies
// a value in code units to give the cgs value.
//   value[cgs] = value[code] * x_code_to_cgs
//   value[code] = value[cgs] / x_code_to_cgs

// L = GM_sun/c^2 [cm].
inline constexpr double length_code_to_cgs = GM_sun_cgs / (c_cgs * c_cgs);

// t = L/c = GM_sun/c^3 [s].
inline constexpr double time_code_to_cgs =
    GM_sun_cgs / (c_cgs * c_cgs * c_cgs);

// M_sun = GM_sun/G [g].
inline constexpr double mass_code_to_cgs = GM_sun_cgs / G_cgs;

// M_sun c^2 = GM_sun c^2/G [erg].
inline constexpr double energy_code_to_cgs =
    GM_sun_cgs * c_cgs * c_cgs / G_cgs;

// The same energy unit expressed in the microphysics energy unit [MeV].
inline constexpr double energy_code_to_MeV = energy_code_to_cgs / MeV_to_erg;

// M_sun/L^3 [g/cm^3].
inline constexpr double mass_density_code_to_cgs =
    mass_code_to_cgs /
    (length_code_to_cgs * length_code_to_cgs * length_code_to_cgs);

// M_sun c^2/L^3 [erg/cm^3]. Energy density and pressure share this factor.
inline constexpr double energy_density_code_to_cgs =
    energy_code_to_cgs /
    (length_code_to_cgs * length_code_to_cgs * length_code_to_cgs);
inline constexpr double pressure_code_to_cgs = energy_density_code_to_cgs;

// Opacity (inverse length), 1/L [cm^-1]. Direction matters and is easy to get
// backwards, because an opacity transforms inversely to a length:
//   kappa[cm^-1] = kappa[code] * opacity_code_to_cgs
//   kappa[code]  = kappa[cm^-1] * length_code_to_cgs
// The second line is the form written in conventions-and-units.md:88.
inline constexpr double opacity_code_to_cgs = 1.0 / length_code_to_cgs;

// Rate (inverse time), 1/t [s^-1]; same inverse-direction caveat as opacity.
inline constexpr double rate_code_to_cgs = 1.0 / time_code_to_cgs;

// Emissivity, M_sun c^2/(L^3 t) [erg cm^-3 s^-1].
inline constexpr double emissivity_code_to_cgs =
    energy_density_code_to_cgs / time_code_to_cgs;

// ---------------------------------------------------------------------------
// Temperature at the WeakLibInterp boundary  [MCNX-CNV-05]
// ---------------------------------------------------------------------------
// Temperature is held in MeV everywhere inside MCNuX. Every temperature
// argument of every WeakLibInterp entry point — EOS, inversion, EmAb/Iso,
// NES/Pair, Brem — receives Kelvin (conventions-and-units.md:104,124). The
// NES/Pair axis is Kelvin too; assuming MeV there is a documented, easy-to-make
// error. This header only defines the factors; applying them at the call
// boundary is owned by the opacity/EOS evaluation code.

// T[K] = T[MeV] * mev_to_kelvin, exact under the SI-2019 k_B and e.
inline constexpr double mev_to_kelvin = MeV_to_erg / k_B_cgs;

// Boltzmann constant in microphysics units [MeV/K], the exact reciprocal of
// the factor above, for Boltzmann factors written with T in Kelvin.
inline constexpr double k_B_MeV_per_K = 1.0 / mev_to_kelvin;

// ---------------------------------------------------------------------------
// Neutrino species enumeration  [MCNX-CNV-06]
// ---------------------------------------------------------------------------
// Exactly three effective species in this binding order
// (conventions-and-units.md:110-114). Index 2 lumps nu_mu, nu_mu_bar, nu_tau,
// nu_tau_bar into one effective species (precedent arXiv:1708.08452).
enum class Species : int {
  NuE = 0,    // electron neutrino
  NuEBar = 1, // electron antineutrino
  NuX = 2     // heavy-lepton species, lumped
};

inline constexpr int NUM_SPECIES = 3;

// Degeneracy g, indexed by species index. Integer-valued and exact.
// NuX's g = 4 is *defined* here but must enter exactly once, at packet
// creation (emission sampling weights) — never in per-species transport
// coefficients (conventions-and-units.md:125).
inline constexpr int species_degeneracy[NUM_SPECIES] = {1, 1, 4};

// Lepton number l, indexed by species index. NuX carries zero net lepton
// number.
inline constexpr int species_lepton_number[NUM_SPECIES] = {1, -1, 0};

constexpr int species_index(Species s) noexcept { return static_cast<int>(s); }

constexpr int degeneracy(Species s) noexcept {
  return species_degeneracy[species_index(s)];
}

constexpr int lepton_number(Species s) noexcept {
  return species_lepton_number[species_index(s)];
}

// ---------------------------------------------------------------------------
// Compile-time verification (conventions-and-units.md:129-131)
// ---------------------------------------------------------------------------
// These static_asserts *are* the T1 test: any translation unit that includes
// this header re-verifies the whole constant set at build time.

// (a) Factor recomputation against the pinned 10-significant-digit table
//     values. Each pinned literal of conventions-and-units.md:82-90,99,103
//     appears in the repository exactly once — here. See detail::rtol_pinned
//     for why the tier is 1e-9 rather than 1e-14.

static_assert(detail::approx_eq(length_code_to_cgs, 1.476625038e5,
                                detail::rtol_pinned),
              "length code->cgs factor disagrees with the pinned value");
static_assert(detail::approx_eq(time_code_to_cgs, 4.925490948e-6,
                                detail::rtol_pinned),
              "time code->cgs factor disagrees with the pinned value");
static_assert(detail::approx_eq(mass_code_to_cgs, 1.988409871e33,
                                detail::rtol_pinned),
              "mass code->cgs factor disagrees with the pinned value");
static_assert(detail::approx_eq(energy_code_to_cgs, 1.787093669e54,
                                detail::rtol_pinned),
              "energy code->cgs factor disagrees with the pinned value");
static_assert(detail::approx_eq(energy_code_to_MeV, 1.115416135e60,
                                detail::rtol_pinned),
              "energy code->MeV factor disagrees with the pinned value");
static_assert(detail::approx_eq(mass_density_code_to_cgs, 6.175828479e17,
                                detail::rtol_pinned),
              "mass density code->cgs factor disagrees with the pinned value");
static_assert(detail::approx_eq(energy_density_code_to_cgs, 5.550557829e38,
                                detail::rtol_pinned),
              "energy density code->cgs factor disagrees with the pinned "
              "value");
static_assert(detail::approx_eq(opacity_code_to_cgs, 6.772199944e-6,
                                detail::rtol_pinned),
              "opacity code->cgs factor disagrees with the pinned value");
static_assert(detail::approx_eq(rate_code_to_cgs, 2.030254467e5,
                                detail::rtol_pinned),
              "rate code->cgs factor disagrees with the pinned value");
static_assert(detail::approx_eq(emissivity_code_to_cgs, 1.126904483e44,
                                detail::rtol_pinned),
              "emissivity code->cgs factor disagrees with the pinned value");
static_assert(detail::approx_eq(mev_to_kelvin, 1.160451812e10,
                                detail::rtol_pinned),
              "MeV->Kelvin factor disagrees with the pinned value");
static_assert(detail::approx_eq(k_B_MeV_per_K, 8.617333262e-11,
                                detail::rtol_pinned),
              "k_B in MeV/K disagrees with the pinned value");

// (b) Identities and round trips among full-precision computed values, at the
//     machine tier of conventions-and-units.md:123,129-130.

static_assert(detail::approx_eq(mev_to_kelvin * k_B_MeV_per_K, 1.0,
                                detail::rtol_machine),
              "MeV->K factor and k_B[MeV/K] are not mutual reciprocals");
static_assert(detail::approx_eq(opacity_code_to_cgs * length_code_to_cgs, 1.0,
                                detail::rtol_machine),
              "opacity factor is not the inverse of the length factor");
static_assert(detail::approx_eq(rate_code_to_cgs * time_code_to_cgs, 1.0,
                                detail::rtol_machine),
              "rate factor is not the inverse of the time factor");
static_assert(detail::approx_eq(time_code_to_cgs * c_cgs, length_code_to_cgs,
                                detail::rtol_machine),
              "t = L/c is violated by the computed factors");
static_assert(detail::approx_eq(mass_code_to_cgs * c_cgs * c_cgs,
                                energy_code_to_cgs, detail::rtol_machine),
              "E = M c^2 is violated by the computed factors");
static_assert(detail::approx_eq(energy_code_to_MeV * MeV_to_erg,
                                energy_code_to_cgs, detail::rtol_machine),
              "the MeV and erg forms of the energy unit disagree");
static_assert(detail::approx_eq(emissivity_code_to_cgs * time_code_to_cgs,
                                energy_density_code_to_cgs,
                                detail::rtol_machine),
              "emissivity is not energy density per unit time");
static_assert(detail::approx_eq(mass_density_code_to_cgs * c_cgs * c_cgs,
                                energy_density_code_to_cgs,
                                detail::rtol_machine),
              "energy density is not c^2 times mass density");

static_assert(detail::round_trips(length_code_to_cgs), "length round trip");
static_assert(detail::round_trips(time_code_to_cgs), "time round trip");
static_assert(detail::round_trips(mass_code_to_cgs), "mass round trip");
static_assert(detail::round_trips(energy_code_to_cgs), "energy round trip");
static_assert(detail::round_trips(energy_code_to_MeV), "energy/MeV round trip");
static_assert(detail::round_trips(mass_density_code_to_cgs),
              "mass density round trip");
static_assert(detail::round_trips(energy_density_code_to_cgs),
              "energy density round trip");
static_assert(detail::round_trips(opacity_code_to_cgs), "opacity round trip");
static_assert(detail::round_trips(rate_code_to_cgs), "rate round trip");
static_assert(detail::round_trips(emissivity_code_to_cgs),
              "emissivity round trip");
static_assert(detail::round_trips(mev_to_kelvin), "T[MeV]->T[K] round trip");
static_assert(detail::round_trips(k_B_MeV_per_K), "T[K]->T[MeV] round trip");
static_assert(detail::approx_eq(1.5 * mev_to_kelvin * k_B_MeV_per_K, 1.5,
                                detail::rtol_machine),
              "T[MeV] -> T[K] -> T[MeV] is not the identity");
static_assert(detail::approx_eq(pressure_code_to_cgs,
                                energy_density_code_to_cgs, 0.0),
              "pressure and energy density must share one factor");

// (c) Species table — exact integer equality, no tolerance.

static_assert(NUM_SPECIES == 3, "the species set is exactly {nu_e, nu_e_bar, "
                                "nu_x}");
static_assert(species_index(Species::NuE) == 0, "nu_e is index 0");
static_assert(species_index(Species::NuEBar) == 1, "nu_e_bar is index 1");
static_assert(species_index(Species::NuX) == 2, "nu_x is index 2");

static_assert(species_degeneracy[0] == 1, "g(nu_e) = 1");
static_assert(species_degeneracy[1] == 1, "g(nu_e_bar) = 1");
static_assert(species_degeneracy[2] == 4, "g(nu_x) = 4");
static_assert(degeneracy(Species::NuE) == 1, "g(nu_e) = 1");
static_assert(degeneracy(Species::NuEBar) == 1, "g(nu_e_bar) = 1");
static_assert(degeneracy(Species::NuX) == 4, "g(nu_x) = 4");

static_assert(species_lepton_number[0] == +1, "l(nu_e) = +1");
static_assert(species_lepton_number[1] == -1, "l(nu_e_bar) = -1");
static_assert(species_lepton_number[2] == 0, "l(nu_x) = 0");
static_assert(lepton_number(Species::NuE) == +1, "l(nu_e) = +1");
static_assert(lepton_number(Species::NuEBar) == -1, "l(nu_e_bar) = -1");
static_assert(lepton_number(Species::NuX) == 0, "l(nu_x) = 0");
static_assert(species_lepton_number[0] + species_lepton_number[1] == 0,
              "nu_e and nu_e_bar carry opposite lepton number");

} // namespace MCNuX

#endif // MCNUX_UNITS_HXX
