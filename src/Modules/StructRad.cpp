#include "Modules/StructRad.hpp"
#include "DiracOperator/include.hpp" //For E1 operator
#include "IO/InputBlock.hpp"
#include "MBPT/FeynmanSR.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Modules/Modules.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

// Declare, register, then define below.
namespace {
const Register r_StructRad{"StructRad", "description", &StructRad};
} // namespace

void StructRad(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check({
    {"stride", "stride [8]"},
    {"n_min_core", "Minimum core state [5]"},
    {"rmax", "Maximum r value [30]"},
    {"r0", "Minimum r value, [1e-3]"},
    {"includeG", "include lower g component [false]"},
    {"includeSR2", "include nd diagrams [false]"},
    {"screening", "Coulomb screening [false]"},
  });

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }
  const auto hE1 = DiracOperator::E1(wf.grid());
  auto r0 = input.get("r0", 1e-3);
  auto rmax = input.get("rmax", 90.0);
  const auto i0 = wf.grid().getIndex(r0);
  const auto imax = wf.grid().getIndex(rmax);
  const std::size_t stride = input.get("stride", 8ul);
  const auto n_min_core = input.get("n_min_core", 5);
  const bool includeG = input.get("includeG", false);
  const bool includeSR2 = input.get("includeSR2", false);
  const auto gsize = (imax - i0) / stride;
  const auto screening = input.get("screening", false);
  MBPT::FeynmanOptions options{};
  options.screening = MBPT::Screening::include;
  options.hole_particle = MBPT::HoleParticle::include;

  auto Fyn = (screening) ?
               MBPT::FeynmanSR(wf.vHF(), i0, stride, gsize, options, n_min_core,
                               includeG, true, "tempqpq1") :
               MBPT::FeynmanSR(wf.vHF(), i0, stride, gsize, {}, n_min_core,
                               includeG, true, "tempqpq1");

  std::cout << "Grid size: " << gsize << ", " << "Stride: " << stride
            << ", Rmax: " << rmax << "\n";
  std::cout << "Coulomb Screening: " << screening << ", Include G:" << includeG
            << "\n";
  std::cout << "Include SR2: " << includeSR2 << "\n";
  // 2) Loop through each pair of valence states, calc E1 matrix elements:
  for (auto a = 0ul; a < wf.valence().size(); ++a) {
    const auto &Fa = wf.valence()[a];
    for (auto b = 0ul; b < wf.valence().size(); ++b) {
      const auto &Fb = wf.valence()[b];

      // Skip the MEs which are zero due to selection rules:
      if (hE1.isZero(Fa.kappa(), Fb.kappa()))
        continue;

      const auto Omega = Fa.en() - Fb.en();

      const auto Sd1 =
        Fyn.Sigma_SR_direct(Fa.kappa(), Fb.en(), Fb.kappa(), Omega);

      const auto e_AB = (Fa) * (Sd1 * Fb);

      std::cout << "<" << Fa.shortSymbol() << "|d|" << Fb.shortSymbol()
                << "> = " << e_AB << "\n"
                << std::flush;

      if (includeSR2) {
        const auto Sd2 =
          Fyn.Sigma_SR_direct2(Fa.kappa(), Fb.en(), Fb.kappa(), Omega);

        const auto e_AB2 = (Fa) * (Sd2 * Fb);

        std::cout << "<" << Fa.shortSymbol() << "|d|" << Fb.shortSymbol()
                  << ">: a = " << e_AB << ", b + c = " << e_AB2 << "\n"
                  << std::flush;
      }
      // const auto Sd2 =
      //   Fyn.Sigma_SR_direct2(Fa.kappa(), Fb.en(), Fb.kappa(), Omega);

      // // const auto Sd = Sd1 + Sd2;
      // const auto e_AB2 = (Fa) * (Sd2 * Fb);

      // std::cout << "Re <" << Fa.shortSymbol() << "|d|" << Fb.shortSymbol()
      //           << "> = " << e_AB << "," << e_AB2 << "\n"
      //           << std::flush;
    }
  }
}

} // namespace Module
