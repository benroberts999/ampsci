#include "Modules/StructRad.hpp"
#include "DiracOperator/include.hpp" //For E1 operator
#include "IO/InputBlock.hpp"
#include "MBPT/Feynman.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

void StructRad(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check({
    {"stride", "stride [8]"},
    {"n_min_core", "Minimum core state [5]"},
    {"rmax", "Maximum r value [30]"},
    {"r0", "Minimum r value, [1e-3]"},
    {"includeG", "include lower g component [false]"},
    {"screening", "Coulomb screening [false]"},
  });

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }
  const auto hE1 = DiracOperator::E1(wf.grid());
  auto r0 = input.get("r0", 1e-3);
  auto rmax = input.get("rmax", 30.0);
  const auto i0 = wf.grid().getIndex(r0);
  const auto imax = wf.grid().getIndex(rmax);
  const std::size_t stride = input.get("stride", 8);
  const auto n_min_core = input.get("n_min_core", 5);
  const bool includeG = input.get("includeG", false);
  const auto gsize = (imax - i0) / stride;
  const auto screening = input.get("screening", false);
  const auto options =
    MBPT::FeynmanOptions{.screening = MBPT::Screening::include,
                         .hole_particle = MBPT::HoleParticle::include};

  const auto Fyn = (screening) ?
                     MBPT::Feynman(wf.vHF(), i0, stride, gsize, options,
                                   n_min_core, includeG, true, "tempqpq1") :
                     MBPT::Feynman(wf.vHF(), i0, stride, gsize, {}, n_min_core,
                                   includeG, true, "tempqpq1");

  std::cout << "Grid size: " << gsize << ", "
            << "Stride: " << stride << "\n";
  std::cout << "Coulomb Screening: " << screening << ", Include G:" << includeG
            << "\n";
  // 2) Loop through each pair of valence states, calc E1 matrix elements:
  for (auto a = 0ul; a < wf.valence().size(); ++a) {
    const auto &Fa = wf.valence()[a];
    for (auto b = 0ul; b < wf.valence().size(); ++b) {
      const auto &Fb = wf.valence()[b];

      // Skip the MEs which are zero due to selection rules:
      if (hE1.isZero(Fa.kappa(), Fb.kappa()))
        continue;

      const auto Omega = Fa.en() - Fb.en();

      const auto Sd =
        Fyn.Sigma_SR_direct(Fa.kappa(), Fb.en(), Fb.kappa(), Omega);
      const auto e_AB = (Fa) * (Sd * Fb);

      std::cout << "Re <" << Fa.shortSymbol() << "|d|" << Fb.shortSymbol()
                << "> = " << e_AB << "\n";
    }
  }
}

} // namespace Module
