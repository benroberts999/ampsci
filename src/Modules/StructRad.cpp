#include "Modules/StructRad.hpp"
#include "DiracOperator/include.hpp" //For E1 operator
#include "IO/InputBlock.hpp"
#include "MBPT/Feynman.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

void StructRad(const IO::InputBlock &input, const Wavefunction &wf) {
  // This is an example module, designed to help you write new modules

  // In this example, we will solve a new wavefunction, assuming a different
  // nuclear charge distribution, and see the difference in the energies and E1
  // matrix elements this produces.

  // Read in some optional input options: A, rms, and type
  // "checkBlock" is optional, but helps prevent errors on input options:

  // Use the same Grid and alpha, but different nuclear parameters (except Z)

  const auto hE1 = DiracOperator::E1(wf.grid());
  const auto i0 = wf.grid().getIndex(1e-3);
  const auto imax = wf.grid().getIndex(30.0);
  const auto gsize = 60;
  const std::size_t stride = (imax - i0) / gsize;
  const auto m_n_min_core = 5;
  const bool m_includeG = false;
  const auto Fyn = MBPT::Feynman(wf.vHF(), i0, stride, gsize, {}, m_n_min_core,
                                 m_includeG, true, "tempqpq");
  std::cout << "Grid size: " << gsize << ", " << "Stride: " << stride << "\n";

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
        Fyn.Sigma_SR_direct(Fa.kappa(), Fa.en(), Fb.kappa(), Omega);
      const auto e_AB = (Fa) * (Sd * Fb);

      std::cout << "Re <" << Fa.shortSymbol() << "|" << Fb.shortSymbol()
                << "> = " << e_AB << "\n";
    }
  }
}

} // namespace Module
