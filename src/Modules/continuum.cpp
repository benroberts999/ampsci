#include "Modules/continuum.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

void continuum(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check({{"ec", "energy for cntm states (>0)"},
                     {"max_l", "maximum l"},
                     {"filename", "filename for output"},
                     {"force_rescale", "force rescale"},
                     {"force_orthog", "force orthogonality"},
                     {"n", "principal quantum number for bound state"},
                     {"kappa", "Dirac quantum number for bound state"}});

  auto en_c = input.get("ec", 0.5);
  auto lmax = input.get("max_l", 0);
  auto fname = input.get("filename", std::string{""});

  // Method options for solveContinuumHF
  auto force_rescale = input.get<bool>("force_rescale", false);
  // auto subtract_self = input.get<bool>("subtract_self", false);
  auto force_orthog = input.get<bool>("force_orthog", false);

  // n=0 invalid -> nullptr
  auto n = input.get<int>("n", 0);
  auto k = input.get<int>("kappa", 0);
  // auto p_psi = wf.getState(n, k);
  auto p_psi = (n == 0) ? nullptr : wf.getState(n, k);

  if (p_psi != nullptr) {
    std::cout << "\n"
              << "State for n = " << n << " and kappa = " << k << ": "
              << p_psi->symbol() << "\n";
    // std::cout << "Bound state: " << p_psi->symbol() << "\n";
  }

  std::cout << "\nContinuum Orbitals:\n";
  std::cout << "energy: " << en_c << "\n";

  auto cntm = ContinuumOrbitals(wf);
  cntm.solveContinuumHF(en_c, lmax, force_rescale, true, force_orthog, p_psi);

  std::cout << "Check orthogonality:\n";
  cntm.check_orthog(true);

  if (fname != "") {
    std::cout << "Writing orbitals (f) to " << fname << "\n";
    std::ofstream of(fname);
    const auto &gr = *wf.rgrid;
    of << "r ";
    for (const auto &Fc : cntm.orbitals) {
      of << "\"" << Fc.symbol(true) << "\" ";
    }
    of << "\n";
    for (std::size_t i = 0; i < gr.num_points(); i++) {
      of << gr.r(i) << " ";
      for (const auto &Fc : cntm.orbitals) {
        of << Fc.f(i) << " ";
      }
      of << "\n";
    }
  }
}

} // namespace Module
