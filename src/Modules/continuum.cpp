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
               {"filename", "filename for output"}});

  auto en_c = input.get("ec", 0.5);
  auto lmax = input.get("max_l", 0);
  auto fname = input.get("filename", std::string{""});

  std::cout << "\nContinuum Orbitals:\n";
  std::cout << "energy: " << en_c << "\n";

  auto cntm = ContinuumOrbitals(wf);
  cntm.solveContinuumHF(en_c, lmax);

  std::cout << "Check orthogonanilty:\n";
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
