#include "Modules/VarAlph.hpp"
#include "DiracOperator/DiracOperator.hpp" //For E1 operator
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

void VarAlph(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check({});

  // Run CI twice, once for regular alpha, one for delta

  // Resturns a vector of CI solutions
  const auto JPsis = wf.CIwfs();

  std::cout << "\nHave CI solutions for:\n";
  for (const auto &JPsi : JPsis) {
    std::cout << JPsi.twoJ() / 2 << ", pi=" << JPsi.parity() << "\n";
  }

  // just first J/Pi, as example
  const auto &JPsi = JPsis.at(0);

  // CI coeficient for the ith CI solution, corresponding to the jth CSF:
  // JPsi.coef(i, j);
  // JPsi.coefs()

  const auto p_Hci = JPsi.Hci();
  if (p_Hci == nullptr) {
    std::cout << "We need to store the matrix!\n";
    return;
  }
  const auto &Hci = *p_Hci;

  for (std::size_t sol = 0; sol < JPsi.num_solutions(); ++sol) {

    double E_mat = 0.0;
    double E_diag = 0.0;
    // square matrix, so rows = cols anyway
    for (std::size_t i = 0; i < Hci.rows(); ++i) {
      for (std::size_t j = 0; j < Hci.cols(); ++j) {

        const auto c_i = JPsi.coef(sol, i);
        const auto c_j = JPsi.coef(sol, j);
        const auto H_ij = Hci(i, j);

        E_mat += c_i * c_j * H_ij;

        if (i == j) {
          E_diag += c_i * c_j * H_ij;
        }
      }
    }

    std::cout << sol << " " << JPsi.energy(sol) << " " << E_mat << " " << E_diag
              << "\n";
  }
}

} // namespace Module
