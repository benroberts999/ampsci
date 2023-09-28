#pragma once
#include "CSF.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/meTable.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <vector>

namespace CI {

//! Runs configuration interaction with options specified by input.
//! Returns a list of PsiJPi (set of CI solutions for list of J/Pi)
std::vector<PsiJPi> configuration_interaction(const IO::InputBlock &input,
                                              const Wavefunction &wf);

//! Performs CI for specified J and Pi
PsiJPi run_CI(const std::vector<DiracSpinor> &ci_sp_basis, int twoJ, int parity,
              int num_solutions, const Coulomb::meTable<double> &h1,
              const Coulomb::QkTable &qk, const Coulomb::LkTable &Sk,
              bool include_Sigma2);

} // namespace CI