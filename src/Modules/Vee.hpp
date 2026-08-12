#pragma once
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace Module {
void Vee(const IO::InputBlock &input, const Wavefunction &wf);
void test_CI(const IO::InputBlock &input, const Wavefunction &wf);
double CI_edm(const Wavefunction &wf, const std::string basis_string,
              const double mu, const bool contact, const int J_V,
              const int pi_V, const int i_V, const bool verbose = true);
double calc_ME(const std::string op, const std::string type, const bool contact,
               const double mu, const bool tdhf, const double omega,
               const Wavefunction &wf, const DiracSpinor &Fw,
               const DiracSpinor &Fv);

} // namespace Module
