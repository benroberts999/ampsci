#pragma once
#include "Coulomb/QkTable.hpp"
#include "Coulomb/meTable.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <string>
#include <vector>

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Just a test: example for playing with VQE
void VQE(const IO::InputBlock &input, const Wavefunction &wf);

//! Takes a subset of input basis according to subset_string.
//! Only includes states not in frozen_core_string
std::vector<DiracSpinor> basis_subset(const std::vector<DiracSpinor> &basis,
                                      const std::string &subset_string,
                                      const std::string &frozen_core_string);

double run_CI(const std::string &atom_name,
              const std::vector<DiracSpinor> &ci_sp_basis, int twoJ, int parity,
              int num_solutions, const Coulomb::meTable<double> &h1,
              const Coulomb::QkTable &qk, double e0, bool write_integrals,
              bool include_Sigma2,
              const std::vector<DiracSpinor> &mbpt_basis = {},
              double E_Fermi = 0.0, int min_n = 1);

class CSF2;

//! Forms list of 2-particle CSFs with given symmetry
std::vector<CSF2> form_CSFs(int twoJ, int parity,
                            const std::vector<DiracSpinor> &ci_sp_basis);

void write_CoulombIntegrals(const std::vector<DiracSpinor> &ci_sp_basis,
                            const Coulomb::QkTable &qk);

void write_CSFs(const std::vector<CSF2> &CSFs, int twoJ,
                const std::string &csf_fname);

//! Calculates the anti-symmetrised Coulomb integral for 2-particle states:
//! C1*C2*(g_abcd-g_abdc), where Cs are C.G. coefs
double CSF2_Coulomb(const Coulomb::QkTable &qk, const DiracSpinor &a,
                    const DiracSpinor &b, const DiracSpinor &c,
                    const DiracSpinor &d, int twoJ);

//! Determines CI Hamiltonian matrix element for two 2-particle CSFs, a and b
double Hab(const CSF2 &A, const CSF2 &B, int twoJ,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk);

//! CI Hamiltonian matrix element for two 2-particle CSFs: diagonal case
double Hab_0(const CSF2 &A, const CSF2 &B, int twoJ,
             const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk);

//! CI Hamiltonian matrix element for two 2-particle CSFs: differ by 1 case
double Hab_1(const CSF2 &A, const CSF2 &B, int twoJ,
             const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk);

//! CI Hamiltonian matrix element for two 2-particle CSFs: differ by 2 case
double Hab_2(const CSF2 &A, const CSF2 &B, int twoJ,
             const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk);

} // namespace Module
