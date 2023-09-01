#pragma once
#include "CSF.hpp"
#include "Coulomb/Coulomb.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <string>
#include <vector>

namespace CI {

//! Calculates the anti-symmetrised Coulomb integral for 2-particle states:
//! C1*C2*(g_abcd-g_abdc), where Cs are C.G. coefficients
double CSF2_Coulomb(const Coulomb::QkTable &qk, const DiracSpinor &v,
                    const DiracSpinor &w, const DiracSpinor &x,
                    const DiracSpinor &y, int twoJ);

double CSF2_Sigma2(const Coulomb::LkTable &Sk, const DiracSpinor &v,
                   const DiracSpinor &w, const DiracSpinor &x,
                   const DiracSpinor &y, int twoJ);

//! Determines CI Hamiltonian matrix element for two 2-particle CSFs, a and b.
//! Does NOT include Sigma(2) matrix, but does include Sigma1 (if it's in h1)
double Hab(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk);

double Sigma2_AB(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
                 const Coulomb::LkTable &Sk);

Coulomb::meTable<double>
calculate_h1_table(const std::vector<DiracSpinor> &ci_basis,
                   const std::vector<DiracSpinor> &s1_basis_core,
                   const std::vector<DiracSpinor> &s1_basis_excited,
                   const Coulomb::QkTable &qk, bool include_Sigma1);

Coulomb::LkTable calculate_Sk(const std::string &filename,
                              const std::vector<DiracSpinor> &cis2_basis,
                              const std::vector<DiracSpinor> &s2_basis_core,
                              const std::vector<DiracSpinor> &s2_basis_excited,
                              const Coulomb::QkTable &qk, int max_delta_n,
                              bool exclude_wrong_parity_box);

//! Calculate reduced matrix elements between two CI states. cA is CI expansion coefficients (row if CI eigenvector matrix)
double ReducedME(const double *cA, const std::vector<CI::CSF2> &CSFAs,
                 int twoJA, const double *cB,
                 const std::vector<CI::CSF2> &CSFBs, int twoJB,
                 const DiracOperator::TensorOperator *h);

//! Overload for diagonal matrix elements
double ReducedME(const double *cA, const std::vector<CI::CSF2> &CSFs, int twoJ,
                 const DiracOperator::TensorOperator *h);

//! Calculate reduce ME between two 2-particle CSFs - XXX not quite right??
double RME_CSF2(const CI::CSF2 &V, int twoJV, const CI::CSF2 &X, int twoJX,
                const DiracOperator::TensorOperator *h);

//! Determines best-fit for S and L for two-electron state by matching g-factor
std::pair<int, int> Term_S_L(int l1, int l2, int twoJ, double gJ_target);

} // namespace CI