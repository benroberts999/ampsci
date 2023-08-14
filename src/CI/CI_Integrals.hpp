#pragma once
#include "CSF.hpp"
#include "Coulomb/Coulomb.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>

namespace CI {

//! Calculates the anti-symmetrised Coulomb integral for 2-particle states:
//! C1*C2*(g_abcd-g_abdc), where Cs are C.G. coefficients
double CSF2_Coulomb(const Coulomb::QkTable &qk, const DiracSpinor &a,
                    const DiracSpinor &b, const DiracSpinor &c,
                    const DiracSpinor &d, int twoJ);

double CSF2_Sigma2(const Coulomb::QkTable &qk, const DiracSpinor &a,
                   const DiracSpinor &b, const DiracSpinor &c,
                   const DiracSpinor &d, int twoJ,
                   const std::vector<DiracSpinor> &core,
                   const std::vector<DiracSpinor> &excited,
                   const Angular::SixJTable &SixJ);

//! Determines CI Hamiltonian matrix element for two 2-particle CSFs, a and b.
//! Does NOT include Sigma(2) matrix, but does include Sigma1 (if it's in h1)
double Hab(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk);

double Sigma2_AB(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
                 const Coulomb::QkTable &qk,
                 const std::vector<DiracSpinor> &core,
                 const std::vector<DiracSpinor> &excited,
                 const Angular::SixJTable &SixJ);

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

} // namespace CI