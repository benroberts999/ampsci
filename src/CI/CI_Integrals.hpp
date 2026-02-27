#pragma once
#include "CSF.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/meTable.hpp"
#include "LinAlg/Matrix.hpp"
#include <string>
#include <vector>
class DiracDiracSpinor;
namespace MBPT {
class CorrelationPotential;
}
namespace HF {
class Breit;
}

//! Functions and classes for Configuration Interaction calculations
namespace CI {

//! Calculates the anti-symmetrised Coulomb integral for 2-particle states:
//! C1*C2*(g_abcd-g_abdc), where Cs are C.G. coefficients
double CSF2_Coulomb(const Coulomb::QkTable &qk, DiracSpinor::Index v,
                    DiracSpinor::Index w, DiracSpinor::Index x,
                    DiracSpinor::Index y, int twoJ);

//! Calculates the correlation [Sigma(2)] correction to CSF2_Coulomb().
//! Sk is the table of two-body single-particle matrix elements of Sigma_2.
double CSF2_Sigma2(const Coulomb::LkTable &Sk, DiracSpinor::Index v,
                   DiracSpinor::Index w, DiracSpinor::Index x,
                   DiracSpinor::Index y, int twoJ);

//! Calculates the anti-symmetrised Breit integral for 2-particle states:
//! C1*C2*(b_abcd-b_abdc), where Cs are C.G. coefficients
double CSF2_Breit(const Coulomb::WkTable &Bk, DiracSpinor::Index v,
                  DiracSpinor::Index w, DiracSpinor::Index x,
                  DiracSpinor::Index y, int twoJ);

//! Determines CI Hamiltonian matrix element for two 2-particle CSFs, a and b.
//! Does NOT include Sigma(2) matrix, but does include Sigma1 (if it's in h1)
double Hab(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk);

//! Calculates the Sigma(2) correction to Hab()
double Sigma2_AB(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
                 const Coulomb::LkTable &Sk);

//! Calculates the Breit correction to Hab()
double Breit_AB(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
                const Coulomb::WkTable &Bk);

//! Calculates table of single-particle matrix elements of one-body Hamiltonian.
//! Note: assumes basis are Hartree-Fock eigenstates!
[[nodiscard]] Coulomb::meTable<double>
calculate_h1_table(const std::vector<DiracSpinor> &ci_basis,
                   const std::vector<DiracSpinor> &s1_basis_core,
                   const std::vector<DiracSpinor> &s1_basis_excited,
                   const Coulomb::QkTable &qk, bool include_Sigma1);

[[nodiscard]] Coulomb::meTable<double>
calculate_h1_table(const std::vector<DiracSpinor> &ci_basis,
                   const MBPT::CorrelationPotential &Sigma,
                   bool include_Sigma1);

//! Calculates table of single-particle matrix elements of two-body Sigma_2 operator.
[[nodiscard]] Coulomb::LkTable calculate_Sk(
  const std::string &filename, const std::vector<DiracSpinor> &cis2_basis,
  const std::vector<DiracSpinor> &s2_basis_core,
  const std::vector<DiracSpinor> &s2_basis_excited, const Coulomb::QkTable &qk,
  int max_k, bool exclude_wrong_parity_box, bool no_new_integrals = false);

//! Calculates table of single-particle matrix elements of two-body Breit operator.
[[nodiscard]] Coulomb::WkTable
calculate_Bk(const std::string &bk_filename, const HF::Breit *const pBr,
             const std::vector<DiracSpinor> &ci_basis, int max_k,
             bool no_new_integralsQ = false);

//! Takes a subset of input basis according to subset_string.
//! Only states *not* included in frozen_core_string are included.
[[nodiscard]] std::vector<DiracSpinor>
basis_subset(const std::vector<DiracSpinor> &basis,
             const std::string &subset_string,
             const std::string &frozen_core_string = "");

//! Calculate reduced matrix elements between two CI states.
//! @details
//! cA is CI expansion coefficients (row if CI eigenvector matrix).
//! h is lookup table of single-particle matrix elements of operator.
double ReducedME(const LinAlg::View<const double> &cA,
                 const std::vector<CI::CSF2> &CSFAs, int twoJA,
                 const LinAlg::View<const double> &cB,
                 const std::vector<CI::CSF2> &CSFBs, int twoJB,
                 const Coulomb::meTable<double> &h, int K_rank, int Parity);

//! Calculate reduced matrix elements between two CI states.
//! @details
//! As is a PsiJPi (set of CI solutions for given J and pi),
//! iA is which solution to calculate for.
//! h is lookup table of single-particle matrix elements of operator.
inline double ReducedME(const PsiJPi &As, std::size_t iA, const PsiJPi &Bs,
                        std::size_t iB, const Coulomb::meTable<double> &h,
                        int K_rank, int Parity) {
  return ReducedME(As.coefs(iA), As.CSFs(), As.twoJ(), Bs.coefs(iB), Bs.CSFs(),
                   Bs.twoJ(), h, K_rank, Parity);
}

//! Calculate reduce ME between two 2-particle CSFs - XXX not quite right??
double RME_CSF2(const CI::CSF2 &X, int twoJX, const CI::CSF2 &V, int twoJV,
                const Coulomb::meTable<double> &h, int K_rank);

//! Determines best-fit for S and L for two-electron state by matching g-factor
std::pair<int, int> Term_S_L(int l1, int l2, int twoJ, double gJ_target);

//! Returns Term_Symbol as string
std::string Term_Symbol(int two_J, int L, int two_S, int parity);
//! Returns Term_Symbol as string, without J part
std::string Term_Symbol(int L, int two_S, int parity);

//! Constructs the CI matrix, optionally including Sigma corrections.
//! @details
//! h1 is table of one-body <a|h1|b> matrix elements, and may include Sigma_1.
//! qk contains two-body Coulomb integrals <vw|q^k|xy>.
//! Sk is optional; if given, constains two-body <vw|Sigma_2|xy> corrections.
LinAlg::Matrix<double> construct_Hci(const PsiJPi &psi,
                                     const Coulomb::meTable<double> &h1,
                                     const Coulomb::QkTable &qk,
                                     const Coulomb::WkTable *Bk = nullptr,
                                     const Coulomb::LkTable *Sk = nullptr);

} // namespace CI