#pragma once
#include "Coulomb/QkTable.hpp"
#include "Coulomb/meTable.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "LinAlg/Matrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <string>
#include <vector>

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

// Stuct, hold n, kappa, and 2*m (for the orbital to index map)
struct nkm {
  int n;
  int kappa;
  int twom;

  //required for std::map (needs unique ordering; order doesn't matter)
  friend bool operator<(const nkm &lhs, const nkm &rhs) {
    if (lhs.n == rhs.n) {
      if (lhs.kappa == rhs.kappa) {
        return lhs.twom < rhs.twom;
      }
      return lhs.kappa < rhs.kappa;
    }
    return lhs.n < rhs.n;
  }
};

//! Just a test: example for playing with VQE
void VQE(const IO::InputBlock &input, const Wavefunction &wf);

//! Takes a subset of input basis according to subset_string.
//! Only includes states not in frozen_core_string
std::vector<DiracSpinor> basis_subset(const std::vector<DiracSpinor> &basis,
                                      const std::string &subset_string,
                                      const std::string &frozen_core_string);

double run_CI(const std::string &atom_name,
              const std::vector<DiracSpinor> &ci_sp_basis,
              const std::map<nkm, int> &orbital_map, int twoJ, int parity,
              int num_solutions, const Coulomb::meTable<double> &h1,
              const Coulomb::QkTable &qk, double e0, bool write_integrals,
              bool include_Sigma2,
              const std::vector<DiracSpinor> &mbpt_basis = {},
              double E_Fermi = 0.0, int min_n = 1,
              const std::string &ci_input = "");

class CSF2;

//! Forms list of 2-particle CSFs with given symmetry
std::vector<CSF2> form_CSFs(int twoJ, int parity,
                            const std::vector<DiracSpinor> &ci_sp_basis);

void write_CoulombIntegrals(const std::vector<DiracSpinor> &ci_sp_basis,
                            const std::map<nkm, int> &orbital_map,
                            const Coulomb::QkTable &qk);

void write_CSFs(const std::vector<CSF2> &CSFs, int twoJ,
                const std::map<nkm, int> &orbital_map,
                const std::string &csf_fname);

//! Calculates the anti-symmetrised Coulomb integral for 2-particle states:
//! C1*C2*(g_abcd-g_abdc), where Cs are C.G. coefficients
double CSF2_Coulomb(const Coulomb::QkTable &qk, const DiracSpinor &a,
                    const DiracSpinor &b, const DiracSpinor &c,
                    const DiracSpinor &d, int twoJ);

//! Determines CI Hamiltonian matrix element for two 2-particle CSFs, a and b
double Hab(const CSF2 &A, const CSF2 &B, int twoJ,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk);

//! Calculate reduced matrix elements between two CI states. cA is CI expansion coefficients (row if CI eigenvector matrix)
double CI_RME(const double *cA, const std::vector<CSF2> &CSFAs, int twoJA,
              const double *cB, const std::vector<CSF2> &CSFBs, int twoJB,
              const DiracOperator::TensorOperator *h);

//! Overload for diagonal matrix elements
double CI_RME(const double *cA, const std::vector<CSF2> &CSFs, int twoJ,
              const DiracOperator::TensorOperator *h);

//! Calculate reduce ME between two 2-particle CSFs - XXX not quite right??
double RME_CSF2(const CSF2 &V, int twoJV, const CSF2 &X, int twoJX,
                const DiracOperator::TensorOperator *h);

} // namespace Module
