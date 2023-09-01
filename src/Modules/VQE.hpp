#pragma once
#include "CI/CSF.hpp"
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

struct CIlevel {
  std::string config;
  int twoJ;
  int parity;
  int index;
  double e;
  double gJ;
  double L;
  double tSp1;
};

//==============================================================================
std::vector<CIlevel>
run_CI(const std::string &atom_name,
       const std::vector<DiracSpinor> &ci_sp_basis,
       const std::map<nkm, int> &orbital_map, int twoJ, int parity,
       int num_solutions, const Coulomb::meTable<double> &h1,
       const Coulomb::QkTable &qk, const Coulomb::LkTable &Sk,
       bool write_integrals, bool include_Sigma2,
       const std::vector<DiracSpinor> &mbpt_basis = {}, double E_Fermi = 0.0,
       int min_n = 1, const std::string &ci_input = "");

void write_CoulombIntegrals(const std::vector<DiracSpinor> &ci_sp_basis,
                            const std::map<nkm, int> &orbital_map,
                            const Coulomb::QkTable &qk);

void write_CSFs(const std::vector<CI::CSF2> &CSFs, int twoJ,
                const std::map<nkm, int> &orbital_map,
                const std::string &csf_fname);

} // namespace Module
