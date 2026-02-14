#pragma once
#include "CI/CSF.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/meTable.hpp"
#include "DiracOperator/include.hpp"
#include "LinAlg/Matrix.hpp"
#include "MBPT/Sigma2.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/ostream.hpp"
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

//------------------------------------------------------------------------------
void write_CSFs(const std::vector<CI::CSF2> &CSFs, int twoJ,
                const std::map<nkm, int> &orbital_map,
                const std::string &csf_fname);

void write_H(const LinAlg::Matrix<double> &Hci, const std::string &csf_fname);

//------------------------------------------------------------------------------
template <class Integrals>
void write_CoulombIntegrals(const std::string &fname,
                            const std::vector<DiracSpinor> &ci_sp_basis,
                            const std::map<nkm, int> &orbital_map,
                            const Integrals &qk) {

  // Modify this to allow for Sigma_2!

  static_assert(std::is_same_v<Integrals, Coulomb::QkTable> ||
                std::is_same_v<Integrals, Coulomb::LkTable>);

  // print all two-particle integrals Rk to file:
  std::ofstream g_file(fname);
  if (std::is_same_v<Integrals, Coulomb::QkTable>) {
    g_file << "# a  b  c  d  g_abcd    ## (nb: abcd = badc = cdab = dcba; only "
              "'smallest' is written)\n";
  } else {
    g_file << "# a  b  c  d  s_abcd    ## (nb: abcd = badc; only "
              "'smallest' is written)\n";
  }

  for (const auto &a : ci_sp_basis) {
    for (const auto &b : ci_sp_basis) {
      for (const auto &c : ci_sp_basis) {
        for (const auto &d : ci_sp_basis) {

          // Selection rules, for g/s
          const auto [k0, k1] = std::is_same_v<Integrals, Coulomb::QkTable> ?
                                  Coulomb::k_minmax_Q(a, b, c, d) :
                                  MBPT::k_minmax_S(a, b, c, d);
          if (k1 < k0)
            continue;

          for (int tma = -a.twoj(); tma <= a.twoj(); tma += 2) {
            for (int tmb = -b.twoj(); tmb <= b.twoj(); tmb += 2) {
              for (int tmc = -c.twoj(); tmc <= c.twoj(); tmc += 2) {
                for (int tmd = -d.twoj(); tmd <= d.twoj(); tmd += 2) {

                  // m = j_z selection rules:
                  if (tmc - tma != tmb - tmd)
                    continue;

                  const auto ia =
                    (uint16_t)orbital_map.at(nkm{a.n(), a.kappa(), tma});
                  const auto ib =
                    (uint16_t)orbital_map.at(nkm{b.n(), b.kappa(), tmb});
                  const auto ic =
                    (uint16_t)orbital_map.at(nkm{c.n(), c.kappa(), tmc});
                  const auto id =
                    (uint16_t)orbital_map.at(nkm{d.n(), d.kappa(), tmd});

                  // Equivilant integrals:
                  // abcd = badc = cdab = dcba for g
                  // abcd = badc               for s
                  // nb: this only works if largest of (ia,ib,ic,id)
                  // is smaller than 2^16, which is always true
                  const auto indexify = [](uint16_t w, uint16_t x, uint16_t y,
                                           uint16_t z) {
                    return ((uint64_t)w << 48) + ((uint64_t)x << 32) +
                           ((uint64_t)y << 16) + (uint64_t)z;
                  };
                  uint64_t i1 = indexify(ia, ib, ic, id);
                  uint64_t i2 = indexify(ib, ia, id, ic);
                  uint64_t i3 = indexify(ic, id, ia, ib);
                  uint64_t i4 = indexify(id, ic, ib, ia);

                  // Only include the unique ones:
                  if (std::is_same_v<Integrals, Coulomb::QkTable>) {
                    // g symmetry
                    if (i1 != std::min({i1, i2, i3, i4}))
                      continue;
                  } else {
                    // s symmetry
                    if (i1 != std::min({i1, i2}))
                      continue;
                  }

                  const auto g = qk.g(a, b, c, d, tma, tmb, tmc, tmd);
                  // Note sure why zero values slip through?
                  // Missing SR? or mistake?
                  if (g == 0.0)
                    continue;

                  fmt::print(g_file, "{} {} {} {} {:.8e}\n", ia, ib, ic, id, g);
                }
              }
            }
          }
        }
      }
    }
  }
}

} // namespace Module
