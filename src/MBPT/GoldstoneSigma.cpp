#include "MBPT/GoldstoneSigma.hpp"
#include "Angular/CkTable.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/SafeProfiler.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "MBPT/GreenMatrix.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <algorithm>
#include <numeric>

namespace MBPT {

//******************************************************************************
GoldstoneSigma::GoldstoneSigma(const HF::HartreeFock *const in_hf,
                               const std::vector<DiracSpinor> &basis,
                               const Sigma_params &sigp,
                               const rgrid_params &subgridp,
                               // const std::vector<DiracSpinor> &valence,
                               const std::string &fname)
    : CorrelationPotential(in_hf, basis, sigp, subgridp) {

  std::cout << "\nCorrelation potential (Sigma^2): Goldstone\n";

  const bool read_ok = read_write(fname, IO::FRW::read);

  print_subGrid();

  if (!read_ok) {
    std::cout << "Form correlation potential: Goldstone method\n";

    if (m_include_G)
      std::cout << "(Including FG/GF and GG)\n";

    if (!m_fk.empty()) {
      std::cout << "Effective screening factors: fk=";
      std::for_each(cbegin(m_fk), cend(m_fk),
                    [](auto x) { std::cout << x << ", "; });
      std::cout << "\n";
    }

    std::cout << "Basis: " << DiracSpinor::state_config(m_holes) << "/"
              << DiracSpinor::state_config(m_excited) << "\n";
  }
} // namespace MBPT

//******************************************************************************
void GoldstoneSigma::formSigma(int kappa, double en, int n) {
  // Calc dir + exchange
  // Print D, X, (D+X) energy shift
  // Add (D+X) to m_Sigma, and (n,k) to lookup_list
  // most of this is the same between each...?

  // already exists:
  const auto index = getSigmaIndex(n, kappa);
  if (index < m_Sigma_kappa.size()) {
    const auto [n2, k2, en2] = m_nk[index];
    // already have this (exact) potential?
    if (n == n2)
      return;
  }

  m_nk.emplace_back(n, kappa, en);
  auto &Sigma = m_Sigma_kappa.emplace_back(m_subgrid_points, m_include_G);

  // if v.kappa > basis, then Ck angular factor won't exist!
  if (Angular::twoj_k(kappa) > m_yeh.Ck().max_tj()) {
    std::cout << "Warning: angular not good\n";
    return;
  }

  printf("k=%2i at en=%8.5f.. ", kappa, en);
  std::cout << std::flush;
  // this->fill_Sigma_k(&Sigma, kappa, en);

  // auto Gmat_D = Sigma; // copy
  auto Gmat_X = Sigma; // copy
  Sigma2(&Sigma, &Gmat_X, kappa, en);

  // find lowest excited state, output <v|S|v> energy shift:
  const auto find_kappa = [kappa, n](const auto &a) {
    return a.k == kappa && (a.n == n || n == 0);
  };
  const auto vk = std::find_if(cbegin(m_excited), cend(m_excited), find_kappa);
  if (vk != cend(m_excited)) {
    auto deD = *vk * act_G_Fv(Sigma, *vk);
    auto deX = *vk * act_G_Fv(Gmat_X, *vk);
    // nb: just approximate (uses splines)
    printf("de= %7.1f + %5.1f = ", deD * PhysConst::Hartree_invcm,
           deX * PhysConst::Hartree_invcm);
    printf("%7.1f", (deD + deX) * PhysConst::Hartree_invcm);
  }

  Sigma += Gmat_X;

  std::cout << "\n";
}

//******************************************************************************
void GoldstoneSigma::Sigma2(GMatrix *Gmat_D, GMatrix *Gmat_X, int kappa,
                            double en) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  // Four second-order diagrams:
  // Diagram (a):
  // |Q^k_amn><Q^k_amn| / de_amn / [k][j]
  // Diagram (b) (exchange):
  // |Q^k_amn><P^k_amn| / de_amn / [k][j]
  // Diagram (c):
  // |Q^k_nba><Q^k_nba| / de_nba / [k][j]
  // Diagram (d) (exchange):
  // |Q^k_nba><P^k_nba| / de_nba / [k][j]
  // where:
  // All indeces are summed over,
  // a & b are core states, n & m are virtual excited states,
  // k is multipolarity [Coloulmb expansion]
  // de_xyz = e_v + e_x - e_y - e_z

  // can do only D/X part if either is null?
  // In which case, this should be in base class?

  const auto &Ck = m_yeh.Ck();

  if (m_holes.empty())
    return;

  std::vector<GMatrix> Gds(m_holes.size(), {m_subgrid_points, m_include_G});
  std::vector<GMatrix> Gxs(m_holes.size(), {m_subgrid_points, m_include_G});
#pragma omp parallel for
  for (auto ia = 0ul; ia < m_holes.size(); ia++) {
    const auto &a = m_holes[ia];
    auto &Ga_d = Gds[ia];
    auto &Ga_x = Gxs[ia];
    auto Qkv = DiracSpinor(0, kappa, p_gr); // re-use to reduce alloc'ns
    auto Pkv = DiracSpinor(0, kappa, p_gr); // re-use to reduce alloc'ns
    for (const auto &n : m_excited) {
      const auto [kmin_nb, kmax_nb] = Coulomb::k_minmax(n, a);
      for (int k = kmin_nb; k <= kmax_nb; ++k) {
        if (Ck(k, a.k, n.k) == 0)
          continue;
        const auto f_kkjj = (2 * k + 1) * (Angular::twoj_k(kappa) + 1);

        // Effective screening parameter:
        const auto fk = get_fk(k);
        if (fk == 0.0)
          continue;

        // Diagrams (a) [direct] and (b) [exchange]
        for (const auto &m : m_excited) {
          if (Ck(k, kappa, m.k) == 0)
            continue;
          Qkv = m_yeh.Qkv_bcd(Qkv.k, a, m, n, k);
          Pkv = m_yeh.Pkv_bcd(Pkv.k, a, m, n, k, m_fk);
          const auto dele = en + a.en() - m.en() - n.en();
          const auto factor = fk / (f_kkjj * dele);
          addto_G(&Ga_d, Qkv, Qkv, factor);
          addto_G(&Ga_x, Qkv, Pkv, factor);
        } // m

        // Diagrams (c) [direct] and (d) [exchange]
        for (const auto &b : m_holes) {
          if (Ck(k, kappa, b.k) == 0)
            continue;
          Qkv = m_yeh.Qkv_bcd(Qkv.k, n, b, a, k);
          Pkv = m_yeh.Pkv_bcd(Pkv.k, n, b, a, k, m_fk);
          const auto dele = en + n.en() - b.en() - a.en();
          const auto factor = fk / (f_kkjj * dele); // XXX
          addto_G(&Ga_d, Qkv, Qkv, factor);
          addto_G(&Ga_x, Qkv, Pkv, factor);
        } // b

      } // k
    }   // n
  }     // a

  // // note: no benefit to sending Gmat in! Just let it be a return value!
  for (const auto &Gd : Gds)
    *Gmat_D += Gd;
  for (const auto &Gx : Gxs)
    *Gmat_X += Gx;
}

} // namespace MBPT
