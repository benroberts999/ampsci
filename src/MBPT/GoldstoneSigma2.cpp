#include "MBPT/GoldstoneSigma2.hpp"
#include "Angular/Angular_tables.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/SafeProfiler.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "MBPT/GreenMatrix.hpp"
#include "Maths/Grid.hpp"
#include "Maths/LinAlg_MatrixVector.hpp"
#include <algorithm>
#include <numeric>

namespace MBPT {

//******************************************************************************
GoldstoneSigma2::GoldstoneSigma2(const HF::HartreeFock *const in_hf,
                                 const std::vector<DiracSpinor> &basis,
                                 const Sigma_params &sigp,
                                 const rgrid_params &subgridp,
                                 const std::vector<double> &en_list,
                                 const std::string &atom)
    : CorrelationPotential(in_hf, basis, sigp, subgridp) {

  std::cout << "\nCorrelation potential (Sigma^2): Goldstone\n";

  const auto fname =
      atom == "" ? ""
                 : atom + "_" + std::to_string(p_gr->num_points) + ".SigmaG";
  const bool read_ok = read_write(fname, IO::FRW::read);

  if (en_list.empty())
    return; //?

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

    form_Sigma(en_list, fname);
  }
} // namespace MBPT

//******************************************************************************
void GoldstoneSigma2::fill_Sigma_k(GMatrix *Gmat, const int kappa,
                                   const double en) {
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

  const auto &Ck = m_yeh.Ck();

  // Just for safety, should already be zero (unless re-calcing G)
  Gmat->zero();
  // This is only so I can print each energy shift:
  auto Sdir = *Gmat; // blank copies!
  auto Sexch = *Gmat;

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
      const auto [kmin_nb, kmax_nb] = m_yeh.k_minmax(n, a);
      for (int k = kmin_nb; k <= kmax_nb; ++k) {
        if (Ck(k, a.k, n.k) == 0)
          continue;
        const auto f_kkjj = (2 * k + 1) * (Angular::twoj_k(kappa) + 1);
        const auto &yknb = m_yeh(k, n, a);

        // Effective screening parameter:
        const auto fk = get_fk(k);

        // Diagrams (a) [direct] and (b) [exchange]
        for (const auto &m : m_excited) {
          if (Ck(k, kappa, m.k) == 0)
            continue;
          Coulomb::Qkv_bcd(&Qkv, a, m, n, k, yknb, Ck);
          // Coulomb::Pkv_bcd(&Pkv, a, m, n, k, m_yeh(m, a), Ck, m_6j);
          // Pkv_bcd_2 allows different screening factor for each 'k2' in exch.
          Coulomb::Pkv_bcd_2(&Pkv, a, m, n, k, m_yeh(m, a), Ck, m_6j, m_fk);
          const auto dele = en + a.en - m.en - n.en;
          const auto factor = fk / (f_kkjj * dele);
          addto_G(&Ga_d, Qkv, Qkv, factor);
          addto_G(&Ga_x, Qkv, Pkv, factor);
        } // m

        // Diagrams (c) [direct] and (d) [exchange]
        for (const auto &b : m_holes) {
          if (Ck(k, kappa, b.k) == 0)
            continue;
          Coulomb::Qkv_bcd(&Qkv, n, b, a, k, yknb, Ck);
          // Coulomb::Pkv_bcd(&Pkv, n, b, a, k, m_yeh(n, b), Ck, m_6j);
          Coulomb::Pkv_bcd_2(&Pkv, n, b, a, k, m_yeh(n, b), Ck, m_6j, m_fk);
          const auto dele = en + n.en - b.en - a.en;
          const auto factor = fk / (f_kkjj * dele); // XXX
          addto_G(&Ga_d, Qkv, Qkv, factor);
          addto_G(&Ga_x, Qkv, Pkv, factor);
        } // b

      } // k
    }   // n
  }     // a

  // // note: no benefit to sending Gmat in! Just let it be a return value!
  // This is only so I can print each energy shift!
  for (const auto &Gd : Gds)
    Sdir += Gd;
  for (const auto &Gx : Gxs)
    Sexch += Gx;

  *Gmat = Sdir + Sexch;

  // XXX This is only so I can print each energy shift!
  auto find_kappa = [=](const auto &f) { return f.k == kappa; };
  const auto Fk = std::find_if(cbegin(m_excited), cend(m_excited), find_kappa);
  if (Fk != cend(m_excited)) {
    auto deD = *Fk * act_G_Fv(Sdir, *Fk);
    auto deX = *Fk * act_G_Fv(Sexch, *Fk);
    printf("de= %.4f + %.5f = ", deD, deX);
  }
}

} // namespace MBPT
