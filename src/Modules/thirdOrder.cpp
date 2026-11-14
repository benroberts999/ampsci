#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "IO/InputBlock.hpp"
#include "LinAlg/include.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "MBPT/Ladder.hpp"
#include "Modules/Modules.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace Module {

// Declare, register, then define below.
void thirdOrder(const IO::InputBlock &input, const Wavefunction &wf);
namespace {
const Register r_thirdOrder{
  "thirdOrder",
  "Third-order correlation corrections to valence energies, by class of "
  "diagram",
  &thirdOrder};
} // namespace

//==============================================================================
/*!
  @brief Third-order correlation correction to the energy of a valence state,
         by class of diagram.

  @details
  Lowest-order (not iterated) evaluation of every third-order diagram for the
  energy of a single valence state \f$v\f$, with Hartree-Fock orbitals (so
  there are no one-body terms); see Blundell, Guo, Johnson, Sapirstein,
  At. Data Nucl. Data Tables 37, 103 (1987). The diagrams are sorted into
  classes, to show the relative size of each.

  Notation: \f$a,b,c\f$ are core orbitals, \f$m,n,r\f$ excited orbitals,
  \f$[k] = 2k+1\f$, and \f$W^k = Q^k + P^k\f$ is the antisymmetrised Coulomb
  integral: \f$Q^k\f$ is the direct part and \f$P^k\f$ the exchange part.
  Functions that return {direct, exchange} do so in that order.
*/
namespace ThirdOrder {

//! Direct and exchange Coulomb integrals: {Q^k_abcd, P^k_abcd}
std::array<double, 2>
direct_exchange(int k, const DiracSpinor &a, const DiracSpinor &b,
                const DiracSpinor &c, const DiracSpinor &d,
                const Coulomb::QkTable &qk, const Angular::SixJTable &sjt);

/*!
  @brief Third-order terms in which the three interactions are linked through
         two particle-hole pairs, by class.
  @details See particle_hole_terms(). The two valence vertices involve the
           valence electron; the core vertex links the two pairs.
*/
struct ParticleHoleTerms {
  //! Screening: all vertices direct
  double S_d{0.0};
  //! Screening: one valence vertex exchange
  double S_x{0.0};
  //! Screening: both valence vertices exchange
  double S_xx{0.0};
  //! Hole-particle interaction (core vertex exchange): valence vertices direct
  double H_d{0.0};
  //! Hole-particle interaction: one valence vertex exchange
  double H_x{0.0};
  //! Crossed ladder: all vertices exchange
  double L_cross{0.0};
};

/*!
  @brief Third-order energy terms in which the three interactions are linked
         through two particle-hole pairs.

  @details
  For each multipolarity \f$k\f$, the particle-hole pairs \f$(a,n)\f$ and
  \f$(b,m)\f$ are linked by the core vertices
  \f[
    G_{an,bm} = W^k_{bnma},
  \f]
  \f[
    T_{an,bm} = \frac{W^k_{abnm}}{\en_a+\en_b-\en_n-\en_m},
  \f]
  and to the valence electron by the valence vertices. With the valence
  electron scattered to an excited orbital \f$x\f$:
  \f[
    t_{x,an} = \frac{W^k_{vaxn}}{\en_v+\en_a-\en_x-\en_n},
  \f]
  \f[
    u_{x,an} = W^k_{vnxa},
  \f]
  and to a core orbital \f$x\f$:
  \f[
    t_{x,an} = \frac{W^k_{vnxa}}{\en_x+\en_a-\en_v-\en_n},
  \f]
  \f[
    u_{x,an} = W^k_{vaxn}.
  \f]
  The energy is
  \f[
    \delta\en_v = \sum_k \frac{(-1)^k}{[k]^2[j_v]} \sum_x (\pm)
      \left( t_x^T\,G\,t_x + 2\,t_x^T\,T\,u_x \right),
  \f]
  with the plus sign for excited \f$x\f$ and the minus sign for core \f$x\f$.

  Each \f$W = Q + P\f$: the eight combinations of direct and exchange vertices
  are sorted into the classes of ParticleHoleTerms.

  @param v        Valence orbital
  @param core     Core orbitals
  @param excited  Excited orbitals
  @param qk       Coulomb integral table
  @param sjt      6j symbol table
*/
ParticleHoleTerms particle_hole_terms(const DiracSpinor &v,
                                      const std::vector<DiracSpinor> &core,
                                      const std::vector<DiracSpinor> &excited,
                                      const Coulomb::QkTable &qk,
                                      const Angular::SixJTable &sjt);

/*!
  @brief Lowest-order ladder integrals: the ladder diagrams with Coulomb
         integrals for both interactions.
  @details See MBPT::calculate_Lk().
  @param qk         Coulomb integral table
  @param excited    Excited orbitals
  @param core       Core orbitals
  @param integrals  The set of integrals L^k_mnij to calculate
  @param diagrams   Which of the four ladder diagrams to include
*/
Coulomb::LkTable lowest_order_ladder(const Coulomb::QkTable &qk,
                                     const std::vector<DiracSpinor> &excited,
                                     const std::vector<DiracSpinor> &core,
                                     const MBPT::LadderIntegralSet &integrals,
                                     const MBPT::LadderDiagrams &diagrams);

/*!
  @brief Energy from ladder integrals with the valence orbital in the third
         index, L^k_mnva. Returns {direct, exchange}.
  @details
  \f[
    \delta\en_v = \sum_{k,a,m,n}
      \frac{F^k_{vamn}\,L^k_{mnva}}{[k][j_v]\,(\en_v+\en_a-\en_m-\en_n)},
  \f]
  with \f$F = Q\f$ (direct) or \f$F = P\f$ (exchange).
*/
std::array<double, 2>
ladder_energy_mnva(const DiracSpinor &v, const std::vector<DiracSpinor> &core,
                   const std::vector<DiracSpinor> &excited,
                   const Coulomb::QkTable &qk, const Coulomb::LkTable &lk,
                   const Angular::SixJTable &sjt);

/*!
  @brief Energy from ladder integrals with the valence orbital in the first
         index, L^k_vnab. Returns {direct, exchange}.
  @details
  \f[
    \delta\en_v = \sum_{k,n,a,b}
      \frac{F^k_{vnab}\,L^k_{vnab}}{[k][j_v]\,(\en_v+\en_n-\en_a-\en_b)},
  \f]
  with \f$F = Q\f$ (direct) or \f$F = P\f$ (exchange).
*/
std::array<double, 2>
ladder_energy_vnab(const DiracSpinor &v, const std::vector<DiracSpinor> &core,
                   const std::vector<DiracSpinor> &excited,
                   const Coulomb::QkTable &qk, const Coulomb::LkTable &lk,
                   const Angular::SixJTable &sjt);

/*!
  @brief Third-order relaxation terms, by class.
  @details See relaxation_terms(). First label: the interaction of the valence
           electron with the core is direct (d) or exchange (x). Second label:
           the two interactions that change the core are both direct (d), or
           one of them is exchange (x).
*/
struct RelaxationTerms {
  //! Valence-core interaction direct; core interactions both direct
  double R_dd{0.0};
  //! Valence-core interaction direct; one core interaction exchange
  double R_dx{0.0};
  //! Valence-core interaction exchange; core interactions both direct
  double R_xd{0.0};
  //! Valence-core interaction exchange; one core interaction exchange
  double R_xx{0.0};
};

/*!
  @brief Third-order relaxation terms: the valence electron interacts with the
         second-order change of the core.

  @details
  \f[
    \delta\en_v = \sum_{pq} \frac{(-1)^{j_v+j_p+1}}{\sqrt{[j_v][j_p]}}\,
                  W^0_{vpvq}\,\rho_{qp},
  \f]
  where \f$p\f$ and \f$q\f$ have the same angular quantum numbers, and
  \f$\rho\f$ is the second-order change of the density (symmetric):
  \f[
    \rho_{ma} = \frac{1}{\en_a-\en_m}\sum_{k}\frac{1}{[k]}\left(
      \sum_{bnr}\frac{Q^k_{mbnr}W^k_{nrab}}{\en_a+\en_b-\en_n-\en_r}
      - \sum_{bcn}\frac{Q^k_{bcan}W^k_{mnbc}}{\en_b+\en_c-\en_m-\en_n}
      \right),
  \f]
  \f[
    \rho_{bc} = -\sum_{k,amn}\frac{Q^k_{abmn}W^k_{mnac}}
      {[k]\,(\en_a+\en_b-\en_m-\en_n)(\en_a+\en_c-\en_m-\en_n)},
  \f]
  \f[
    \rho_{mr} = \sum_{k,abn}\frac{Q^k_{abmn}W^k_{nrba}}
      {[k]\,(\en_a+\en_b-\en_m-\en_n)(\en_a+\en_b-\en_n-\en_r)}.
  \f]
  Each \f$W = Q + P\f$: sorted into the classes of RelaxationTerms.

  @param v        Valence orbital
  @param core     Core orbitals
  @param excited  Excited orbitals
  @param qk       Coulomb integral table
  @param sjt      6j symbol table
*/
RelaxationTerms relaxation_terms(const DiracSpinor &v,
                                 const std::vector<DiracSpinor> &core,
                                 const std::vector<DiracSpinor> &excited,
                                 const Coulomb::QkTable &qk,
                                 const Angular::SixJTable &sjt);

} // namespace ThirdOrder
using namespace ThirdOrder;

//==============================================================================

void thirdOrder(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
    {{"", "Third-order corrections to valence energies:\n"
          "each class of diagram separately (not iterated).\n"
          " - S_d, S_x, S_xx: screening\n"
          " - H_d, H_x: hole-particle interaction\n"
          " - L_cross: crossed ladder\n"
          " - L_par, L_parx: ladder\n"
          " - R_dd, R_dx, R_xd, R_xx: core relaxation\n"
          " - d: direct; x: exchange\n"
          "Also prints the energy of the Ladder{} diagrams:\n"
          "3rd order, and iterated (if Lk_file exists);\n"
          "and <v|Sigma|v> (if Correlations included)"},
     {"min_n_core", "Lowest core n to include [1]"},
     {"basis", "Basis string: which excited states to include, e.g., "
               "'30spdf' [entire basis]"},
     {"max_k", "Maximum k of the Coulomb and ladder integrals. -1 for all "
               "[8]"},
     {"Qk_file", "Filename for Qk Coulomb integrals: read if exists, "
                 "otherwise calculated and written. 'false' to not read or "
                 "write [<identity>.qk.abf]"},
     {"Lk_file", "Filename of Lk ladder integrals, from Ladder{} with "
                 "include_L4 and each_valence true "
                 "[<identity>_l4.lk.abf]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  using namespace std::string_literals;
  const auto min_n_core = input.get("min_n_core", 1);
  const auto basis_string = input.get("basis", ""s);
  const auto max_k = input.get("max_k", 8);
  const auto Qk_file = input.get("Qk_file", wf.identity() + ".qk.abf"s);
  const auto Lk_file = input.get("Lk_file", wf.identity() + "_l4.lk.abf"s);

  // Core and excited orbitals: from the basis
  const auto core_excited =
    DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel(), min_n_core);
  const auto &core = core_excited.first;
  const auto excited =
    basis_string.empty() ?
      core_excited.second :
      DiracSpinor::subset(core_excited.second, basis_string, false);

  // Valence states: the basis versions
  std::vector<DiracSpinor> valence;
  for (const auto &Fv : wf.valence()) {
    const auto basis_v = std::find(excited.cbegin(), excited.cend(), Fv);
    if (basis_v == excited.cend()) {
      std::cout << "Warning: Basis missing valence state: " << Fv << "\n";
    } else {
      valence.push_back(*basis_v);
    }
  }
  if (valence.empty() || core.empty())
    return;

  std::cout << "\nThird-order correlation corrections to valence energies\n";
  std::cout << "Core   : " << DiracSpinor::state_config(core) << "\n";
  std::cout << "Excited: " << DiracSpinor::state_config(excited) << "\n";

  // Coulomb integrals
  const auto core_and_excited = qip::merge(core, excited);
  const Coulomb::YkTable yk(core_and_excited);
  const auto &sjt = yk.SixJ();
  Coulomb::QkTable qk;
  std::cout << "\nFill (or read) Qk table:\n";
  qk.read(Qk_file);
  const auto num_read = qk.count();
  qk.fill(core_and_excited, yk, max_k);
  if (qk.count() != num_read) {
    qk.write(Qk_file);
  }

  // Iterated ladder integrals, if they have been calculated
  Coulomb::LkTable lk;
  const auto have_ladder = lk.read(Lk_file);

  // Lowest-order ladder integrals with the valence orbital in the third
  // index, L_mnva (for all valence states): the particle-particle (L1),
  // particle-hole (L2 + L3), and hole-hole (L4) diagrams
  const MBPT::LadderDiagrams diagram_L1{true, false, false, false};
  const MBPT::LadderDiagrams diagrams_L2_L3{false, true, true, false};
  const MBPT::LadderDiagrams diagram_L4{false, false, false, true};
  const auto integrals_mnva =
    MBPT::ladder_integrals_mnib(excited, core, valence);
  const auto L1_mnva =
    lowest_order_ladder(qk, excited, core, integrals_mnva, diagram_L1);
  const auto L2_L3_mnva =
    lowest_order_ladder(qk, excited, core, integrals_mnva, diagrams_L2_L3);
  const auto L4_mnva =
    lowest_order_ladder(qk, excited, core, integrals_mnva, diagram_L4);

  // rows[term][valence state]
  const std::vector<std::string> terms{
    "de(2)", "S_d",    "S_x",  "S_xx", "H_d",  "H_x",  "L_cross",
    "L_par", "L_parx", "R_dd", "R_dx", "R_xd", "R_xx", "de(3)"};
  std::vector<std::vector<double>> rows(terms.size());
  std::vector<double> ladder_third_order, ladder_iterated, sigma;

  std::cout << "\nCalculating:" << std::flush;
  for (const auto &v : valence) {
    std::cout << " " << v.shortSymbol() << std::flush;

    const auto de2 = MBPT::de_valence(v, qk, qk, core, excited);

    const auto particle_hole = particle_hole_terms(v, core, excited, qk, sjt);

    // Ladder integrals with the valence orbital in the first index, L_vnab
    const auto integrals_vnab = MBPT::external_m_slot_integrals(
      MBPT::ladder_orbitals(core), MBPT::ladder_orbitals({v}));
    const auto L1_vnab =
      lowest_order_ladder(qk, excited, core, integrals_vnab, diagram_L1);
    const auto L2_L3_vnab =
      lowest_order_ladder(qk, excited, core, integrals_vnab, diagrams_L2_L3);
    const auto L4_vnab =
      lowest_order_ladder(qk, excited, core, integrals_vnab, diagram_L4);

    // {direct, exchange} energy of each ladder diagram
    const auto de_L1_mnva =
      ladder_energy_mnva(v, core, excited, qk, L1_mnva, sjt);
    const auto de_L2_L3_mnva =
      ladder_energy_mnva(v, core, excited, qk, L2_L3_mnva, sjt);
    const auto de_L4_mnva =
      ladder_energy_mnva(v, core, excited, qk, L4_mnva, sjt);
    const auto de_L1_vnab =
      ladder_energy_vnab(v, core, excited, qk, L1_vnab, sjt);
    const auto de_L2_L3_vnab =
      ladder_energy_vnab(v, core, excited, qk, L2_L3_vnab, sjt);
    const auto de_L4_vnab =
      ladder_energy_vnab(v, core, excited, qk, L4_vnab, sjt);

    // L4_mnva and L1_vnab each have a mirror-image diagram, equal in value,
    // that the ladder integrals do not contain: counted twice.
    // (The particle-hole ladders, L2 and L3, are parts of L_cross and H_x)
    const auto L_par =
      de_L1_mnva[0] + 2.0 * de_L4_mnva[0] + 2.0 * de_L1_vnab[0] + de_L4_vnab[0];
    const auto L_parx =
      de_L1_mnva[1] + 2.0 * de_L4_mnva[1] + 2.0 * de_L1_vnab[1] + de_L4_vnab[1];

    const auto relaxation = relaxation_terms(v, core, excited, qk, sjt);

    const std::vector<double> column{de2,
                                     particle_hole.S_d,
                                     particle_hole.S_x,
                                     particle_hole.S_xx,
                                     particle_hole.H_d,
                                     particle_hole.H_x,
                                     particle_hole.L_cross,
                                     L_par,
                                     L_parx,
                                     relaxation.R_dd,
                                     relaxation.R_dx,
                                     relaxation.R_xd,
                                     relaxation.R_xx};
    for (std::size_t i = 0; i < column.size(); ++i) {
      rows[i].push_back(column[i]);
    }
    // de(3): sum of all the third-order terms (all but the first entry)
    rows.back().push_back(
      std::accumulate(column.cbegin() + 1, column.cend(), 0.0));

    // The part of the third order that the ladder integrals contain: all
    // diagrams, direct and exchange
    double de_ladder = 0.0;
    for (const auto &de : {de_L1_mnva, de_L2_L3_mnva, de_L4_mnva, de_L1_vnab,
                           de_L2_L3_vnab, de_L4_vnab}) {
      de_ladder += de[0] + de[1];
    }
    ladder_third_order.push_back(de_ladder);
    if (have_ladder) {
      ladder_iterated.push_back(MBPT::de_valence(v, qk, lk, core, excited));
    }
    if (wf.Sigma() != nullptr) {
      sigma.push_back(v * wf.Sigma()->SigmaFv(v));
    }
  }
  std::cout << "\n";

  const auto print_row = [](const std::string &name,
                            const std::vector<double> &row) {
    fmt::print("{:9s}", name);
    for (const auto value : row) {
      fmt::print(" {:10.2f}", value * PhysConst::Hartree_invcm);
    }
    std::cout << "\n";
  };

  std::cout << "\nThird-order corrections, by class (cm^-1):\n";
  fmt::print("{:9s}", "");
  for (const auto &v : valence) {
    fmt::print(" {:>10s}", v.shortSymbol());
  }
  std::cout << "\n";
  for (std::size_t i = 0; i < terms.size(); ++i) {
    print_row(terms[i], rows[i]);
  }

  std::cout << "\nDiagrams of Ladder{}, with include_L4 (cm^-1):\n";
  print_row("3rd order", ladder_third_order);
  if (have_ladder) {
    print_row("iterated", ladder_iterated);
  }
  if (wf.Sigma() != nullptr) {
    std::cout << "\nCorrelation potential, <v|Sigma|v> (cm^-1):\n";
    print_row("Sigma", sigma);
  }
}

//==============================================================================
namespace ThirdOrder {

std::array<double, 2>
direct_exchange(int k, const DiracSpinor &a, const DiracSpinor &b,
                const DiracSpinor &c, const DiracSpinor &d,
                const Coulomb::QkTable &qk, const Angular::SixJTable &sjt) {
  return {qk.Q(k, a, b, c, d), qk.P(k, a, b, c, d, &sjt)};
}

//==============================================================================

ParticleHoleTerms particle_hole_terms(const DiracSpinor &v,
                                      const std::vector<DiracSpinor> &core,
                                      const std::vector<DiracSpinor> &excited,
                                      const Coulomb::QkTable &qk,
                                      const Angular::SixJTable &sjt) {

  // sum[number of valence vertices that are exchange][core vertex is exchange]
  std::array<std::array<double, 2>, 3> sum{};

  const auto max_k =
    (DiracSpinor::max_tj(core) + DiracSpinor::max_tj(excited)) / 2;
  for (int k = 0; k <= max_k; ++k) {

    // Particle-hole pairs (a, n) that can couple to k: {index of a, of n}
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    for (std::size_t ia = 0; ia < core.size(); ++ia) {
      for (std::size_t in = 0; in < excited.size(); ++in) {
        if (Angular::triangle(core[ia].twoj(), excited[in].twoj(), 2 * k)) {
          pairs.emplace_back(ia, in);
        }
      }
    }
    const auto num_pairs = pairs.size();
    if (num_pairs == 0)
      continue;

    // Core vertices, {direct, exchange}. Symmetric matrices
    std::array<LinAlg::Matrix<double>, 2> G{
      LinAlg::Matrix<double>(num_pairs, num_pairs),
      LinAlg::Matrix<double>(num_pairs, num_pairs)};
    auto T = G;
#pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < num_pairs; ++i) {
      const auto &a = core[pairs[i].first];
      const auto &n = excited[pairs[i].second];
      for (std::size_t j = i; j < num_pairs; ++j) {
        const auto &b = core[pairs[j].first];
        const auto &m = excited[pairs[j].second];
        const auto W_bnma = direct_exchange(k, b, n, m, a, qk, sjt);
        const auto W_abnm = direct_exchange(k, a, b, n, m, qk, sjt);
        const auto de = a.en() + b.en() - n.en() - m.en();
        for (const auto part : {0ul, 1ul}) {
          G[part](i, j) = W_bnma[part];
          G[part](j, i) = W_bnma[part];
          T[part](i, j) = W_abnm[part] / de;
          T[part](j, i) = W_abnm[part] / de;
        }
      }
    }

    // The valence electron is scattered to x: excited or core
    for (const bool x_in_core : {false, true}) {
      const auto &x_orbitals = x_in_core ? core : excited;
      const auto num_x = x_orbitals.size();

      // Valence vertices, {direct, exchange}: one column for each x
      std::array<LinAlg::Matrix<double>, 2> t{
        LinAlg::Matrix<double>(num_pairs, num_x),
        LinAlg::Matrix<double>(num_pairs, num_x)};
      auto u = t;
#pragma omp parallel for schedule(dynamic)
      for (std::size_t ix = 0; ix < num_x; ++ix) {
        const auto &x = x_orbitals[ix];
        if (Angular::triangle(v.twoj(), x.twoj(), 2 * k) == 0)
          continue;
        for (std::size_t i = 0; i < num_pairs; ++i) {
          const auto &a = core[pairs[i].first];
          const auto &n = excited[pairs[i].second];
          const auto W_vaxn = direct_exchange(k, v, a, x, n, qk, sjt);
          const auto W_vnxa = direct_exchange(k, v, n, x, a, qk, sjt);
          const auto de = x_in_core ? x.en() + a.en() - v.en() - n.en() :
                                      v.en() + a.en() - x.en() - n.en();
          for (const auto part : {0ul, 1ul}) {
            t[part](i, ix) = (x_in_core ? W_vnxa[part] : W_vaxn[part]) / de;
            u[part](i, ix) = x_in_core ? W_vaxn[part] : W_vnxa[part];
          }
        }
      }

      // t^T (G t + 2 T u), summed over x, for each choice of direct or
      // exchange in the left valence vertex, core vertex, right valence vertex
      const auto factor = (x_in_core ? -1.0 : 1.0) * Angular::neg1pow(k) /
                          ((2 * k + 1) * (2 * k + 1) * v.twojp1());
      for (const auto core_vertex : {0ul, 1ul}) {
        for (const auto right : {0ul, 1ul}) {
          const auto Gt_2Tu =
            G[core_vertex] * t[right] + 2.0 * (T[core_vertex] * u[right]);
          for (const auto left : {0ul, 1ul}) {
            const auto chain = std::inner_product(
              t[left].cbegin(), t[left].cend(), Gt_2Tu.cbegin(), 0.0);
            sum[left + right][core_vertex] += factor * chain;
          }
        }
      }
    }
  }

  return {sum[0][0], sum[1][0], sum[2][0], sum[0][1], sum[1][1], sum[2][1]};
}

//==============================================================================

Coulomb::LkTable lowest_order_ladder(const Coulomb::QkTable &qk,
                                     const std::vector<DiracSpinor> &excited,
                                     const std::vector<DiracSpinor> &core,
                                     const MBPT::LadderIntegralSet &integrals,
                                     const MBPT::LadderDiagrams &diagrams) {
  const auto excited_orbitals = MBPT::ladder_orbitals(excited);
  const auto core_orbitals = MBPT::ladder_orbitals(core);
  Coulomb::LkTable lk;
  MBPT::add_Lk_entries(&lk, excited_orbitals, integrals, qk.max_k());
  MBPT::calculate_Lk(&lk, qk, excited_orbitals, core_orbitals, integrals,
                     diagrams, nullptr, nullptr, qk.max_k(), false, 0.0);
  return lk;
}

//==============================================================================

std::array<double, 2>
ladder_energy_mnva(const DiracSpinor &v, const std::vector<DiracSpinor> &core,
                   const std::vector<DiracSpinor> &excited,
                   const Coulomb::QkTable &qk, const Coulomb::LkTable &lk,
                   const Angular::SixJTable &sjt) {
  double direct = 0.0;
  double exchange = 0.0;
#pragma omp parallel for reduction(+ : direct, exchange)
  for (std::size_t in = 0; in < excited.size(); ++in) {
    const auto &n = excited[in];
    for (const auto &m : excited) {
      for (const auto &a : core) {
        const auto de = v.en() + a.en() - m.en() - n.en();
        const auto [k0, kI] = Coulomb::k_minmax_W(v, a, m, n);
        for (int k = k0; k <= kI; ++k) {
          const auto L_mnva = lk.Q(k, m, n, v, a);
          if (L_mnva == 0.0)
            continue;
          const auto F_vamn = direct_exchange(k, v, a, m, n, qk, sjt);
          direct += F_vamn[0] * L_mnva / ((2 * k + 1) * de);
          exchange += F_vamn[1] * L_mnva / ((2 * k + 1) * de);
        }
      }
    }
  }
  return {direct / v.twojp1(), exchange / v.twojp1()};
}

//==============================================================================

std::array<double, 2>
ladder_energy_vnab(const DiracSpinor &v, const std::vector<DiracSpinor> &core,
                   const std::vector<DiracSpinor> &excited,
                   const Coulomb::QkTable &qk, const Coulomb::LkTable &lk,
                   const Angular::SixJTable &sjt) {
  double direct = 0.0;
  double exchange = 0.0;
#pragma omp parallel for reduction(+ : direct, exchange)
  for (std::size_t in = 0; in < excited.size(); ++in) {
    const auto &n = excited[in];
    for (const auto &a : core) {
      for (const auto &b : core) {
        const auto de = v.en() + n.en() - a.en() - b.en();
        const auto [k0, kI] = Coulomb::k_minmax_W(v, n, a, b);
        for (int k = k0; k <= kI; ++k) {
          const auto L_vnab = lk.Q(k, v, n, a, b);
          if (L_vnab == 0.0)
            continue;
          const auto F_vnab = direct_exchange(k, v, n, a, b, qk, sjt);
          direct += F_vnab[0] * L_vnab / ((2 * k + 1) * de);
          exchange += F_vnab[1] * L_vnab / ((2 * k + 1) * de);
        }
      }
    }
  }
  return {direct / v.twojp1(), exchange / v.twojp1()};
}

//==============================================================================

RelaxationTerms relaxation_terms(const DiracSpinor &v,
                                 const std::vector<DiracSpinor> &core,
                                 const std::vector<DiracSpinor> &excited,
                                 const Coulomb::QkTable &qk,
                                 const Angular::SixJTable &sjt) {

  double R_dd = 0.0, R_dx = 0.0, R_xd = 0.0, R_xx = 0.0;

  // In each block: rho = {direct, exchange} change of the density, then its
  // interaction with the valence electron

  // Core-excited: rho_ma (and the equal rho_am: factor of two)
#pragma omp parallel for reduction(+ : R_dd, R_dx, R_xd, R_xx) schedule(dynamic)
  for (std::size_t im = 0; im < excited.size(); ++im) {
    const auto &m = excited[im];
    for (const auto &a : core) {
      if (a.kappa() != m.kappa())
        continue;
      std::array<double, 2> rho{0.0, 0.0};
      for (const auto &b : core) {
        for (const auto &n : excited) {
          for (const auto &r : excited) {
            const auto de = a.en() + b.en() - n.en() - r.en();
            const auto [k0, kI] = Coulomb::k_minmax_Q(m, b, n, r);
            for (int k = k0; k <= kI; k += 2) {
              const auto Q_mbnr = qk.Q(k, m, b, n, r);
              if (Q_mbnr == 0.0)
                continue;
              const auto W_nrab = direct_exchange(k, n, r, a, b, qk, sjt);
              rho[0] += Q_mbnr * W_nrab[0] / ((2 * k + 1) * de);
              rho[1] += Q_mbnr * W_nrab[1] / ((2 * k + 1) * de);
            }
          }
          for (const auto &c : core) {
            const auto de = b.en() + c.en() - m.en() - n.en();
            const auto [k0, kI] = Coulomb::k_minmax_Q(b, c, a, n);
            for (int k = k0; k <= kI; k += 2) {
              const auto Q_bcan = qk.Q(k, b, c, a, n);
              if (Q_bcan == 0.0)
                continue;
              const auto W_mnbc = direct_exchange(k, m, n, b, c, qk, sjt);
              rho[0] -= Q_bcan * W_mnbc[0] / ((2 * k + 1) * de);
              rho[1] -= Q_bcan * W_mnbc[1] / ((2 * k + 1) * de);
            }
          }
        }
      }
      const auto W_vavm = direct_exchange(0, v, a, v, m, qk, sjt);
      const auto factor =
        2.0 * Angular::neg1pow_2(v.twoj() + a.twoj() + 2) /
        (std::sqrt(v.twojp1() * a.twojp1()) * (a.en() - m.en()));
      R_dd += factor * W_vavm[0] * rho[0];
      R_dx += factor * W_vavm[0] * rho[1];
      R_xd += factor * W_vavm[1] * rho[0];
      R_xx += factor * W_vavm[1] * rho[1];
    }
  }

  // Core-core: rho_bc (here, the part of it from each excited m)
#pragma omp parallel for reduction(+ : R_dd, R_dx, R_xd, R_xx) schedule(dynamic)
  for (std::size_t im = 0; im < excited.size(); ++im) {
    const auto &m = excited[im];
    for (const auto &b : core) {
      for (const auto &c : core) {
        if (c.kappa() != b.kappa())
          continue;
        std::array<double, 2> rho{0.0, 0.0};
        for (const auto &a : core) {
          for (const auto &n : excited) {
            const auto de = (a.en() + b.en() - m.en() - n.en()) *
                            (a.en() + c.en() - m.en() - n.en());
            const auto [k0, kI] = Coulomb::k_minmax_Q(a, b, m, n);
            for (int k = k0; k <= kI; k += 2) {
              const auto Q_abmn = qk.Q(k, a, b, m, n);
              if (Q_abmn == 0.0)
                continue;
              const auto W_mnac = direct_exchange(k, m, n, a, c, qk, sjt);
              rho[0] -= Q_abmn * W_mnac[0] / ((2 * k + 1) * de);
              rho[1] -= Q_abmn * W_mnac[1] / ((2 * k + 1) * de);
            }
          }
        }
        const auto W_vcvb = direct_exchange(0, v, c, v, b, qk, sjt);
        const auto factor = Angular::neg1pow_2(v.twoj() + b.twoj() + 2) /
                            std::sqrt(v.twojp1() * b.twojp1());
        R_dd += factor * W_vcvb[0] * rho[0];
        R_dx += factor * W_vcvb[0] * rho[1];
        R_xd += factor * W_vcvb[1] * rho[0];
        R_xx += factor * W_vcvb[1] * rho[1];
      }
    }
  }

  // Excited-excited: rho_mr
#pragma omp parallel for reduction(+ : R_dd, R_dx, R_xd, R_xx) schedule(dynamic)
  for (std::size_t im = 0; im < excited.size(); ++im) {
    const auto &m = excited[im];
    for (const auto &r : excited) {
      if (r.kappa() != m.kappa())
        continue;
      std::array<double, 2> rho{0.0, 0.0};
      for (const auto &a : core) {
        for (const auto &b : core) {
          for (const auto &n : excited) {
            const auto de = (a.en() + b.en() - m.en() - n.en()) *
                            (a.en() + b.en() - n.en() - r.en());
            const auto [k0, kI] = Coulomb::k_minmax_Q(a, b, m, n);
            for (int k = k0; k <= kI; k += 2) {
              const auto Q_abmn = qk.Q(k, a, b, m, n);
              if (Q_abmn == 0.0)
                continue;
              const auto W_nrba = direct_exchange(k, n, r, b, a, qk, sjt);
              rho[0] += Q_abmn * W_nrba[0] / ((2 * k + 1) * de);
              rho[1] += Q_abmn * W_nrba[1] / ((2 * k + 1) * de);
            }
          }
        }
      }
      const auto W_vmvr = direct_exchange(0, v, m, v, r, qk, sjt);
      const auto factor = Angular::neg1pow_2(v.twoj() + m.twoj() + 2) /
                          std::sqrt(v.twojp1() * m.twojp1());
      R_dd += factor * W_vmvr[0] * rho[0];
      R_dx += factor * W_vmvr[0] * rho[1];
      R_xd += factor * W_vmvr[1] * rho[0];
      R_xx += factor * W_vmvr[1] * rho[1];
    }
  }

  return {R_dd, R_dx, R_xd, R_xx};
}

} // namespace ThirdOrder

} // namespace Module
