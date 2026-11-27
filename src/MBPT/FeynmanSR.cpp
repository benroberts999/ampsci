#include "MBPT/FeynmanSR.hpp"
#include "Angular/Wigner369j.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "ExternalField/CorePolarisation.hpp"
#include "MBPT/Feynman.hpp"
#include "MBPT/RadialMatrix.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/omp.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

namespace MBPT {

//==============================================================================
// Structure Radiation
//==============================================================================

ComplexGMatrix FeynmanSR::perturbed_greens(ComplexGMatrix gA,
                                           ComplexGMatrix gB) const {

  // This the generalised perturbed greens function, which is
  // defined as dG = gA * h * gB, where,
  // At the moment, h is the dipole operator, which is defined as h = r.
  ComplexGMatrix dG(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);
  std::size_t smax = m_include_G ? 2 : 1;

  for (auto k = 0ul; k < m_subgrid_points; ++k) {
    const auto rk_dr = gA.r(k) * gA.dr(k);
    for (auto i = 0ul; i < m_subgrid_points; ++i) {
      for (auto j = 0ul; j < m_subgrid_points; ++j) {
        for (auto a = 0ul; a < smax; a++) {
          for (auto b = 0ul; b < smax; b++) {
            for (auto c = 0ul; c < smax; c++) {
              dG.sp(a, b)(i, j) +=
                gA.sp(a, c)(i, k) * rk_dr * gB.sp(c, b)(k, j);
            }
          }
        }
      }
    }
  }

  return dG;
}

//==============================================================================

double L_ang_dSR(int K, int kw, int kv, int kB, int kA, int k) {
  const auto tj_w = Angular::twoj_k(kw);
  const auto tj_v = Angular::twoj_k(kv);
  const auto tj_B = Angular::twoj_k(kB);
  const auto tj_A = Angular::twoj_k(kA);

  const auto sjs = Angular::sixj_2(tj_w, 2 * K, tj_v, tj_B, 2 * k, tj_A);
  const auto s = Angular::neg1pow_2(tj_A + tj_B + k * 2 + K * 2);
  return s * sjs;
}

//==============================================================================

GMatrix FeynmanSR::Sigma_SR_direct(int kv, double env, int kw, double Omega) {
  // Coded for dipole operator,
  // kappa_v, en, kappa_w, Omega-> external field
  //<v|Sigma(env)|w > = <v|gA pi gB| w >
  // This is the diagram:
  //                   /
  //                  /
  //                 /
  //   -----.-------X-------.-----
  //         \             /
  //          \    .-.    /
  //           '--(   )--'
  //               '-'

  if (!has_qpiq()) {
    form_qpiq(m_ident);
  }
  GMatrix Sigma(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

  constexpr std::complex<double> I{0.0, 1.0};
  const auto num_kappas = m_max_ki + 1;

// Tell OpenMP how to reduce GMatrix
#pragma omp declare reduction(+ : GMatrix : omp_out += omp_in)                 \
  initializer(omp_priv = omp_orig)

  // #pragma omp parallel for collapse(2) reduction(+ : Sigma)
  for (auto iw = 0ul; iw < m_wgrid_points.size(); iw++) {
    for (auto iB = 0ul; iB < num_kappas; ++iB) {
      const auto kB = Angular::kindex_to_kappa(iB);
      for (auto iA = 0ul; iA < num_kappas; ++iA) {
        const auto kA = Angular::kindex_to_kappa(iA);
        const auto omega = std::complex{m_omre, m_wgrid_points[iw]};
        // Simpson's rule: Implicit ends (integrand zero at w=0 and w>wmax)
        // const auto weight = iw % 2 == 0 ? 4.0 / 3 : 2.0 / 3;
        // I, since dw is on imag. grid; 2 from symmetric +/- w
        const auto dw = I * m_wgrid_weights[iw];
        const auto gA = m_include_G ? green(kA, env + omega + Omega) :
                                      green(kA, env + omega + Omega).drop_g();
        const auto gB = m_include_G ? green(kB, env + omega) :
                                      green(kB, env + omega).drop_g();

        const auto dG_AB = perturbed_greens(gA, gB);
        for (auto k = 0ul; int(k) <= m_max_k; k++) {
          // - 1 for dipole constant factor
          const auto K = 1;
          const auto h_AB = -1 * Angular::Ck_kk(K, kA, kB);
          if (h_AB == 0.0)
            continue;

          const auto ck_Bw = Angular::Ck_kk(int(k), kB, kw);
          const auto ck_vA = Angular::Ck_kk(int(k), kv, kA);

          if (ck_Bw == 0.0 or ck_vA == 0.0)
            continue;

          const auto &qpq_dw = m_qpiq_wk[iw][k];

          const auto LAng = L_ang_dSR(K, kv, kw, kB, kA, int(k));

          if (LAng == 0.0)
            continue;

          const auto c_ang_dw = LAng * h_AB * dw * ck_Bw *
                                ck_vA; // / double(Angular::twoj_k(kv) + 1);

          Sigma += (c_ang_dw * mult_elements(dG_AB, qpq_dw)).real();

          // }
        }
      }
    }
  }
  // Extra 2 from symmetric +/-w
  // Sigma *= (m_wgrid_points.du() / (M_PI));
  Sigma *= 1.0 / (M_PI);

  return Sigma;
}

//==============================================================================

ComplexGMatrix
FeynmanSR::perturbed_greens(const DiracOperator::TensorOperator &h, int kappa_A,
                            int kappa_B, const ComplexGMatrix &gA,
                            const ComplexGMatrix &gB,
                            const std::optional<ComplexGMatrix> &dV_AB) const {

  // gA * <A||h + dV||B> * gB: angular factor of h times its radial part
  // (local in r), plus the core polarisation matrix (non-local), if given

  // Radial part of the operator between partial waves A and B, t^{mu nu}(r):
  // its action on unit spinors, (1, 0) and (0, 1). Valid for any operator
  // that is local in r
  DiracSpinor unit_f(0, kappa_B, m_grid);
  DiracSpinor unit_g(0, kappa_B, m_grid);
  std::fill(unit_f.f().begin(), unit_f.f().end(), 1.0);
  std::fill(unit_g.g().begin(), unit_g.g().end(), 1.0);
  // t_f = (t^ff, t^gf), t_g = (t^fg, t^gg)
  const auto h_AB = h.angularF(kappa_A, kappa_B);
  const auto t_f = h_AB * h.radial_rhs(kappa_A, unit_f);
  const auto t_g = h_AB * h.radial_rhs(kappa_A, unit_g);

  // t*gB, with the integration measure: (t gB)^{mu nu}_ij =
  // dr_i * sum_lambda t^{mu lambda}(r_i) * gB^{lambda nu}_ij
  ComplexGMatrix t_gB(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto full_index = gB.index_to_fullgrid(i);
    const auto dr = gB.dr(i);
    const auto t_ff = t_f.f(full_index) * dr;
    const auto t_gf = t_f.g(full_index) * dr;
    const auto t_fg = t_g.f(full_index) * dr;
    const auto t_gg = t_g.g(full_index) * dr;
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      t_gB.ff(i, j) = t_ff * gB.ff(i, j);
      if (m_include_G) {
        t_gB.ff(i, j) += t_fg * gB.gf(i, j);
        t_gB.fg(i, j) = t_ff * gB.fg(i, j) + t_fg * gB.gg(i, j);
        t_gB.gf(i, j) = t_gf * gB.ff(i, j) + t_gg * gB.gf(i, j);
        t_gB.gg(i, j) = t_gf * gB.fg(i, j) + t_gg * gB.gg(i, j);
      }
    }
  }

  if (dV_AB) {
    t_gB += (*dV_AB) * gB;
  }

  return gA * t_gB;
}

//==============================================================================

LinAlg::Matrix<std::optional<ComplexGMatrix>>
FeynmanSR::dV_matrices(const DiracOperator::TensorOperator &h,
                       const ExternalField::CorePolarisation *dV, bool conj) {

  // dV_AB(r1,r2) = sum_n [dV*F_n]_A(r1) F_n^dag(r2), with n a complete set
  // of states of partial wave B, and [..]_A the projection onto partial
  // wave A: the reduced matrix element is <a||dV||b> = F_a * (dV_AB * F_b)

  const auto num_kappas = m_max_ki + 1;
  LinAlg::Matrix<std::optional<ComplexGMatrix>> dV_kappa(num_kappas,
                                                         num_kappas);
  if (dV == nullptr)
    return dV_kappa;

  // The complete set of states is the same for every operator: formed once
  if (m_coordinate_basis.empty()) {
    m_coordinate_basis = coordinate_basis();
  }

#pragma omp parallel for collapse(2) schedule(dynamic)
  for (auto i_A = 0ul; i_A < num_kappas; ++i_A) {
    for (auto i_B = 0ul; i_B < num_kappas; ++i_B) {
      const auto kappa_A = Angular::kindex_to_kappa(i_A);
      const auto kappa_B = Angular::kindex_to_kappa(i_B);
      if (h.isZero(kappa_A, kappa_B))
        continue;

      GMatrix dV_AB(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);
      for (const auto &Fn : m_coordinate_basis) {
        if (Fn.kappa() != kappa_B)
          continue;
        dV_AB.add(dV->dV_rhs(kappa_A, Fn, conj), Fn, 1.0);
      }
      // includes both integration measures
      dV_AB.dri_in_place();
      dV_AB.drj_in_place();
      dV_kappa(i_A, i_B) = dV_AB.complex();
    }
  }

  return dV_kappa;
}

//==============================================================================

GMatrix FeynmanSR::dSigma_a(const DiracOperator::TensorOperator &h, int kappa_w,
                            int kappa_v, double en_v, double omega,
                            const ExternalField::CorePolarisation *dV) {

  // dSigma^(a)_12(e_v) = int de/2pi dG_12(e_v + e) Q_1i Pi_ij(e) Q_j2,
  // dG_12(e) = G_1i(e + omega) h_i G_i2(e)

  if (!has_qpiq()) {
    form_qpiq(m_ident);
  }

  // Core polarisation of the external field vertex. dV is solved at
  // |omega|: for negative omega the field is its conjugate component
  const auto dV_kappa = dV_matrices(h, dV, omega < 0.0);

  const GMatrix zero(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

  constexpr std::complex<double> I{0.0, 1.0};
  const auto num_kappas = m_max_ki + 1;
  const auto K = h.rank();

  // BLAS must run single-threaded inside the omp region below
  const qip::SingleThreadBlas single_thread_blas{};

  std::vector<GMatrix> dSigma_ts(std::size_t(omp_get_max_threads()), zero);

#pragma omp parallel for
  for (auto ie = 0ul; ie < m_wgrid_points.size(); ie++) {

    auto &dSigma_t = dSigma_ts[std::size_t(omp_get_thread_num())];

    // Integration energy, along the imaginary axis
    const auto epsilon = std::complex{m_omre, m_wgrid_points[ie]};
    const auto d_epsilon = I * m_wgrid_weights[ie];

    // Green's functions either side of the external field, each partial wave:
    // G_alpha(e_v + e + omega) h G_beta(e_v + e)
    std::vector<ComplexGMatrix> G_alpha, G_beta;
    G_alpha.reserve(num_kappas);
    G_beta.reserve(num_kappas);
    for (auto ik = 0ul; ik < num_kappas; ++ik) {
      const auto kappa = Angular::kindex_to_kappa(ik);
      G_alpha.push_back(green(kappa, en_v + epsilon + omega));
      G_beta.push_back(green(kappa, en_v + epsilon));
    }

    for (auto i_beta = 0ul; i_beta < num_kappas; ++i_beta) {
      const auto kappa_beta = Angular::kindex_to_kappa(i_beta);
      for (auto i_alpha = 0ul; i_alpha < num_kappas; ++i_alpha) {
        const auto kappa_alpha = Angular::kindex_to_kappa(i_alpha);

        if (h.isZero(kappa_alpha, kappa_beta))
          continue;

        // Q*Pi*Q, summed over multipolarity k with the angular factors
        ComplexRMatrix qpiq(m_i0, m_stride, m_subgrid_points, m_grid);
        bool non_zero = false;
        for (auto k = 0ul; int(k) <= m_max_k; k++) {

          const auto ck_w_alpha = Angular::Ck_kk(int(k), kappa_w, kappa_alpha);
          const auto ck_beta_v = Angular::Ck_kk(int(k), kappa_beta, kappa_v);
          if (ck_w_alpha == 0.0 || ck_beta_v == 0.0)
            continue;

          const auto L =
            L_ang_dSR(K, kappa_w, kappa_v, kappa_beta, kappa_alpha, int(k));
          if (L == 0.0)
            continue;

          const std::complex<double> c_ang{L * ck_w_alpha * ck_beta_v, 0.0};
          qpiq += c_ang * m_qpiq_wk[ie][k];
          non_zero = true;
        }
        if (!non_zero)
          continue;

        const auto dG =
          perturbed_greens(h, kappa_alpha, kappa_beta, G_alpha[i_alpha],
                           G_beta[i_beta], dV_kappa(i_alpha, i_beta));

        dSigma_t += (d_epsilon * mult_elements(dG, qpiq)).real();
      }
    }
  }

  GMatrix dSigma = zero;
  for (const auto &dSigma_t : dSigma_ts) {
    dSigma += dSigma_t;
  }

  // 1/pi = (1/2pi)*2; the 2 from the symmetric +/- Im(e) contour halves
  // (integrand at -Im(e) is conjugate of that at +Im(e), and we take Re part)
  dSigma *= (1.0 / M_PI);

  return dSigma;
}

//==============================================================================

ComplexRMatrix FeynmanSR::dQ_screening(int k, std::size_t ie) const {
  if (!m_screen_Coulomb) {
    return ComplexRMatrix(m_i0, m_stride, m_subgrid_points, m_grid);
  }
  // q~ = q - i q*pi*q~: the stored Q*Pi*Q is formed with the screened line
  constexpr std::complex<double> I{0.0, 1.0};
  return (-I) * m_qpiq_wk[ie][std::size_t(k)];
}

//==============================================================================

ComplexGMatrix FeynmanSR::green_loop(int kappa, std::complex<double> en,
                                     const DiracSpinor *Fc_hp) const {
  // G^e + G^c, with G^c over the polarised core states only
  auto G_loop = green_excited(kappa, en, Fc_hp);
  const auto &core = m_HF->core();
  for (const auto ia : m_core_index) {
    if (core[ia].kappa() != kappa)
      continue;
    const auto inv_de = 1.0 / (en - std::complex<double>{core[ia].en()});
    G_loop += inv_de * m_pa[ia];
  }
  return G_loop;
}

//==============================================================================

ComplexRMatrix trace_product(const ComplexGMatrix &A, const ComplexGMatrix &X) {
  // Same sub-grid as A
  auto trace = A.radial(0, 0);
  const auto num_spins = std::min(A.num_spins(), X.num_spins());
  for (auto i = 0ul; i < A.size(); ++i) {
    for (auto j = 0ul; j < A.size(); ++j) {
      std::complex<double> sum{0.0, 0.0};
      for (auto mu = 0ul; mu < num_spins; ++mu) {
        for (auto nu = 0ul; nu < num_spins; ++nu) {
          sum += A.sp(mu, nu)(i, j) * X.sp(nu, mu)(j, i);
        }
      }
      trace(i, j) = sum;
    }
  }
  return trace;
}

//==============================================================================

double L_ang_dPi(int K, int k, int l, int kappa_A, int kappa_B, int kappa_C) {
  const auto tj_A = Angular::twoj_k(kappa_A);
  const auto tj_B = Angular::twoj_k(kappa_B);
  const auto tj_C = Angular::twoj_k(kappa_C);
  const auto sjs = Angular::sixj_2(2 * k, 2 * l, 2 * K, tj_B, tj_C, tj_A);
  if (sjs == 0.0)
    return 0.0;
  return Angular::neg1pow(k + l + K + 1) * sjs *
         Angular::Ck_kk(k, kappa_C, kappa_A) *
         Angular::Ck_kk(l, kappa_A, kappa_B);
}

//==============================================================================

double L_ang_dSigma_bc(int K, int kappa_w, int kappa_v, int kappa_alpha, int k,
                       int l) {
  const auto tj_w = Angular::twoj_k(kappa_w);
  const auto tj_v = Angular::twoj_k(kappa_v);
  const auto tj_alpha = Angular::twoj_k(kappa_alpha);
  const auto sjs = Angular::sixj_2(tj_w, 2 * K, tj_v, 2 * l, tj_alpha, 2 * k);
  return Angular::neg1pow_2(tj_w + tj_v) * sjs;
}

//==============================================================================

LinAlg::Matrix<ComplexRMatrix> FeynmanSR::dPi_b(
  const DiracOperator::TensorOperator &h, std::complex<double> epsilon,
  double omega,
  const LinAlg::Matrix<std::optional<ComplexGMatrix>> &dV_kappa) const {

  // dPi^(b)_12(e) = Tr int de'/2pi dG_21(e' + e) G_12(e' + omega)
  //
  // The integral is the sum of the residues at the core poles of each of
  // the three Green's functions. Each residue is a term of the form
  // Tr[ A_12 (B h C)_21 ], with one of A, B, C the core projector P^a, and
  // the other two the loop Green's functions, G = G^e + G^c:
  //   P^a_12                  [G(e_a + e) h G(e_a + e - omega)]_21
  //   G_12(e_a - e)           [P^a h G(e_a - omega)]_21
  //   G_12(e_a + omega - e)   [G(e_a + omega) h P^a]_21

  const auto num_ks = std::size_t(m_max_k + 1);
  const auto num_kappas = m_max_ki + 1;
  const auto K = h.rank();

  LinAlg::Matrix<ComplexRMatrix> dPi_kl(
    num_ks, num_ks, ComplexRMatrix{m_i0, m_stride, m_subgrid_points, m_grid});

  // Adds one term, Tr[ A_12 (B h C)_21 ], to each dPi^{kl}, with its
  // angular factor. BhC = B*<B||h||C>*C (see perturbed_greens)
  const auto add_term = [&](int kappa_A, int kappa_B, int kappa_C,
                            const ComplexGMatrix &A,
                            const ComplexGMatrix &BhC) {
    std::vector<std::tuple<std::size_t, std::size_t, double>> kl_coefficients;
    for (auto k = 0ul; k < num_ks; ++k) {
      for (auto l = 0ul; l < num_ks; ++l) {
        const auto L = L_ang_dPi(K, int(k), int(l), kappa_A, kappa_B, kappa_C);
        if (L != 0.0) {
          kl_coefficients.emplace_back(k, l, L);
        }
      }
    }
    if (kl_coefficients.empty())
      return;
    const auto trace = trace_product(A, BhC);
    for (const auto &[k, l, coefficient] : kl_coefficients) {
      dPi_kl[k][l] += std::complex<double>{coefficient} * trace;
    }
  };

  const auto &core = m_HF->core();
  for (const auto ia : m_core_index) {
    const auto &Fa = core[ia];
    const auto kappa_a = Fa.kappa();
    const auto en_a = std::complex<double>{Fa.en()};
    const auto i_a = std::size_t(Fa.k_index());
    const auto &Pa = m_pa[ia];
    // Hole-particle: the excited electron feels the hole in a
    const auto *Fa_hp = m_hole_particle ? &Fa : nullptr;

    // Loop Green's functions, for each partial wave
    std::vector<ComplexGMatrix> G_plus, G_plus_omega, G_minus, G_minus_omega;
    for (auto ik = 0ul; ik < num_kappas; ++ik) {
      const auto kappa = Angular::kindex_to_kappa(ik);
      G_plus.push_back(green_loop(kappa, en_a + epsilon, Fa_hp));
      G_plus_omega.push_back(green_loop(kappa, en_a + epsilon - omega, Fa_hp));
      G_minus.push_back(green_loop(kappa, en_a - epsilon, Fa_hp));
      G_minus_omega.push_back(green_loop(kappa, en_a + omega - epsilon, Fa_hp));
    }

    // P^a_12 [G(e_a + e) h G(e_a + e - omega)]_21
    for (auto i_beta = 0ul; i_beta < num_kappas; ++i_beta) {
      const auto kappa_beta = Angular::kindex_to_kappa(i_beta);
      for (auto i_gamma = 0ul; i_gamma < num_kappas; ++i_gamma) {
        const auto kappa_gamma = Angular::kindex_to_kappa(i_gamma);
        if (h.isZero(kappa_beta, kappa_gamma))
          continue;
        const auto GhG =
          perturbed_greens(h, kappa_beta, kappa_gamma, G_plus[i_beta],
                           G_plus_omega[i_gamma], dV_kappa(i_beta, i_gamma));
        add_term(kappa_a, kappa_beta, kappa_gamma, Pa, GhG);
      }
    }

    for (auto i_gamma = 0ul; i_gamma < num_kappas; ++i_gamma) {
      const auto kappa_gamma = Angular::kindex_to_kappa(i_gamma);

      // G_12(e_a - e) [P^a h G(e_a - omega)]_21
      if (!h.isZero(kappa_a, kappa_gamma)) {
        const auto PhG = perturbed_greens(
          h, kappa_a, kappa_gamma, Pa,
          green_loop(kappa_gamma, en_a - omega, Fa_hp), dV_kappa(i_a, i_gamma));
        for (auto i_beta = 0ul; i_beta < num_kappas; ++i_beta) {
          const auto kappa_beta = Angular::kindex_to_kappa(i_beta);
          add_term(kappa_beta, kappa_a, kappa_gamma, G_minus[i_beta], PhG);
        }
      }

      // G_12(e_a + omega - e) [G(e_a + omega) h P^a]_21
      if (!h.isZero(kappa_gamma, kappa_a)) {
        const auto GhP = perturbed_greens(
          h, kappa_gamma, kappa_a, green_loop(kappa_gamma, en_a + omega, Fa_hp),
          Pa, dV_kappa(i_gamma, i_a));
        for (auto i_beta = 0ul; i_beta < num_kappas; ++i_beta) {
          const auto kappa_beta = Angular::kindex_to_kappa(i_beta);
          add_term(kappa_beta, kappa_gamma, kappa_a, G_minus_omega[i_beta],
                   GhP);
        }
      }
    }
  }

  const auto Iunit = std::complex<double>{0.0, 1.0};
  for (auto &dPi : dPi_kl) {
    dPi *= Iunit;
  }
  return dPi_kl;
}

//==============================================================================

std::pair<GMatrix, GMatrix>
FeynmanSR::dSigma_bc(const DiracOperator::TensorOperator &h, int kappa_w,
                     int kappa_v, double en_v, double omega,
                     const ExternalField::CorePolarisation *dV) {

  // dSigma^(b+c)_12(e_v) = int de/2pi G_12(e_v + e) Q_1i dPi_ij(e) Q_j2,
  // dPi = dPi^(b) + dPi^(c), with dPi^(c)_12(e) = dPi^(b)_21(omega - e)
  //
  // With screening, both Coulomb lines are the screened ones, at the
  // integration energy e. (The line next to w is at e - omega: the
  // difference is neglected; none for omega = 0.)

  if (m_screen_Coulomb && !has_qpiq()) {
    form_qpiq(m_ident);
  }

  // Core polarisation of the external field vertex. dV is solved at
  // |omega|: for negative omega the field is its conjugate component
  const auto dV_kappa = dV_matrices(h, dV, omega < 0.0);

  const GMatrix zero(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

  constexpr std::complex<double> I{0.0, 1.0};
  const auto num_kappas = m_max_ki + 1;
  const auto num_ks = std::size_t(m_max_k + 1);
  const auto K = h.rank();

  // BLAS must run single-threaded inside the omp region below
  const qip::SingleThreadBlas single_thread_blas{};

  const auto num_threads = std::size_t(omp_get_max_threads());
  std::vector<GMatrix> dSigma_b_ts(num_threads, zero);
  std::vector<GMatrix> dSigma_c_ts(num_threads, zero);

#pragma omp parallel for
  for (auto ie = 0ul; ie < m_wgrid_points.size(); ie++) {

    auto &dSigma_b_t = dSigma_b_ts[std::size_t(omp_get_thread_num())];
    auto &dSigma_c_t = dSigma_c_ts[std::size_t(omp_get_thread_num())];

    // Integration energy, along the imaginary axis
    const auto epsilon = std::complex{m_omre, m_wgrid_points[ie]};
    const auto d_epsilon = I * m_wgrid_weights[ie];

    // dPi^(b)(e), and dPi^(b)(omega - e) for dPi^(c)
    const auto dPi_b_kl = dPi_b(h, epsilon, omega, dV_kappa);
    const auto dPi_b_swapped_kl = dPi_b(h, omega - epsilon, omega, dV_kappa);

    // Outer Green's function, each partial wave
    std::vector<ComplexGMatrix> G_alpha;
    G_alpha.reserve(num_kappas);
    for (auto i_alpha = 0ul; i_alpha < num_kappas; ++i_alpha) {
      G_alpha.push_back(
        green(Angular::kindex_to_kappa(i_alpha), en_v + epsilon));
    }

    for (auto k = 0ul; k < num_ks; ++k) {
      for (auto l = 0ul; l < num_ks; ++l) {
        if (Angular::triangle(int(k), int(l), K) == 0)
          continue;

        // dPi^(c)_kl(r1, r2; e) = (-1)^(k+l+K) dPi^(b)_lk(r2, r1; omega - e)
        const std::complex<double> sign_c{
          double(Angular::neg1pow(int(k) + int(l) + K)), 0.0};
        const auto dPi_c_kl = sign_c * dPi_b_swapped_kl[l][k].transpose();

        // Screened Coulomb lines, Q + dQ
        const auto qk_screened = get_qk(int(k)) + dQ_screening(int(k), ie);
        const auto ql_screened = get_qk(int(l)) + dQ_screening(int(l), ie);

        // Coulomb lines: q^k * dPi_kl * q^l (with integration measures)
        const auto ql_dri = ql_screened.dri();
        const auto q_dPi_b_q = qk_screened * dPi_b_kl[k][l] * ql_dri;
        const auto q_dPi_c_q = qk_screened * dPi_c_kl * ql_dri;

        for (auto i_alpha = 0ul; i_alpha < num_kappas; ++i_alpha) {
          const auto kappa_alpha = Angular::kindex_to_kappa(i_alpha);

          const auto ck_w_alpha = Angular::Ck_kk(int(k), kappa_w, kappa_alpha);
          const auto cl_alpha_v = Angular::Ck_kk(int(l), kappa_alpha, kappa_v);
          if (ck_w_alpha == 0.0 || cl_alpha_v == 0.0)
            continue;

          const auto L =
            L_ang_dSigma_bc(K, kappa_w, kappa_v, kappa_alpha, int(k), int(l));
          if (L == 0.0)
            continue;

          const auto c_ang_de = d_epsilon * (L * ck_w_alpha * cl_alpha_v);
          dSigma_b_t +=
            (c_ang_de * mult_elements(G_alpha[i_alpha], q_dPi_b_q)).real();
          dSigma_c_t +=
            (c_ang_de * mult_elements(G_alpha[i_alpha], q_dPi_c_q)).real();
        }
      }
    }
  }

  GMatrix dSigma_b = zero;
  GMatrix dSigma_c = zero;
  for (auto i = 0ul; i < num_threads; ++i) {
    dSigma_b += dSigma_b_ts[i];
    dSigma_c += dSigma_c_ts[i];
  }

  // 1/pi = (1/2pi)*2; the 2 from the symmetric +/- Im(e) contour halves
  // (integrand at -Im(e) is conjugate of that at +Im(e), and we take Re part)
  dSigma_b *= (1.0 / M_PI);
  dSigma_c *= (1.0 / M_PI);

  return {dSigma_b, dSigma_c};
}

//==============================================================================

double L_ang_dSigma_ax(int kappa_v, int kappa_beta, int kappa_gamma,
                       int kappa_delta, int k, int l) {
  const auto cl_beta_gamma = Angular::Ck_kk(l, kappa_beta, kappa_gamma);
  const auto ck_gamma_delta = Angular::Ck_kk(k, kappa_gamma, kappa_delta);
  const auto cl_delta_v = Angular::Ck_kk(l, kappa_delta, kappa_v);
  if (cl_beta_gamma == 0.0 || ck_gamma_delta == 0.0 || cl_delta_v == 0.0)
    return 0.0;
  const auto tj_v = Angular::twoj_k(kappa_v);
  const auto tj_beta = Angular::twoj_k(kappa_beta);
  const auto tj_gamma = Angular::twoj_k(kappa_gamma);
  const auto tj_delta = Angular::twoj_k(kappa_delta);
  const auto sjs =
    Angular::sixj_2(tj_beta, 2 * l, tj_gamma, tj_delta, 2 * k, tj_v);
  return Angular::neg1pow_2(tj_gamma + tj_delta + 2 * k + 2 * l) * sjs *
         cl_beta_gamma * ck_gamma_delta * cl_delta_v;
}

//==============================================================================

LinAlg::Matrix<std::optional<ComplexGMatrix>>
FeynmanSR::perturbed_greens_summed(
  const DiracOperator::TensorOperator &h, int kappa_w, int kappa_v,
  const std::vector<ComplexGMatrix> &G_alpha,
  const std::vector<ComplexGMatrix> &G_beta,
  const LinAlg::Matrix<std::optional<ComplexGMatrix>> &dV_kappa) const {

  // dG_beta^k = sum_alpha L <alpha||h||beta> C^k_{w alpha} * G_alpha h G_beta

  const auto num_kappas = m_max_ki + 1;
  const auto num_ks = std::size_t(m_max_k + 1);
  const auto K = h.rank();

  LinAlg::Matrix<std::optional<ComplexGMatrix>> dG_beta_k(num_kappas, num_ks);

  for (auto i_beta = 0ul; i_beta < num_kappas; ++i_beta) {
    const auto kappa_beta = Angular::kindex_to_kappa(i_beta);
    for (auto i_alpha = 0ul; i_alpha < num_kappas; ++i_alpha) {
      const auto kappa_alpha = Angular::kindex_to_kappa(i_alpha);

      if (h.isZero(kappa_alpha, kappa_beta))
        continue;

      // Formed when first needed: most k vanish by selection rules
      std::optional<ComplexGMatrix> dG;
      for (auto k = 0ul; k < num_ks; ++k) {
        const auto ck_w_alpha = Angular::Ck_kk(int(k), kappa_w, kappa_alpha);
        const auto L =
          L_ang_dSR(K, kappa_w, kappa_v, kappa_beta, kappa_alpha, int(k));
        const std::complex<double> c_ang{L * ck_w_alpha, 0.0};
        if (c_ang == 0.0)
          continue;
        if (!dG) {
          dG = perturbed_greens(h, kappa_alpha, kappa_beta, G_alpha[i_alpha],
                                G_beta[i_beta], dV_kappa(i_alpha, i_beta));
        }
        auto &dG_k = dG_beta_k(i_beta, k);
        if (dG_k) {
          *dG_k += c_ang * (*dG);
        } else {
          dG_k = c_ang * (*dG);
        }
      }
    }
  }

  return dG_beta_k;
}

//==============================================================================

GMatrix FeynmanSR::dSigma_ax_core_term(
  const LinAlg::Matrix<std::optional<ComplexGMatrix>> &dG_beta_k, int kappa_v,
  const DiracSpinor &Fa, const std::vector<ComplexGMatrix> &G_delta,
  const std::vector<ComplexGMatrix> &G_gamma, std::size_t ie,
  std::complex<double> i_d_epsilon) const {

  // Re[ i*d_epsilon * dG_1i Q_1j Gamma^a_ij2 Q_i2 ],
  // Gamma^a_ij2 = P^a_ij G^e_j2(e_a - e) + G^e_ij(e_a + e) P^a_j2
  //
  // The two lines of Gamma are gamma (i to j) and delta (j to 2): gamma = a
  // in the first term, delta = a in the second.
  // P^a_12 = F_a(r_1) F_a(r_2) is rank one, so each term is a product of
  // matrices in two coordinates.
  //
  // With screening, Q_1j (multipolarity k) is the screened line: it is at
  // the integration energy e. Q_i2 is at the energy of the integral inside
  // Gamma, and stays bare.

  const ComplexRMatrix radial_zero(m_i0, m_stride, m_subgrid_points, m_grid);
  const auto num_kappas = m_max_ki + 1;
  const auto num_ks = std::size_t(m_max_k + 1);
  // Spinor components in the contractions: f only, unless g is included
  // in the exchange diagrams
  const auto num_spinor = m_include_G_exchange ? num_spins() : 1ul;
  const auto kappa_a = Fa.kappa();

  GMatrix dSigma(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

  // First term of Gamma, P^a_ij G^e_j2(e_a - e) (gamma = a):
  // dSigma_12 = [dG F_a q^l]_12 * [q^k F_a G_delta]_12
  for (auto k = 0ul; k < num_ks; ++k) {
    for (auto i_delta = 0ul; i_delta < num_kappas; ++i_delta) {
      const auto kappa_delta = Angular::kindex_to_kappa(i_delta);
      if (Angular::Ck_kk(int(k), kappa_a, kappa_delta) == 0.0)
        continue;

      // [q^k F_a G_delta]^nu_12
      //   = sum_t q^k_1j F_a^t(r_j) G_delta^{t nu}_j2
      // Formed when first needed
      std::vector<ComplexRMatrix> qk_Fa_G;

      for (auto l = 0ul; l < num_ks; ++l) {

        // [dG F_a]^mu_1i = sum_beta L sum_s dG_beta^{mu s}_1i F_a^s(r_i)
        std::vector<ComplexRMatrix> dG_Fa(num_spinor, radial_zero);
        bool non_zero = false;
        for (auto i_beta = 0ul; i_beta < num_kappas; ++i_beta) {
          const auto &dG = dG_beta_k(i_beta, k);
          if (!dG)
            continue;
          const auto L =
            L_ang_dSigma_ax(kappa_v, Angular::kindex_to_kappa(i_beta), kappa_a,
                            kappa_delta, int(k), int(l));
          if (L == 0.0)
            continue;
          for (auto mu = 0ul; mu < num_spinor; ++mu) {
            for (auto s = 0ul; s < num_spinor; ++s) {
              dG_Fa[mu] +=
                L * mult_cols_full(dG->radial(mu, s), Fa.component(s));
            }
          }
          non_zero = true;
        }
        if (!non_zero)
          continue;

        if (qk_Fa_G.empty()) {
          // Screened Coulomb line, Q + dQ
          const auto qk_screened = get_qk(int(k)) + dQ_screening(int(k), ie);
          for (auto nu = 0ul; nu < num_spinor; ++nu) {
            ComplexRMatrix Fa_G(radial_zero);
            for (auto t = 0ul; t < num_spinor; ++t) {
              Fa_G +=
                mult_rows_full(G_delta[i_delta].radial(t, nu), Fa.component(t));
            }
            qk_Fa_G.push_back(qk_screened * Fa_G);
          }
        }

        const auto ql = get_qk(int(l)).dri();
        for (auto mu = 0ul; mu < num_spinor; ++mu) {
          const auto dG_Fa_ql = dG_Fa[mu] * ql;
          for (auto nu = 0ul; nu < num_spinor; ++nu) {
            dSigma.sp(mu, nu) +=
              (i_d_epsilon * mult_elements(dG_Fa_ql, qk_Fa_G[nu]))
                .real()
                .Rmatrix();
          }
        }
      }
    }
  }

  // Second term of Gamma, G^e_ij(e_a + e) P^a_j2 (delta = a):
  // dSigma_12 = [ (dG_1i * [q^k F_a G_gamma]_1i) q^l ]_12 F_a(r_2)
  // The part in round brackets is summed over k, beta, gamma, for each
  // l and spinor component mu: [l][mu]. Formed when first needed
  std::vector<std::vector<ComplexRMatrix>> dG_qk_Fa_G(num_ks);
  for (auto k = 0ul; k < num_ks; ++k) {
    for (auto i_gamma = 0ul; i_gamma < num_kappas; ++i_gamma) {
      const auto kappa_gamma = Angular::kindex_to_kappa(i_gamma);
      if (Angular::Ck_kk(int(k), kappa_gamma, kappa_a) == 0.0)
        continue;

      // [q^k F_a G_gamma]^s_1i
      //   = sum_t q^k_1j F_a^t(r_j) G_gamma^{s t}_ij
      // (uses G^{s t}_ij = G^{t s}_ji). Formed when first needed
      std::vector<ComplexRMatrix> qk_Fa_G;

      for (auto l = 0ul; l < num_ks; ++l) {
        for (auto i_beta = 0ul; i_beta < num_kappas; ++i_beta) {
          const auto &dG = dG_beta_k(i_beta, k);
          if (!dG)
            continue;
          const auto L =
            L_ang_dSigma_ax(kappa_v, Angular::kindex_to_kappa(i_beta),
                            kappa_gamma, kappa_a, int(k), int(l));
          if (L == 0.0)
            continue;

          if (qk_Fa_G.empty()) {
            // Screened Coulomb line, Q + dQ
            const auto qk_screened = get_qk(int(k)) + dQ_screening(int(k), ie);
            for (auto s = 0ul; s < num_spinor; ++s) {
              ComplexRMatrix Fa_G(radial_zero);
              for (auto t = 0ul; t < num_spinor; ++t) {
                Fa_G += mult_rows_full(G_gamma[i_gamma].radial(t, s),
                                       Fa.component(t));
              }
              qk_Fa_G.push_back(qk_screened * Fa_G);
            }
          }
          if (dG_qk_Fa_G[l].empty()) {
            dG_qk_Fa_G[l].assign(num_spinor, radial_zero);
          }
          for (auto mu = 0ul; mu < num_spinor; ++mu) {
            for (auto s = 0ul; s < num_spinor; ++s) {
              dG_qk_Fa_G[l][mu] +=
                L * mult_elements(dG->radial(mu, s), qk_Fa_G[s]);
            }
          }
        }
      }
    }
  }
  for (auto l = 0ul; l < num_ks; ++l) {
    if (dG_qk_Fa_G[l].empty())
      continue;
    const auto ql = get_qk(int(l)).dri();
    for (auto mu = 0ul; mu < num_spinor; ++mu) {
      const auto dG_qk_Fa_G_ql = dG_qk_Fa_G[l][mu] * ql;
      for (auto nu = 0ul; nu < num_spinor; ++nu) {
        dSigma.sp(mu, nu) +=
          (i_d_epsilon * mult_cols_full(dG_qk_Fa_G_ql, Fa.component(nu)))
            .real()
            .Rmatrix();
      }
    }
  }

  return dSigma;
}

//==============================================================================

std::pair<GMatrix, GMatrix>
FeynmanSR::dSigma_ax_cx(const DiracOperator::TensorOperator &h, int kappa_w,
                        int kappa_v, double en_v, double omega,
                        const ExternalField::CorePolarisation *dV) {

  // dSigma^(ax)_12(e_v) = -int de/2pi dG_1i(e_v + e) Q_1j Gamma_ij2(e) Q_i2,
  // Gamma_ij2(e) = i sum_a [ P^a_ij G^e_j2(e_a - e) + G^e_ij(e_a + e) P^a_j2 ]
  //
  // (c_x) is the mirror image of (a_x): the same diagram with w and v
  // swapped and omega -> -omega (so at energy e_w = e_v + omega),
  // transposed, times the symmetry sign of the operator. Its perturbed
  // line is G(e_v + e) h G(e_v + e + omega): the same Green's functions as
  // for (a_x), as are those of Gamma, so the two are formed together.

  if (m_screen_Coulomb && !has_qpiq()) {
    form_qpiq(m_ident);
  }

  // Core polarisation of the external field vertex. dV is solved at
  // |omega|: for negative omega the field is its conjugate component.
  // For the mirror image the field has frequency -omega
  const auto dV_kappa = dV_matrices(h, dV, omega < 0.0);
  const auto dV_kappa_mirror = dV_matrices(h, dV, -omega < 0.0);

  const GMatrix zero(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

  constexpr std::complex<double> I{0.0, 1.0};
  const auto num_kappas = m_max_ki + 1;
  const auto &core = m_HF->core();

  // BLAS must run single-threaded inside the omp region below
  const qip::SingleThreadBlas single_thread_blas{};

  const auto num_threads = std::size_t(omp_get_max_threads());
  std::vector<GMatrix> dSigma_ax_ts(num_threads, zero);
  std::vector<GMatrix> mirror_ts(num_threads, zero);

#pragma omp parallel for
  for (auto ie = 0ul; ie < m_wgrid_points.size(); ie++) {

    auto &dSigma_ax_t = dSigma_ax_ts[std::size_t(omp_get_thread_num())];
    auto &mirror_t = mirror_ts[std::size_t(omp_get_thread_num())];

    // Integration energy, along the imaginary axis
    const auto epsilon = std::complex{m_omre, m_wgrid_points[ie]};
    const auto d_epsilon = I * m_wgrid_weights[ie];
    // With the factor of i from Gamma
    const auto i_d_epsilon = I * d_epsilon;

    // Green's functions of the perturbed line, each partial wave
    std::vector<ComplexGMatrix> G_omega, G_0;
    G_omega.reserve(num_kappas);
    G_0.reserve(num_kappas);
    for (auto ik = 0ul; ik < num_kappas; ++ik) {
      const auto kappa = Angular::kindex_to_kappa(ik);
      G_omega.push_back(green(kappa, en_v + epsilon + omega));
      G_0.push_back(green(kappa, en_v + epsilon));
    }

    // (a_x): G(e_v + e + omega) h G(e_v + e), between w and v
    const auto dG_ax =
      perturbed_greens_summed(h, kappa_w, kappa_v, G_omega, G_0, dV_kappa);
    // mirror image: G(e_v + e) h G(e_v + e + omega), between v and w
    const auto dG_mirror = perturbed_greens_summed(h, kappa_v, kappa_w, G_0,
                                                   G_omega, dV_kappa_mirror);

    for (const auto ia : m_core_index) {
      const auto &Fa = core[ia];
      const auto en_a = std::complex<double>{Fa.en()};
      // Hole-particle: the excited electron feels the hole in a
      const auto *Fa_hp = m_hole_particle ? &Fa : nullptr;

      // Excited Green's functions of Gamma, each partial wave
      std::vector<ComplexGMatrix> G_delta, G_gamma;
      G_delta.reserve(num_kappas);
      G_gamma.reserve(num_kappas);
      for (auto ik = 0ul; ik < num_kappas; ++ik) {
        const auto kappa = Angular::kindex_to_kappa(ik);
        G_delta.push_back(green_excited(kappa, en_a - epsilon, Fa_hp));
        G_gamma.push_back(green_excited(kappa, en_a + epsilon, Fa_hp));
      }

      dSigma_ax_t += dSigma_ax_core_term(dG_ax, kappa_v, Fa, G_delta, G_gamma,
                                         ie, i_d_epsilon);
      mirror_t += dSigma_ax_core_term(dG_mirror, kappa_w, Fa, G_delta, G_gamma,
                                      ie, i_d_epsilon);
    }
  }

  GMatrix dSigma_ax = zero;
  GMatrix mirror = zero;
  for (auto i = 0ul; i < num_threads; ++i) {
    dSigma_ax += dSigma_ax_ts[i];
    mirror += mirror_ts[i];
  }

  // -1/pi = -(1/2pi)*2; minus sign of the exchange diagram, and 2 from the
  // symmetric +/- Im(e) contour halves (integrand at -Im(e) is conjugate of
  // that at +Im(e), and we take Re part)
  dSigma_ax *= (-1.0 / M_PI);
  mirror *= (-1.0 / M_PI);

  // Relative sign between <w||h||v> and <v||h||w>
  const auto symm_sign = h.symm_sign(DiracSpinor(0, kappa_w, m_grid),
                                     DiracSpinor(0, kappa_v, m_grid));

  return {dSigma_ax, double(symm_sign) * mirror.transpose_drj()};
}

//==============================================================================
// E1 operator only
//==============================================================================

double L_ang_dSR_QdPiQ(int K, int k, int l, int ka, int kb, int kc) {
  // Angular Factor for
  // dPi ~ G^a G^b h G^c
  //        a
  //k -<-\-->--/--- l
  //   c  \   / b
  //       \./
  //        |
  //        K
  const auto tj_ka = Angular::twoj_k(ka);
  const auto tj_kb = Angular::twoj_k(kb);
  const auto tj_kc = Angular::twoj_k(kc);
  const auto sjs = Angular::sixj_2(2 * K, 2 * l, 2 * k, tj_ka, tj_kc, tj_kb);
  const auto Ck_Ca = Angular::Ck_kk(k, kc, ka);
  const auto Cl_aB = Angular::Ck_kk(l, ka, kb);
  const auto s = Angular::neg1pow(l + k + K + 1);
  return s * Ck_Ca * Cl_aB * sjs;
}

//==============================================================================

ComplexRMatrix FeynmanSR::perturbed_polarisation_kl(int k, int l,
                                                    std::complex<double> omega,
                                                    double Omega_external,
                                                    bool hole_particle) const {

  // Perturbed polarisation operator, dPi^(b)_kl(omega), for the E1 operator:
  // Tr int G_12(e) [G(e + omega) h G(e + omega - Omega)]_21 de/2pi
  //
  // Sum of the residues at the core poles of each Green's function:
  // dPi ~ G G h P + G P h G + P G h G
  // with G = G^e + G^c, G^c over the polarised core states only
  // (the pieces with three core states sum to zero over the three terms)
  ComplexRMatrix dpi_k(m_i0, m_stride, m_subgrid_points, m_grid);
  // Tensor Rank of External field h
  const auto K = 1;

  // Generalised arguments for the triple Green's function integral
  const auto alpha = Omega_external - omega;
  const auto beta = Omega_external;
  const auto gamma = 0.0;
  const auto D_ab = alpha - beta;
  const auto D_bc = beta - gamma;
  const auto D_ac = alpha - gamma;
  const auto D_ba = -D_ab;
  const auto D_cb = -D_bc;
  const auto D_ca = -D_ac;

  const auto Iunit = std::complex<double>{0.0, 1.0};
  const auto &core = m_HF->core();

  // core state summation
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &Fa = core[ia];
    if (Fa.n() < m_min_core_n)
      continue;
    auto const &Pa = m_pa[ia];
    const auto ea = std::complex<double>{Fa.en()};
    // not m_hole_particle, as need both for "screen only"
    const auto *Fa_hp = hole_particle ? &Fa : nullptr;
    for (auto iB = 0ul; iB <= m_max_ki; ++iB) {
      const auto kB = Angular::kindex_to_kappa(iB);
      for (auto iC = 0ul; iC <= m_max_ki; ++iC) {
        const auto kC = Angular::kindex_to_kappa(iC);

        // G^B(alpha-gamma) [G^C(beta-gamma) h P^a]

        // angular factor for dipole operator
        const auto h_Ca = -Angular::Ck_kk(K, kC, Fa.kappa());

        const auto c_ang_ggp =
          h_Ca * L_ang_dSR_QdPiQ(K, k, l, kB, kC, Fa.kappa());

        if (c_ang_ggp != 0.0) {

          auto g_B_ac = green_loop(kB, ea + D_ac, Fa_hp);
          auto g_C_bc = green_loop(kC, ea + D_bc, Fa_hp);

          if (!m_include_G) {
            g_B_ac.drop_g();
            g_C_bc.drop_g();
          }
          const auto gC_h_Pa = perturbed_greens(g_C_bc, Pa);

          dpi_k += c_ang_ggp * trace_product(g_B_ac, gC_h_Pa);
        }

        // G^B(alpha-beta) [P^a h G^C(gamma-beta)]

        // angular factor for dipole operator
        const auto h_aC = -Angular::Ck_kk(K, Fa.kappa(), kC);
        const auto c_ang_gpg =
          h_aC * L_ang_dSR_QdPiQ(K, k, l, kB, Fa.kappa(), kC);

        if (c_ang_gpg != 0) {
          auto g_B_ab = green_loop(kB, ea + D_ab, Fa_hp);
          auto g_C_cb = green_loop(kC, ea + D_cb, Fa_hp);

          if (!m_include_G) {
            g_B_ab.drop_g();
            g_C_cb.drop_g();
          }

          const auto Pa_h_gC = perturbed_greens(Pa, g_C_cb);

          dpi_k += c_ang_gpg * trace_product(g_B_ab, Pa_h_gC);
        }

        // P^a [G^B(beta-alpha) h G^C(gamma-alpha)]

        // angular factor for dipole
        auto h_BC = -Angular::Ck_kk(K, kB, kC);
        const auto c_ang_pgg =
          h_BC * L_ang_dSR_QdPiQ(K, k, l, Fa.kappa(), kB, kC);

        if (c_ang_pgg != 0) {
          auto g_B_ba = green_loop(kB, ea + D_ba, Fa_hp);
          auto g_C_ca = green_loop(kC, ea + D_ca, Fa_hp);

          if (!m_include_G) {
            g_B_ba.drop_g();
            g_C_ca.drop_g();
          }

          const auto gB_h_gC = perturbed_greens(g_B_ba, g_C_ca);

          dpi_k += c_ang_pgg * trace_product(Pa, gB_h_gC);
        }
      }
    }
  }
  dpi_k *= Iunit;

  return dpi_k;
}

GMatrix FeynmanSR::Sigma_SR_direct2(int kv, double env, int kw, double Omega) {

  // Calculates the diagrams:
  //   .--------/------.
  //    '.   _ /__   .'
  //      '.(_____).'
  // +
  //   .-----------------.
  //    '.   _______   .'
  //      '.(___.___).'
  //           /
  //          /

  if (!has_qpiq()) {
    form_qpiq(m_ident);
  }
  GMatrix Sigma(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

  constexpr std::complex<double> I{0.0, 1.0};
  const auto num_kappas = m_max_ki + 1;
  const auto K = 1; // Rank of the interaction (dipole)

  // BLAS must run single-threaded inside the omp region below
  const qip::SingleThreadBlas single_thread_blas{};

// Tell OpenMP how to reduce GMatrix
#pragma omp declare reduction(+ : GMatrix : omp_out += omp_in)                 \
  initializer(omp_priv = omp_orig)

#pragma omp parallel for collapse(2) reduction(+ : Sigma)
  for (auto iw = 0ul; iw < m_wgrid_points.size(); iw++) {
    for (auto iB = 0ul; iB < num_kappas; ++iB) {

      const auto omega = std::complex{m_omre, m_wgrid_points[iw]};

      // I, since dw is on imag. grid; 2 from symmetric +/- w
      // const auto dw = I * weight * m_wgrid_points.drdu(iw);
      // const auto weight = iw % 2 == 0 ? 4.0 / 3 : 2.0 / 3;
      const auto dw = I * m_wgrid_weights[iw];

      const auto kB = Angular::kindex_to_kappa(iB);

      // Silly, but ig gB includes G, then so will gB_QPQ
      const auto gB =
        m_include_G ? green(kB, env + omega) : green(kB, env + omega).drop_g();

      for (auto k = 0ul; int(k) <= m_max_k; k++) {
        for (auto l = 0ul; int(l) <= m_max_k; l++) {

          // Angular factors
          const auto ck_Bw = Angular::Ck_kk(int(l), kB, kw);
          const auto ck_vB = Angular::Ck_kk(int(k), kv, kB);
          const auto sjs_sf = L_ang_dSigma_bc(K, kv, kw, kB, int(k), int(l));

          const auto c_ang_dw = dw * sjs_sf * ck_vB * ck_Bw;

          if (c_ang_dw == 0.0)
            continue;

          // Calculates the top and bottom parts
          const auto dpi_bottom = perturbed_polarisation_kl(
            int(k), int(l), omega, Omega, m_hole_particle);

          // Field on the other loop line: k and l swapped, r1 and r2
          // swapped, at Omega - omega
          const std::complex<double> sign_top{
            double(Angular::neg1pow(int(k) + int(l) + K)), 0.0};
          const auto dpi_top =
            sign_top * perturbed_polarisation_kl(int(l), int(k), Omega - omega,
                                                 Omega, m_hole_particle);

          // Coulomb lines: q^k * dPi_kl * q^l (with integration measures)
          const auto qdpq = get_qk(int(k)) *
                            (dpi_bottom + dpi_top.transpose()) *
                            get_qk(int(l)).dri();
          Sigma += mult_elements(c_ang_dw * gB, qdpq).real();
        }
      }
    }
  }
  Sigma *= (1.0 / M_PI);
  return Sigma;
}

//======================================================================
//======================== EXCHANGE SR DIAGRAMS ========================
//======================================================================

// Unfinished: not declared in the class, and does not compile
// ComplexGMatrix FeynmanSR::GammaFw(int kv, double env, int kw,
//                                 const DiracSpinor *Fw) const {
//
//   ComplexGMatrix Gma_Fw(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);
//   const auto Iunit = std::complex<double>{0.0, 1.0};
//   const auto &core = m_HF->core();
//   for (auto ia = 0ul; ia < core.size(); ++ia) {
//     const auto &Fa = core[m_core_index[ia]];
//     const auto ea_minus_w = std::complex<double>{Fa.en()} - omega;
//     const auto ea_plus_w = std::complex<double>{Fa.en()} + omega;
//     const auto *Fa_hp = m_hole_particle ? &Fa : nullptr;
//     const auto Pa = m_pa[ia];
//     const auto Gex_minus = m_include_G ?
//                              green_excited(kv, env - omega, Fa_hp) :
//                              green_excited(kv, env - omega, Fa_hp).drop_g();
//     const auto Gex_plus = m_include_G ? green_excited(kv, env + omega, Fa_hp) :
//                                         green_excited(kv, env + omega).drop_g();
//     const auto Gex_plusFw = Gex_plus * Fw;
//     const auto PaFw = Pa * Fw;
//
//     Gma_Fw += mult_cols_by(Pa, Gex_plusFw) + mult_cols_by(Gex_minus, PaFw);
//   }
//
//   return Gma_Fw;
// }

// GMatrix FeynmanSR::Sigma_SR_exchange_BC(int kv, double env, int kw,
//                                    double Omega) const {
//   // Coded for dipole operator,
//   // kappa_v, en, kappa_w, Omega-> external field
//   //<v|Sigma(env)|w > = <v|gA pi gB| w >
//   // This is the diagram:
//   //        ..
//   //      ,'  ',
//   // --,-x-'--,--'--
//   //    '...'
//   //        ..
//   //      ,'  ',
//   // --,--'--,-x-'--
//   //    '...'

//   GMatrix Sigma(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

//   constexpr std::complex<double> I{0.0, 1.0};
//   const auto num_kappas = std::size_t(m_max_ki + 1);

// // Tell OpenMP how to reduce GMatrix
// #pragma omp declare reduction(+ : GMatrix : omp_out += omp_in)
//   initializer(omp_priv = omp_orig)

// #pragma omp parallel for collapse(2) reduction(+ : Sigma)
//   for (auto iw = 0ul; iw < m_wgrid.num_points(); iw++) {
//     for (auto iB = 0ul; iB < num_kappas; ++iB) {
//       const auto kB = Angular::kindex_to_kappa(int(iB));
//       for (auto iA = 0ul; iA < num_kappas; ++iA) {
//         const auto kA = Angular::kindex_to_kappa(int(iA));
//         // for (auto s : {-1.0, 1.0}) {
//         const auto omega = std::complex{m_omre, m_wgrid(iw)};

//         // Simpson's rule: Implicit ends (integrand zero at w=0 and w>wmax)
//         const auto weight = iw % 2 == 0 ? 4.0 / 3 : 2.0 / 3;
//         // I, since dw is on imag. grid; 2 from symmetric +/- w
//         const auto dw = I * weight * m_wgrid.drdu(iw);
//         const auto gA = m_include_G ? green(kA, env + omega + Omega) :
//                                       green(kA, env + omega + Omega).drop_g();
//         const auto gB = m_include_G ? green(kB, env + omega) :
//                                       green(kB, env + omega).drop_g();

//         const auto dG_AB = perturbed_greens(gA, gB);

//         for (auto k = 0ul; int(k) <= m_max_k; k++) {

//           // For doing single k (tests and for fk factors)
//           // if (in_k && *in_k != int(k))
//           // continue;

//           // - 1 for dipole constant factor
//           const auto h_AB = -1 * Angular::Ck_kk(int(1), kA, kB);
//           if (h_AB == 0.0)
//             continue;

//           const auto ck_Bw = Angular::Ck_kk(int(k), kB, kw);
//           const auto ck_vA = Angular::Ck_kk(int(k), kv, kA);

//           if (ck_Bw == 0.0 or ck_vA == 0.0)
//             continue;

//           const auto &qpq_dw = m_qpiq_wk[iw][k];

//           const auto LAng = L_ang_dSR(kv, kw, kB, kA, int(k));

//           if (LAng == 0.0)
//             continue;

//           const auto c_ang_dw = LAng * h_AB * dw * ck_Bw *
//                                 ck_vA; // / double(Angular::twoj_k(kv) + 1);

//           Sigma += (c_ang_dw * mult_elements(dG_AB, qpq_dw)).real();
//           // }
//         }
//       }
//     }
//   }
//   // Extra 2 from symmetric +/-w
//   Sigma *= (m_wgrid.du() / (M_PI));

//   return Sigma;
// }

// ComplexRMatrix FeynmanSR::perturbed_exchange_kl(int k, int l,
//                                               std::complex<double> omega,
//                                               double Omega,
//                                               bool hole_particle) const {

//   // modifiedPolarisation operator is Tr G_12(omega + Omega) dG_21(omega + en) d_omega
//   // Modified/perturbed polarisation operator
//   // For dPi^(b) ~ GGhP + GPhG + PGhG
//   // Should be the sameish for dPi^(c)
//   // dPi ~ G^e G_e h P + G^e G^c h P + G^c G^e h P
//   //     + G^e P h G^e + G^e P h G^c + G^c P h G^e
//   //     + P G^e h G^e + P G^e h G^c + P G^c h G^e

//   ComplexRMatrix dpi_k(m_i0, m_stride, m_subgrid_points, m_grid);
//   // Tensor Rank of External field h
//   const auto K = 1;

//   // Generalised arguments for the triple Green's function integral
//   const auto alpha = Omega;
//   const auto beta = omega + Omega;
//   const auto gamma = omega;
//   const auto D_ab = alpha - beta;
//   const auto D_bc = beta - gamma;
//   const auto D_ac = alpha - gamma;

//   // auto green1 = core ? green_core : green_excited;

//   const auto Iunit = std::complex<double>{0.0, 1.0};
//   const auto &core = m_HF->core();
//   // const auto ck_an = Angular::Ck_kk(k, Fa.kappa(), kA);

//   // core state summation
//   for (auto ia = 0ul; ia < core.size(); ++ia) {
//     const auto &Fa = core[ia];
//     if (Fa.n() < m_min_core_n)
//       continue;
//     const auto &Pa = m_pa[ia];
//     const auto ea = std::complex<double>{Fa.en()};
//     // not m_hole_particle, as need both for "screen only"
//     const auto *Fa_hp = hole_particle ? &Fa : nullptr;
//     for (int iB = 0; iB <= m_max_ki; ++iB) {
//       const auto kB = Angular::kindex_to_kappa(iB);
//       for (int iC = 0; iC <= m_max_ki; ++iC) {
//         const auto kC = Angular::kindex_to_kappa(iC);
//         // Gex^B(beta - alpha) Gex^C (gamma-alpha)) h P^A

//         // angular factor for dipole operator
//         const auto h_Ca = Angular::Ck_kk(k, kC, Fa.kappa());

//         const auto c_ang_ggp =
//           h_Ca * L_ang_dSR_GGP(K, k, l, Fa.kappa(), kB, kC);

//         if (c_ang_ggp != 0.0) {
//           auto gex_B_ba = green_excited(kB, ea - D_ab, Fa_hp);
//           auto gex_C_ca = green_excited(kC, ea - D_ac, Fa_hp);
//           auto gcor_B_ba = green_core(kB, ea - D_ab);
//           auto gcor_C_ca = green_core(kB, ea - D_ac);

//           if (!m_include_G) {
//             gex_B_ba.drop_g();
//             gex_C_ca.drop_g();
//             gcor_B_ba.drop_g();
//             gcor_C_ca.drop_g();
//           }
//           const auto dG_gexC_h_Pa = perturbed_greens(gex_C_ca, Pa);
//           const auto dG_gcorC_h_Pa = perturbed_greens(gcor_C_ca, Pa);

//           const auto dPi1 = mult_elements(gex_B_ba, dG_gexC_h_Pa) +
//                             mult_elements(gex_B_ba, dG_gcorC_h_Pa) +
//                             mult_elements(gcor_B_ba, dG_gexC_h_Pa);
//           dpi_k += c_ang_ggp * dPi1.tr();
//         }

//         // Gex^B(alpha-beta) P^A h Gex^C(gamma-beta)

//         // angular factor for dipole operator
//         const auto h_aC = Angular::Ck_kk(k, Fa.kappa(), kC);

//         const auto c_ang_gpg =
//           h_aC * L_ang_dSR_GPG(K, k, l, Fa.kappa(), kB, kC);

//         if (c_ang_gpg != 0) {
//           auto gex_B_ab = green_excited(kB, ea + D_ab, Fa_hp);
//           auto gex_C_cb = green_excited(kC, ea - D_bc, Fa_hp);
//           auto gcor_B_ab = green_core(kB, ea + D_ab);
//           auto gcor_C_cb = green_core(kB, ea - D_bc);

//           if (!m_include_G) {
//             gex_B_ab.drop_g();
//             gex_C_cb.drop_g();
//             gcor_B_ab.drop_g();
//             gcor_C_cb.drop_g();
//           }

//           auto dG_Pa_h_gexC = perturbed_greens(Pa, gex_C_cb);
//           auto dG_Pa_h_gcorC = perturbed_greens(Pa, gcor_C_cb);
//           auto dPi2 = mult_elements(gex_B_ab, dG_Pa_h_gexC) +
//                       mult_elements(gex_B_ab, dG_Pa_h_gcorC) +
//                       mult_elements(gcor_B_ab, dG_Pa_h_gexC);
//           dpi_k += c_ang_gpg * dPi2.tr();
//         }
//         //  P^A Gex^B(alpha-gamma) h Gex^C(beta-gamma)

//         // angular factor for dipole
//         auto h_BC = Angular::Ck_kk(k, kB, kC);

//         const auto c_ang_pgg =
//           h_BC * L_ang_dSR_PGG(K, k, l, Fa.kappa(), kB, kC);

//         if (c_ang_pgg != 0) {
//           auto gex_B_ac = green_excited(kB, ea + D_ac, Fa_hp);
//           auto gex_C_bc = green_excited(kC, ea + D_bc, Fa_hp);
//           auto gcor_B_ac = green_core(kB, ea + D_ac);
//           auto gcor_C_bc = green_core(kB, ea + D_bc);

//           if (!m_include_G) {
//             gex_B_ac.drop_g();
//             gex_C_bc.drop_g();
//             gcor_B_ac.drop_g();
//             gcor_C_bc.drop_g();
//           }

//           const auto dG_gexB_h_gexC = perturbed_greens(gex_B_ac, gex_C_bc);
//           const auto dG_gexB_h_gcorC = perturbed_greens(gex_B_ac, gcor_C_bc);
//           const auto dG_gcorB_h_gexC = perturbed_greens(gcor_B_ac, gex_C_bc);

//           const auto dPi3 = mult_elements(Pa, dG_gexB_h_gexC + dG_gexB_h_gcorC +
//                                                 dG_gcorB_h_gexC);
//           dpi_k += c_ang_pgg * dPi3.tr();
//         }
//       }
//     }
//   }
//   dpi_k *= Iunit;
//   return dpi_k;
// }

// GMatrix FeynmanSR::Sigma_SR_Exchange2(int kv, double env, int kw,
//                                     double Omega) const {
//   // Calculates the diagram:
//   //        ..
//   //      ,'  ',
//   // --,--'-x-,--'--
//   //    '...'
//   //  Requires findign dGamma,
//   // NOTE: NOT IMPLEMENTED FULLY

//   GMatrix Sigma(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

//   constexpr std::complex<double> I{0.0, 1.0};
//   const auto num_kappas = std::size_t(m_max_ki + 1);

//   // Tell OpenMP how to reduce GMatrix
//   // #pragma omp declare reduction(+ : GMatrix : omp_out += omp_in)
//   //   initializer(omp_priv = omp_orig)

//   // #pragma omp parallel for collapse(2) reduction(+ : Sigma)
//   for (auto iw = 0ul; iw < m_wgrid.num_points(); iw++) {
//     for (auto iB = 0ul; iB < num_kappas; ++iB) {

//       const auto omega = std::complex{m_omre, m_wgrid(iw)};

//       // Simpson's rule: Implicit ends (integrand zero at w=0 and w>wmax)
//       const auto weight = iw % 2 == 0 ? 4.0 / 3 : 2.0 / 3;

//       // I, since dw is on imag. grid; 2 from symmetric +/- w
//       const auto dw = I * weight * m_wgrid.drdu(iw);

//       const auto kB = Angular::kindex_to_kappa(int(iB));

//       // Silly, but ig gB includes G, then so will gB_QPQ
//       const auto gB =
//         m_include_G ? green(kB, env + omega) : green(kB, env + omega).drop_g();

//       for (auto k = 0ul; int(k) <= m_max_k; k++) {
//         for (auto l = 0ul; int(l) <= m_max_k; l++) {

//           const auto ck_Bw = Angular::Ck_kk(int(l), kB, kw);
//           const auto ck_vB = Angular::Ck_kk(int(k), kv, kB);
//           if (ck_vB * ck_Bw == 0.0)
//             continue;

//           // Calculates the top and bottom parts
//           const auto qdpq_top = perturbed_polarisation_kl(
//             int(k), int(l), omega, Omega, m_hole_particle);
//           const auto qdpq_bottom = perturbed_polarisation_kl(
//             int(k), int(l), -omega, Omega, m_hole_particle);

//           const auto qdpq = qdpq_top + qdpq_bottom.transpose();
//           const auto c_ang_dw = dw * ck_vB * ck_vB;

//           Sigma += mult_elements(c_ang_dw * gB, qdpq).real();
//         }
//       }
//     }
//   }

//   // Extra 2 from symmetric + / -w
//   Sigma *= (m_wgrid.du() / M_PI);

//   return Sigma;
// }

} // namespace MBPT
