#include "MBPT/FeynmanSR.hpp"
#include "Angular/Wigner369j.hpp"
#include "MBPT/Feynman.hpp"
#include "MBPT/RadialMatrix.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <cmath>
#include <complex>

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

double L_ang_dSR(int kw, int kv, int kB, int kA, int k) {
  const auto K = 1; //rank of operator
  const auto tj_w = Angular::twoj_k(kw);
  const auto tj_v = Angular::twoj_k(kv);
  const auto tj_B = Angular::twoj_k(kB);
  const auto tj_A = Angular::twoj_k(kA);

  const auto sjs = Angular::sixj_2(tj_w, 2 * K, tj_v, tj_B, 2 * k, tj_A);
  const auto s = Angular::neg1pow_2(tj_A + tj_B + k * 2 + K * 2);
  return s * sjs;
}

GMatrix FeynmanSR::Sigma_SR_direct(int kv, double env, int kw, double Omega) {
  // Coded for dipole operator,
  // kappa_v, en, kappa_w, Omega-> external field
  //<v|Sigma(env)|w > = <v|gA pi gB| w >
  // This is the diagram:
  //      /
  // -,--/--,-
  //   '.O.'

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

          const auto LAng = L_ang_dSR(kv, kw, kB, kA, int(k));

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

// ================WORK IN PROGRESS================================

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
  const auto s = Angular::neg1pow(l + k + 1);
  return s * Ck_Ca * Cl_aB * sjs;
}

//==============================================================================

ComplexGMatrix mult_elements_spinor(ComplexGMatrix a, const ComplexGMatrix &b) {
  a.ff().mult_elements_by(b.ff());
  if (a.includes_g()) {
    if (b.includes_g()) {
      a.fg().mult_elements_by(b.fg());
      a.gf().mult_elements_by(b.gf());
      a.gg().mult_elements_by(b.gg());
    } else {
      // b has zero g parts
      a.fg().zero();
      a.gf().zero();
      a.gg().zero();
    }
  }
  return a;
}

ComplexRMatrix FeynmanSR::perturbed_polarisation_kl(int k, int l,
                                                    std::complex<double> omega,
                                                    double Omega_external,
                                                    bool hole_particle) const {

  // modifiedPolarisation operator is Tr G_12(omega + Omega) dG_21(omega + en) d_omega
  // Modified/perturbed polarisation operator
  // For dPi^(b) ~ GGhP + GPhG + PGhG
  // Should be the sameish for dPi^(c)
  // dPi ~ G^e G^e h P + G^e G^c h P + G^c G^e h P
  //     + G^e P h G^e + G^e P h G^c + G^c P h G^e
  //     + P G^e h G^e + P G^e h G^c + P G^c h G^e
  ComplexRMatrix dpi_k(m_i0, m_stride, m_subgrid_points, m_grid);
  // Tensor Rank of External field h
  const auto K = 1;

  // Generalised arguments for the triple Green's function integral
  const auto alpha = omega;
  const auto beta = Omega_external;
  const auto gamma = 0.0;
  const auto D_ab = alpha - beta;
  const auto D_bc = beta - gamma;
  const auto D_ac = alpha - gamma;
  const auto D_ba = -D_ab;
  const auto D_cb = -D_bc;
  const auto D_ca = -D_ac;

  // auto green1 = core ? green_core : green_excited;

  const auto Iunit = std::complex<double>{0.0, 1.0};
  const auto &core = m_HF->core();
  // const auto ck_an = Angular::Ck_kk(k, Fa.kappa(), kA);

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

        // angular factor for dipole operator
        const auto h_Ca = -Angular::Ck_kk(K, kC, Fa.kappa());

        const auto c_ang_ggp =
          h_Ca * L_ang_dSR_QdPiQ(K, k, l, kB, kC, Fa.kappa());

        if (c_ang_ggp != 0.0) {

          auto gex_B_ac = green_excited(kB, ea + D_ac, Fa_hp);
          auto gcor_B_ac = green_core(kB, ea + D_ac);

          auto gcor_C_bc = green_core(kC, ea + D_bc);
          auto gex_C_bc = green_excited(kC, ea + D_bc, Fa_hp);

          if (!m_include_G) {
            gex_B_ac.drop_g();
            gex_C_bc.drop_g();
            gcor_B_ac.drop_g();
            gcor_C_bc.drop_g();
          }
          const auto gexC_h_Pa = perturbed_greens(gex_C_bc, Pa);
          const auto gcorC_h_Pa = perturbed_greens(gcor_C_bc, Pa);

          const auto dPi1 = mult_elements_spinor(gex_B_ac, gexC_h_Pa) +
                            mult_elements_spinor(gex_B_ac, gcorC_h_Pa) +
                            mult_elements_spinor(gcor_B_ac, gexC_h_Pa);
          dpi_k += c_ang_ggp * dPi1.tr();
        }

        // angular factor for dipole operator
        const auto h_aC = -Angular::Ck_kk(K, Fa.kappa(), kC);
        const auto c_ang_gpg =
          h_aC * L_ang_dSR_QdPiQ(K, k, l, kB, Fa.kappa(), kC);

        if (c_ang_gpg != 0) {
          auto gex_B_ab = green_excited(kB, ea + D_ab, Fa_hp);
          auto gcor_B_ab = green_core(kB, ea + D_ab);

          auto gcor_C_cb = green_core(kC, ea + D_cb);
          auto gex_C_cb = green_excited(kC, ea + D_cb, Fa_hp);

          if (!m_include_G) {
            gex_B_ab.drop_g();
            gex_C_cb.drop_g();
            gcor_B_ab.drop_g();
            gcor_C_cb.drop_g();
          }

          auto Pa_h_gexC = perturbed_greens(Pa, gex_C_cb);
          auto Pa_h_gcorC = perturbed_greens(Pa, gcor_C_cb);
          auto dPi2 = mult_elements_spinor(gex_B_ab, Pa_h_gexC) +
                      mult_elements_spinor(gex_B_ab, Pa_h_gcorC) +
                      mult_elements_spinor(gcor_B_ab, Pa_h_gexC);

          dpi_k += c_ang_gpg * dPi2.tr();
        }
        //  P^A Gex^B(alpha-gamma) h Gex^C(beta-gamma)

        // angular factor for dipole
        auto h_BC = -Angular::Ck_kk(K, kB, kC);
        const auto c_ang_pgg =
          h_BC * L_ang_dSR_QdPiQ(K, k, l, Fa.kappa(), kB, kC);

        if (c_ang_pgg != 0) {
          auto gex_B_ba = green_excited(kB, ea + D_ba, Fa_hp);
          auto gcor_B_ba = green_core(kB, ea + D_ba);

          auto gex_C_ca = green_excited(kC, ea + D_ca, Fa_hp);
          auto gcor_C_ca = green_core(kC, ea + D_ca);

          if (!m_include_G) {
            gex_B_ba.drop_g();
            gex_C_ca.drop_g();
            gcor_B_ba.drop_g();
            gcor_C_ca.drop_g();
          }

          const auto gexB_h_gexC = perturbed_greens(gex_B_ba, gex_C_ca);
          const auto gexB_h_gcorC = perturbed_greens(gex_B_ba, gcor_C_ca);
          const auto gcorB_h_gexC = perturbed_greens(gcor_B_ba, gex_C_ca);

          const auto dPi3 =
            mult_elements_spinor(Pa, gexB_h_gexC + gexB_h_gcorC + gcorB_h_gexC);
          dpi_k += c_ang_pgg * dPi3.tr();
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
          const auto sf = Angular::neg1pow(K + kv + kw);

          const auto tj_v = Angular::twoj_k(kv);
          const auto tj_w = Angular::twoj_k(kw);
          const auto tj_B = Angular::twoj_k(kB);

          const auto sjs_sf = sf * Angular::sixj_2(2 * K, tj_w, tj_v, tj_B,
                                                   2 * int(l), 2 * int(k));

          const auto c_ang_dw = dw * sjs_sf * ck_vB * ck_Bw;

          if (c_ang_dw == 0.0)
            continue;

          // Calculates the top and bottom parts
          const auto qdpq_bottom = perturbed_polarisation_kl(
            int(k), int(l), omega, Omega, m_hole_particle);

          const auto fudge = 1.0;
          const auto qdpq_top =
            fudge * perturbed_polarisation_kl(int(k), int(l), Omega - omega,
                                              Omega, m_hole_particle);

          const auto qdpq = (qdpq_bottom + qdpq_top.transpose());
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
