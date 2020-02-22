#include "ExternalField.hpp"
#include "Angular/Angular.hpp"
#include "Angular/Wigner_369j.hpp"
#include "Dirac/DiracOperator.hpp"
#include "Dirac/DiracSpinor.hpp"
#include "Dirac/Wavefunction.hpp"
#include "HF/HartreeFockClass.hpp"
#include "IO/ChronoTimer.hpp"
#include <algorithm>
#include <fstream>
#include <vector>
//
#include "Dirac/BSplineBasis.hpp" // XXX add to make
#include <memory>

//******************************************************************************
ExternalField::ExternalField(const DiracOperator *const h,
                             const std::vector<DiracSpinor> &core,
                             const std::vector<double> &vl, const double alpha)
    : m_h(h), p_core(&core), m_vl(vl), m_alpha(alpha), m_rank(h->rank()),
      m_pi(h->parity()), m_imag(h->imaginaryQ()) //
// w>0 typically. Allowed to be -ve for tests?
{

  m_X.resize(core.size());
  for (auto ic = 0u; ic < core.size(); ic++) {
    const auto &phic = core[ic];
    const auto pi_ch = phic.parity() * m_pi;
    const auto tj_c = phic.twoj();
    const auto tjmin_tmp = tj_c - 2 * m_rank;
    const auto tjmin = tjmin_tmp < 1 ? 1 : tjmin_tmp;
    const auto tjmax = tj_c + 2 * m_rank;
    for (int tj = tjmin; tj <= tjmax; tj += 2) {
      const auto l_minus = (tj - 1) / 2;
      const auto pi_chla = Wigner::parity_l(l_minus) * pi_ch;
      const auto l = (pi_chla == 1) ? l_minus : l_minus + 1;
      // XXX Add option to restrict l = lc ?
      const auto kappa = Wigner::kappa_twojl(tj, l);
      m_X[ic].emplace_back(0, kappa, *(phic.p_rgrid));
      m_X[ic].back().pinf = phic.pinf;
    }
  }
  m_Y = m_X;

  // Fill 6s symbol look-up table
  auto max_twoj = [](const auto &a, const auto &b) {
    return a.twoj() < b.twoj();
  };
  auto max_tj_core = std::max_element(core.cbegin(), core.cend(), max_twoj);
  auto max_tj_dPsi = max_tj_core->twoj() + 2 * m_rank;
  m_6j.fill(m_rank, max_tj_dPsi);
}

//******************************************************************************
void ExternalField::reZero() {
  for (auto &mx : m_X) {
    for (auto &m : mx) {
      m *= 0.0;
    }
  }
  m_Y = m_X;
}

//******************************************************************************
std::size_t ExternalField::core_index(const DiracSpinor &phic) {
  // XXX Note: no bounds checking here!
  return static_cast<std::size_t>(
      std::find(p_core->cbegin(), p_core->cend(), phic) - p_core->cbegin());
  // Better to just return itorator/pointer??
}

//******************************************************************************
std::vector<DiracSpinor> &ExternalField::get_dPsis(const DiracSpinor &phic,
                                                   dPsiType XorY) {
  const auto index = core_index(phic);
  return XorY == dPsiType::X ? m_X[index] : m_Y[index];
}

//******************************************************************************
const DiracSpinor &ExternalField::get_dPsi_x(const DiracSpinor &phic,
                                             dPsiType XorY, const int kappa_x) {
  const auto &dPsis = get_dPsis(phic, XorY);
  auto match_kappa_x = [=](const auto &phi) { return phi.k == kappa_x; };
  return *std::find_if(dPsis.cbegin(), dPsis.cend(), match_kappa_x);
}

//******************************************************************************
void ExternalField::solve_TDHFcore(const double omega, const int max_its,
                                   const bool print) {

  // nb: dV_ab_rhs takes 5-10x longer than any other part of this routine!

  // ChronoTimer timer("solve_TDHFcore");

  const double converge_targ = 1.0e-9;
  const auto damper = rampedDamp(0.5, 0.25, 1, 10);
  // const auto damper = rampedDamp(0.5, 0.5, 4, 10);

  const bool staticQ = std::abs(omega) < 1.0e-10;

  const auto imag = m_h->imaginaryQ();

  // The h Psi_c terms don't change, so calculate them just once?
  auto hPsi = m_X;
  for (auto ic = 0ul; ic < p_core->size(); ic++) {
    const auto &phic = (*p_core)[ic];
    for (auto j = 0ul; j < hPsi[ic].size(); j++) {
      const auto &Xx = m_X[ic][j];
      hPsi[ic][j] = m_h->reduced_rhs(Xx.k, phic);
    }
  }

  auto eps = 0.0;
  for (int it = 1; it <= max_its; it++) {
    eps = 0.0;
    const auto a_damp = (it == 1) ? 0.0 : damper(it);

    // eps for solveMixedState - doesn't need to be small!
    const auto eps_ms = (it == 1) ? 1.0e-6 : 1.0e-2;

    auto tmp_X = m_X;
    auto tmp_Y = m_Y;
#pragma omp parallel for
    for (auto ic = 0ul; ic < p_core->size(); ic++) {
      const auto &phic = (*p_core)[ic];
      auto eps_c = 0.0;

      // delta_en: always same, usually zero..
      const auto de0 = m_h->reducedME(phic, phic);
      const auto de1 = dV_ab(phic, phic);
      const auto de1_dag = dV_ab(phic, phic, true);

      const auto dePsic = (de0 + de1) * phic;
      const auto dePsic_dag = (de0 + de1_dag) * phic;

      for (auto j = 0ul; j < tmp_X[ic].size(); j++) {
        auto &Xx = tmp_X[ic][j];
        const auto &oldX = m_X[ic][j];
        const auto &hPsic = hPsi[ic][j];
        const auto dVpsic = dV_ab_rhs(Xx, phic, false);
        auto rhs = hPsic + dVpsic;
        if (Xx.k == phic.k && !imag)
          rhs -= dePsic;
        auto s = imag ? -1 : 1; // why is this needed??
        HartreeFock::solveMixedState(Xx, phic, omega, m_vl, m_alpha, *p_core,
                                     s * rhs, eps_ms);
        Xx = a_damp * oldX + (1.0 - a_damp) * Xx;
        const auto delta = (Xx - oldX) * (Xx - oldX) / (Xx * Xx);
        if (delta > eps_c)
          eps_c = delta;
      }
      if (!staticQ) {
        for (auto j = 0ul; j < tmp_Y[ic].size(); j++) {
          auto &Yx = tmp_Y[ic][j];
          const auto &oldY = m_Y[ic][j];
          const auto &hPsic = hPsi[ic][j];
          const auto dVpsic = dV_ab_rhs(Yx, phic, true);
          auto s = m_h->imaginaryQ() ? -1 : 1;
          auto rhs = (hPsic + s * dVpsic); // why here only second?
          if (Yx.k == phic.k && !imag)
            rhs -= dePsic_dag;
          HartreeFock::solveMixedState(Yx, phic, -omega, m_vl, m_alpha, *p_core,
                                       rhs, eps_ms);
          Yx = a_damp * oldY + (1.0 - a_damp) * Yx;
        }
      } else {
        auto s = m_h->imaginaryQ() ? -1 : 1;
        for (auto j = 0ul; j < tmp_Y[ic].size(); j++) {
          tmp_Y[ic][j] = s * tmp_X[ic][j];
        }
      }
#pragma omp critical(compare)
      if (eps_c > eps)
        eps = eps_c;
    }
    m_X = tmp_X;
    m_Y = tmp_Y;
    if (print)
      printf("TDHF (w=%.3f): %2i  %.1e\r", omega, it, eps);
    std::cout << std::flush;
    if (it > 0 && eps < converge_targ)
      break;
  }
  if (print)
    std::cout << "\n";
}

//******************************************************************************
// does it matter if a or b is in the core?
double ExternalField::dV_ab(const DiracSpinor &phi_n, const DiracSpinor &phi_m,
                            bool conj) {
  // conj = phi_m.en >= phi_n.en; // auto conj! XXX "correct" but broken?
  conj = phi_m.en <= phi_n.en; // auto conj! XXX "wrong"? but works?
  auto s = conj && m_h->imaginaryQ() ? -1 : 1;
  return s * phi_n * dV_ab_rhs(phi_n, phi_m, conj);
}

//******************************************************************************
DiracSpinor ExternalField::dV_ab_rhs(const DiracSpinor &phi_n,
                                     const DiracSpinor &phi_m, bool conj) {

  const auto ChiType = conj ? dPsiType::X : dPsiType::Y;
  const auto EtaType = conj ? dPsiType::Y : dPsiType::X;

  const auto k = m_h->rank();
  const auto tkp1 = double(2 * k + 1);

  const auto tjn = phi_n.twoj();
  const auto tjm = phi_m.twoj();
  const auto Ckala = Wigner::Ck_kk(k, phi_n.k, phi_m.k);

  auto dVFm = DiracSpinor(0, phi_n.k, *(phi_n.p_rgrid));
  dVFm.pinf = phi_m.pinf;

#pragma omp parallel for
  for (auto ib = 0ul; ib < p_core->size(); ib++) {
    const auto &phi_b = (*p_core)[ib];
    const auto tjb = phi_b.twoj();
    const auto &X_betas = get_dPsis(phi_b, ChiType);
    const auto &Y_betas = get_dPsis(phi_b, EtaType);
    auto dVFm_c = DiracSpinor(0, phi_n.k, *(phi_n.p_rgrid));
    dVFm_c.pinf = phi_m.pinf;

    for (auto ibeta = 0ul; ibeta < X_betas.size(); ++ibeta) {
      const auto &X_beta = X_betas[ibeta];
      const auto &Y_beta = Y_betas[ibeta];
      const auto tjbeta = X_beta.twoj();

      const auto Ckbeb = Wigner::Ck_kk(k, X_beta.k, phi_b.k);
      if (Ckala != 0 && Ckbeb != 0) {
        const auto Rkabcd =
            Coulomb::Rk_abcd_rhs(phi_n, phi_b, phi_m, X_beta + Y_beta, k);
        dVFm_c += (Ckala * Ckbeb / tkp1) * Rkabcd;
      }

      const auto s = Wigner::evenQ_2(tjbeta - tjm) ? 1 : -1;

      // exchange part (X):
      const auto l_min_X =
          std::max(std::abs(tjn - tjbeta), std::abs(tjm - tjb)) / 2;
      auto l_max_X = std::min((tjn + tjbeta), (tjm + tjb)) / 2;

      for (int l = l_min_X; l <= l_max_X; ++l) {
        const auto sixj = Wigner::sixj_2(tjm, tjn, 2 * k, tjbeta, tjb, 2 * l);
        if (sixj == 0)
          continue;
        const auto m1kpl = Wigner::evenQ(k + l) ? 1 : -1;
        const auto Ckba = Wigner::Ck_kk(l, phi_m.k, phi_b.k);
        const auto Ckalbe = Wigner::Ck_kk(l, phi_n.k, X_beta.k);
        if (Ckba == 0 || Ckalbe == 0)
          continue;
        const auto Rk = Coulomb::Rk_abcd_rhs(phi_n, phi_m, X_beta, phi_b, l);
        dVFm_c += (s * m1kpl * Ckba * Ckalbe * sixj) * Rk;
      }

      // exchange part (Y):
      const auto l_min_Y =
          std::max(std::abs(tjn - tjb), std::abs(tjm - tjbeta)) / 2;
      auto l_max_Y = std::min((tjn + tjb), (tjm + tjbeta)) / 2;

      for (int l = l_min_Y; l <= l_max_Y; ++l) {
        const auto sixj = Wigner::sixj_2(tjm, tjn, 2 * k, tjb, tjbeta, 2 * l);
        if (sixj == 0)
          continue;
        const auto m1kpl = Wigner::evenQ(k + l) ? 1 : -1;
        const auto Ckbea = Wigner::Ck_kk(l, phi_m.k, Y_beta.k);
        const auto Ckbal = Wigner::Ck_kk(l, phi_n.k, phi_b.k);
        if (Ckbea == 0 || Ckbal == 0)
          continue;
        const auto Rk = Coulomb::Rk_abcd_rhs(phi_n, phi_m, phi_b, Y_beta, l);
        dVFm_c += (s * m1kpl * Ckbea * Ckbal * sixj) * Rk;
      }
    }
#pragma omp critical(sum_X_core)
    { dVFm += dVFm_c; }
  }

  return dVFm;
}

//******************************************************************************

//******************************************************************************

//******************************************************************************

//******************************************************************************

//******************************************************************************
void ExternalField::solve_TDHFcore_matrix(const Wavefunction &wf,
                                          const double omega,
                                          const int max_its) {
  // static bool first = true;
  // This is just for testing?? Very slow. Should give same as reg method!

  ChronoTimer timer("solve_TDHFcore_matrix");
  // const bool staticQ = omega == 0;
  const bool staticQ = std::abs(omega) < 1.0e-10;

  const std::size_t nspl = 90;
  const std::size_t kspl = 7;
  const double rmin = 2.0e-6; //?
  const double rmax = 45.0;   //?

  // A := H - (e-w)S
  // b := -h -dV +deS_c

  auto max_ki = 0;
  for (const auto &Px : m_X) {
    for (const auto &x : Px) {
      auto ki = Angular::indexFromKappa(x.k);
      if (ki > max_ki)
        max_ki = ki;
    }
  }
  std::vector<std::vector<DiracSpinor>> basis_kappa;
  basis_kappa.reserve(max_ki + 1);
  for (int ki = 0; ki <= max_ki; ki++) {
    auto k = Angular::kappaFromIndex(ki);
    basis_kappa.push_back(SplineBasis::form_spline_basis(
        k, nspl, kspl, rmin, rmax, *((*p_core)[0].p_rgrid), m_alpha));
  }

  const auto imag = m_h->imaginaryQ();

  const double converge_targ = 1.0e-4;
  const auto damper = rampedDamp(0.5, 0.5, 4, 15);

  auto eps = 0.0;
  for (int it = 0; it < max_its; it++) {
    ChronoTimer timer2("solve_TDHFcore: iterations");
    eps = 0.0;
    const auto a_damp = (it == 0) ? 0.0 : damper(it);

    auto tmp_X = m_X;
    auto tmp_Y = m_Y;
#pragma omp parallel for
    for (std::size_t ic = 0; ic < p_core->size(); ic++) {
      const auto &phic = (*p_core)[ic];
      auto &dPsiX = tmp_X[ic];
      auto &dPsiY = tmp_Y[ic];
      const auto de0 = m_h->reducedME(phic, phic);
      const auto de1 = dV_ab(phic, phic);
      const auto de1_dag = dV_ab(phic, phic, true);
      for (auto ibeta = 0ul; ibeta < dPsiX.size(); ++ibeta) {
        auto &Xx = dPsiX[ibeta];
        auto &Yx = dPsiY[ibeta];
        const auto ki = Angular::indexFromKappa(Xx.k);
        const auto &basis = basis_kappa[ki];

        LinAlg::Vector bi_X((int)basis.size());
        LinAlg::Vector bi_Y((int)basis.size());
        for (int i = 0; i < (int)basis.size(); ++i) {
          const auto &xi = basis[i];
          // fill LHS vector, b
          const auto hi = m_h->reducedME(xi, phic);
          const auto hidag = m_h->reducedME(phic, xi); //??
          const auto dV = dV_ab(xi, phic);
          const auto dV_dag = dV_ab(xi, phic, true);
          const auto s = imag ? -1 : 1;
          const auto Sic = (xi.k == phic.k && !imag) ? (xi * phic) : 0.0;
          const auto deS = (de0 + de1) * Sic;
          const auto deS_dag = (de0 + de1_dag) * Sic;
          bi_X[i] = -s * hi - dV + deS; // why s here? check above??
          // bi_Y[i] = -s * (s * hi + dV_dag) + deS_dag;
          bi_Y[i] = -s * hidag - s * dV_dag + deS_dag; //???
        }
        const auto [Hij, Sij] = SplineBasis::fill_Hamiltonian_matrix(basis, wf);

        auto Aij_X = Hij - (phic.en + omega) * Sij;
        auto Aij_Y = Hij - (phic.en - omega) * Sij;

        auto s = m_h->imaginaryQ() ? -1 : 1;
        const auto c_X = LinAlg::solve_Axeqb(Aij_X, bi_X);
        const auto c_Y = staticQ ? s * c_X : LinAlg::solve_Axeqb(Aij_Y, bi_Y);

        Xx.scale(a_damp);
        Yx.scale(a_damp);
        for (int i = 0; i < (int)basis.size(); ++i) {
          Xx += (1.0 - a_damp) * c_X[i] * basis[i];
          Yx += (1.0 - a_damp) * c_Y[i] * basis[i];
        }
      }
    }

    for (std::size_t ic = 0; ic < p_core->size(); ic++) {
      for (auto ibeta = 0ul; ibeta < tmp_X[ic].size(); ++ibeta) {
        const auto &dphi = tmp_X[ic][ibeta];
        const auto &dphi0 = m_X[ic][ibeta];
        const auto eps_c = (dphi - dphi0) * (dphi - dphi0) / (dphi * dphi);
        if (eps_c > eps)
          eps = eps_c;
      }
    }

    m_Y = tmp_Y;
    m_X = tmp_X;
    printf("TDHF [matrix] (w=%.3f): %2i  %.1e\n", omega, it, eps);
    if (it > 0 && eps < converge_targ)
      break;
  }
}

//******************************************************************************
void ExternalField::print() const {
  std::ofstream of("dPsi.txt");
  const auto &gr = *((p_core->front()).p_rgrid);
  for (auto i = 0ul; i < gr.num_points; ++i) {
    // std::cout << i << std::endl;
    of << gr.r[i] << " ";
    for (auto ic = 0ul; ic < p_core->size(); ic++) {
      // std::cout << ic << std::endl;
      const auto &phic = (*p_core)[ic];
      of << phic.f[i] << " ";
      for (auto j = 0ul; j < m_X[ic].size(); j++) {
        // std::cout << j << std::endl;
        const auto &Xx = m_X[ic][j];
        const auto &Yx = m_Y[ic][j];
        of << Xx.f[i] << " " << Yx.f[i] << " ";
      }
    }
    of << "\n";
  }
}

//******************************************************************************

//******************************************************************************

//******************************************************************************
//******************************************************************************

//******************************************************************************

//******************************************************************************

//******************************************************************************
//******************************************************************************

//******************************************************************************
double ExternalField::dX_nm_bbe_rhs(const DiracSpinor &phi_n,
                                    const DiracSpinor &phi_m,
                                    const DiracSpinor &phi_b,
                                    const DiracSpinor &X_beta) {

  const auto k = m_h->rank();
  const auto tkp1 = double(2 * k + 1);

  double dX_nm_bbe = 0.0;

  const auto tjn = phi_n.twoj();
  const auto tjm = phi_m.twoj();
  const auto Ckala = Wigner::Ck_kk(k, phi_n.k, phi_m.k);

  const auto tjb = phi_b.twoj();
  const auto tjbeta = X_beta.twoj();
  const auto Ckbeb = Wigner::Ck_kk(k, X_beta.k, phi_b.k);

  if (Ckala != 0 && Ckbeb != 0) {
    const auto Rkabcd = Coulomb::Rk_abcd_any(phi_n, phi_b, phi_m, X_beta, k);
    dX_nm_bbe += (Ckala * Ckbeb / tkp1) * Rkabcd;
  }

  auto s = Wigner::evenQ_2(tjn + tjbeta + 2) ? 1 : -1;

  // exchange part (X):
  auto l_min_X = std::max(std::abs(tjn - tjbeta), std::abs(tjm - tjb)) / 2;
  auto l_max_X = std::min((tjn + tjbeta), (tjm + tjb)) / 2;

  l_min_X = 0;
  l_max_X = 20;

  for (int l = l_min_X; l <= l_max_X; ++l) {
    const auto sixj = Wigner::sixj_2(tjm, tjn, 2 * k, tjbeta, tjb, 2 * l);
    if (sixj == 0)
      continue;
    const auto m1kpl = Wigner::evenQ(k + l) ? 1 : -1;
    const auto Ckba = Wigner::Ck_kk(l, phi_m.k, phi_b.k);
    const auto Ckalbe = Wigner::Ck_kk(l, phi_n.k, X_beta.k);
    if (Ckba == 0 || Ckalbe == 0)
      continue;
    const auto Rk = Coulomb::Rk_abcd_any(phi_n, phi_m, X_beta, phi_b, l);
    dX_nm_bbe += (s * m1kpl * Ckba * Ckalbe * sixj) * Rk;
  }

  return dX_nm_bbe;
}

//******************************************************************************
double ExternalField::dY_nm_bbe_rhs(const DiracSpinor &phi_n,
                                    const DiracSpinor &phi_m,
                                    const DiracSpinor &phi_b,
                                    const DiracSpinor &Y_beta) {

  const auto k = m_h->rank();
  const auto tkp1 = double(2 * k + 1);

  double dY_nm_bbe = 0.0;

  const auto tjn = phi_n.twoj();
  const auto tjm = phi_m.twoj();
  const auto Ckala = Wigner::Ck_kk(k, phi_n.k, phi_m.k);

  const auto tjb = phi_b.twoj();
  const auto tjbeta = Y_beta.twoj();
  const auto Ckbeb = Wigner::Ck_kk(k, Y_beta.k, phi_b.k);

  if (Ckala != 0 && Ckbeb != 0) {
    const auto Rkabcd = Coulomb::Rk_abcd_any(phi_n, phi_b, phi_m, Y_beta, k);
    dY_nm_bbe += (Ckala * Ckbeb / tkp1) * Rkabcd;
  }

  auto s = Wigner::evenQ_2(tjn + tjbeta + 2) ? 1 : -1;

  // exchange part (Y):
  auto l_min_Y = std::max(std::abs(tjn - tjb), std::abs(tjm - tjbeta)) / 2;
  auto l_max_Y = std::min((tjn + tjb), (tjm + tjbeta)) / 2;

  l_min_Y = 0;
  l_max_Y = 20;

  for (int l = l_min_Y; l <= l_max_Y; ++l) {
    const auto sixj = Wigner::sixj_2(tjm, tjn, 2 * k, tjb, tjbeta, 2 * l);
    if (sixj == 0)
      continue;
    const auto m1kpl = Wigner::evenQ(k + l) ? 1 : -1;
    const auto Ckbea = Wigner::Ck_kk(l, phi_m.k, Y_beta.k);
    const auto Ckbal = Wigner::Ck_kk(l, phi_n.k, phi_b.k);
    if (Ckbea == 0 || Ckbal == 0)
      continue;
    const auto Rk = Coulomb::Rk_abcd_any(phi_n, phi_m, phi_b, Y_beta, l);
    dY_nm_bbe += (s * m1kpl * Ckbea * Ckbal * sixj) * Rk;
  }

  return dY_nm_bbe;
}
