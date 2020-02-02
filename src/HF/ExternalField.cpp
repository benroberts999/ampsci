#include "ExternalField.hpp"
#include "Angular/Angular.hpp"
#include "Angular/Wigner_369j.hpp"
#include "Dirac/DiracOperator.hpp"
#include "Dirac/DiracSpinor.hpp"
#include "Dirac/Wavefunction.hpp"
#include "HF/HartreeFockClass.hpp"
#include "IO/ChronoTimer.hpp"
#include <algorithm>
#include <vector>
//
#include "Dirac/BSplineBasis.hpp" // XXX add to make

//******************************************************************************
ExternalField::ExternalField(const DiracOperator *const h,
                             const std::vector<DiracSpinor> &core,
                             const std::vector<double> &vl, const double alpha,
                             const double omega)
    : m_h(h), p_core(&core), m_vl(vl), m_alpha(alpha), m_omega(omega), //
      static_fieldQ(std::abs(omega) > 0 ? false : false),              //
      m_rank(h->rank()), m_pi(h->parity()), m_imag(h->imaginaryQ())    //
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
    }
  }
  // if (!static_fieldQ)
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
void ExternalField::solve_TDHFcore(int max_its) {

  ChronoTimer timer("solve_TDHFcore");

  const double converge_targ = 1.0e-3;
  const auto damper = rampedDamp(0.75, 0.5, 4, 15);

  const auto imag = m_h->imaginaryQ();

  // the h Psi_c terms don't change, so calculate them just once?
  // No quicker, and much messier...
  // auto hPsi = m_X;
  // for (auto ic = 0ul; ic < p_core->size(); ic++) {
  //   const auto &phic = (*p_core)[ic];
  //   for (auto j = 0ul; j < hPsi[ic].size(); j++) {
  //     const auto &Xx = m_X[ic][j];
  //     hPsi[ic][j] = m_h->reduced_rhs(Xx.k, phic);
  //   }
  // }

  auto eps = 0.0;
  for (int it = 0; it < max_its; it++) {
    eps = 0.0;
    const auto a_damp = (it == 0) ? 0.0 : damper(it);

    auto tmp_X = m_X;
    auto tmp_Y = m_Y;
#pragma omp parallel for
    for (auto ic = 0ul; ic < p_core->size(); ic++) {
      const auto &phic = (*p_core)[ic];

      auto eps_c = 0.0;

      // these also a) always the same
      // b) always 0 in many cases!
      const auto de0 = m_h->reducedME(phic, phic);
      const auto de1 = dV_ab(phic, phic);
      const auto de1_dag = dV_ab(phic, phic, true);

      const auto dePsic = (de0 + de1) * phic;
      const auto dePsic_dag = (de0 + de1_dag) * phic;

      for (auto j = 0ul; j < tmp_X[ic].size(); j++) {
        auto &Xx = tmp_X[ic][j];
        const auto &oldX = m_X[ic][j];
        const auto hPsic = m_h->reduced_rhs(Xx.k, phic); // always same!! XX
        // const auto &hPsic = hPsi[ic][j];
        const auto dVpsic = dV_ab_rhs(Xx, phic);
        auto rhs = hPsic + dVpsic;
        if (Xx.k == phic.k && !imag)
          rhs -= dePsic;
        HartreeFock::solveMixedState(Xx, phic, m_omega, m_vl, m_alpha, *p_core,
                                     rhs);
        Xx = a_damp * oldX + (1.0 - a_damp) * Xx;
        const auto delta = std::abs((oldX * Xx) / (Xx * Xx) - 1.0);
        if (delta > eps_c)
          eps_c = delta;
      }
      for (auto j = 0ul; j < tmp_Y[ic].size(); j++) {
        auto &Yx = tmp_Y[ic][j];
        const auto &oldY = m_Y[ic][j];
        const auto hPsic = m_h->reduced_rhs(Yx.k, phic); // always same!! XX
        // const auto &hPsic = hPsi[ic][j];
        const auto dVpsic = dV_ab_rhs(Yx, phic, true);
        const auto s = m_h->imaginaryQ() ? -1 : 1;
        auto rhs = s * (hPsic + dVpsic);
        if (Yx.k == phic.k && !imag)
          rhs -= dePsic_dag;
        HartreeFock::solveMixedState(Yx, phic, -m_omega, m_vl, m_alpha, *p_core,
                                     rhs);
        Yx = a_damp * oldY + (1.0 - a_damp) * Yx;
        const auto delta = std::abs((oldY * Yx) / (Yx * Yx) - 1.0);
        if (delta > eps_c)
          eps_c = delta;
      }
#pragma omp critical(compare)
      if (eps_c > eps) {
        eps = eps_c; // worst epsilon in core
      }
    }
    m_X = tmp_X;
    m_Y = tmp_Y;
    printf("TDHF (w=%.3f): %2i  %.1e\n", m_omega, it, eps);
    if (it > 0 && eps < converge_targ)
      break;
  }
}

//******************************************************************************
// does it matter if a or b is in the core?
double ExternalField::dV_ab(const DiracSpinor &phi_n, const DiracSpinor &phi_m,
                            bool conj) {
  return phi_n * dV_ab_rhs(phi_n, phi_m, conj);
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

  auto rme_sum_dir = 0.0 * phi_n;
  auto rme_sum_exc = 0.0 * phi_n;

#pragma omp parallel for
  for (auto ib = 0ul; ib < p_core->size(); ib++) {
    const auto &phi_b = (*p_core)[ib];
    const auto tjb = phi_b.twoj();
    const auto &X_betas = get_dPsis(phi_b, ChiType);
    const auto &Y_betas = get_dPsis(phi_b, EtaType);
    auto rme_sum_dir_c = 0.0 * phi_n;
    auto rme_sum_exc_c = 0.0 * phi_n;

    for (auto ibeta = 0ul; ibeta < X_betas.size(); ++ibeta) {
      const auto &X_beta = X_betas[ibeta];
      const auto &Y_beta = Y_betas[ibeta];
      const auto tjbeta = X_beta.twoj();

      const auto Ckbeb = Wigner::Ck_kk(k, X_beta.k, phi_b.k);
      if (Ckala != 0 && Ckbeb != 0) {
        const auto Rkabcd =
            Coulomb::Rk_abcd_rhs(phi_n, phi_b, phi_m, X_beta + Y_beta, k);
        rme_sum_dir_c += (Ckala * Ckbeb) * Rkabcd;
      }

      auto s = Wigner::evenQ_2(tjn + tjb) ? 1 : -1;

      // exchange part (X):
      const auto l_min_X =
          std::max(std::abs(tjn - tjbeta), std::abs(tjm - tjb)) / 2;
      const auto l_max_X = std::min((tjn + tjbeta), (tjm + tjb)) / 2;

      for (int l = l_min_X; l <= l_max_X; ++l) {
        auto ll = 1.0 / (2 * l + 1); // XXX ???
        const auto sixj = m_6j.get_6j_mutable(tjm, tjn, tjbeta, tjb, k, l);
        if (sixj == 0)
          continue;
        const auto m1kpl = Wigner::evenQ(k + l) ? 1 : -1;
        const auto Ckba = Wigner::Ck_kk(l, phi_m.k, phi_b.k);
        const auto Ckalbe = Wigner::Ck_kk(l, phi_n.k, X_beta.k);
        if (Ckba == 0 || Ckalbe == 0)
          continue;
        const auto Rk = Coulomb::Rk_abcd_rhs(phi_n, phi_m, X_beta, phi_b, l);
        rme_sum_exc_c += (ll * s * m1kpl * Ckba * Ckalbe * sixj) * Rk;
      }

      // exchange part (Y):
      const auto l_min_Y =
          std::max(std::abs(tjn - tjb), std::abs(tjm - tjbeta)) / 2;
      const auto l_max_Y = std::min((tjn + tjb), (tjm + tjbeta)) / 2;

      for (int l = l_min_Y; l <= l_max_Y; ++l) {
        auto ll = 1.0 / (2 * l + 1); // XXX ???
        const auto sixj = m_6j.get_6j(tjm, tjn, tjb, tjbeta, k, l);
        if (sixj == 0)
          continue;
        const auto m1kpl = Wigner::evenQ(k + l) ? 1 : -1;
        const auto Ckbea = Wigner::Ck_kk(l, phi_m.k, Y_beta.k);
        const auto Ckbal = Wigner::Ck_kk(l, phi_n.k, phi_b.k);
        if (Ckbea == 0 || Ckbal == 0)
          continue;
        const auto Rk = Coulomb::Rk_abcd_rhs(phi_n, phi_m, phi_b, Y_beta, l);
        rme_sum_exc_c += (ll * s * m1kpl * Ckbea * Ckbal * sixj) * Rk;
      }
    }
#pragma omp critical(sum_X_core)
    {
      rme_sum_dir += rme_sum_dir_c;
      rme_sum_exc += rme_sum_exc_c;
    }
  }

  rme_sum_dir *= 1.0 / tkp1;

  return rme_sum_dir - rme_sum_exc;
  // return rme_sum_dir - 0.8314 * rme_sum_exc;
  // XXX Fudge factor!?
}

//******************************************************************************

//******************************************************************************

//******************************************************************************

//******************************************************************************

//******************************************************************************
void ExternalField::solve_TDHFcore_matrix(const Wavefunction &wf) {
  static bool first = true;
  // This is just for testing?? Very slow. Should give same as reg method!

  ChronoTimer timer("solve_TDHFcore_matrix");

  const std::size_t nspl = 30;
  const std::size_t kspl = 4;
  const double rmin = 1.0e-4; //?
  const double rmax = 30.0;   //?

  // A := H - (e-w)S
  // b := -h -dV +deS_c

  auto a_damp = first ? 0.0 : 0.5;

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
        const auto dV = first ? 0.0 : dV_ab(xi, phic);
        const auto dV_dag = first ? 0.0 : dV_ab(xi, phic, true);
        const auto s = imag ? -1 : 1;
        const auto Sic = (xi.k == phic.k && !imag) ? (xi * phic) : 0.0;
        const auto deS = (de0 + de1) * Sic;
        const auto deS_dag = (de0 + de1_dag) * Sic;
        bi_X[i] = -hi - dV + deS;
        bi_Y[i] = -s * (hi + dV_dag) + deS_dag;
      }
      const auto [Hij, Sij] = SplineBasis::fill_Hamiltonian_matrix(basis, wf);
      auto Aij_X = Hij - (phic.en - m_omega) * Sij;
      auto Aij_Y = Hij - (phic.en + m_omega) * Sij;

      const auto c_X = LinAlg::solve_Axeqb(Aij_X, bi_X);
      const auto c_Y = LinAlg::solve_Axeqb(Aij_Y, bi_Y);

      Xx.scale(a_damp); // or 0.5, for DAMP!
      Yx.scale(a_damp); // or 0.5, for DAMP!
      for (int i = 0; i < (int)basis.size(); ++i) {
        Xx += (1.0 - a_damp) * c_X[i] * basis[i];
        Yx += (1.0 - a_damp) * c_Y[i] * basis[i];
      }
    }
  }

  m_Y = tmp_Y;
  m_X = tmp_X;
  first = false;
}
