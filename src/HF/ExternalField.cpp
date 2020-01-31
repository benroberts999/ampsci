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
void ExternalField::solve_TDHFcore() {
  // XXX THIS should have omega inside! ?
  //
  ChronoTimer timer("solve_TDHFcore");

  auto tmp_X = m_X;
  auto tmp_Y = m_Y;

  static bool first = true;
  auto a_damp = first ? 0.0 : 0.5;

  // auto damper = rampedDamp(0.5, 0.3, 4, 10);

  // XXX Doesn't include de yet

#pragma omp parallel for
  for (auto ic = 0ul; ic < p_core->size(); ic++) {
    const auto &phic = (*p_core)[ic];
    auto &dPsis_X = tmp_X[ic];
    for (auto &Xx : dPsis_X) {
      const auto hPsic = m_h->reduced_rhs(Xx.k, phic); // always same!?
      const auto dVpsic = dV_ab_rhs(Xx, phic);
      const auto newX = HartreeFock::solveMixedState(
          Xx.k, phic, m_omega, m_vl, m_alpha, *p_core, hPsic + dVpsic);
      Xx = a_damp * Xx + (1.0 - a_damp) * newX;
    }
    auto &dPsis_Y = tmp_Y[ic];
    for (auto &Yx : dPsis_Y) {
      const auto hPsic = m_h->reduced_rhs(Yx.k, phic); // always same!?
      const auto dVpsic = dV_ab_rhs(Yx, phic, true);
      const auto s = m_h->imaginaryQ() ? -1 : 1;
      const auto newY = HartreeFock::solveMixedState(
          Yx.k, phic, -m_omega, m_vl, m_alpha, *p_core, s * (hPsic + dVpsic));
      Yx = a_damp * Yx + (1.0 - a_damp) * newY;
    }
  }

  m_X = tmp_X;
  m_Y = tmp_Y;
  first = false;
}

//******************************************************************************
// does it matter if a or b is in the core?
double ExternalField::dV_ab(const DiracSpinor &phi_alpha,
                            const DiracSpinor &phi_a, bool conj) {

  return phi_alpha * dV_ab_rhs(phi_alpha, phi_a, conj);
}

//******************************************************************************
DiracSpinor ExternalField::dV_ab_rhs(const DiracSpinor &phi_alpha,
                                     const DiracSpinor &phi_a, bool conj) {

  const auto ChiType = conj ? dPsiType::X : dPsiType::Y;
  const auto EtaType = conj ? dPsiType::Y : dPsiType::X;

  const auto k = m_h->rank();
  const auto tkp1 = double(2 * k + 1);

  const auto tjalpha = phi_alpha.twoj();
  const auto tja = phi_a.twoj();
  const auto Ckala = Wigner::Ck_kk(k, phi_alpha.k, phi_a.k);

  auto rme_sum_dir = 0.0 * phi_alpha;
  auto rme_sum_exc = 0.0 * phi_alpha;

#pragma omp parallel for
  for (auto ib = 0ul; ib < p_core->size(); ib++) {
    const auto &phi_b = (*p_core)[ib];
    const auto tjb = phi_b.twoj();
    const auto &X_betas = get_dPsis(phi_b, ChiType);
    const auto &Y_betas = get_dPsis(phi_b, EtaType);
    auto rme_sum_dir_c = 0.0 * phi_alpha;
    auto rme_sum_exc_c = 0.0 * phi_alpha;

    for (auto ibeta = 0ul; ibeta < X_betas.size(); ++ibeta) {
      const auto &X_beta = X_betas[ibeta];
      const auto &Y_beta = Y_betas[ibeta];
      const auto tjbeta = X_beta.twoj();

      const auto Ckbeb = Wigner::Ck_kk(k, X_beta.k, phi_b.k);
      if (Ckala != 0 && Ckbeb != 0) {
        const auto Rkabcd =
            Coulomb::Rk_abcd_rhs(phi_alpha, phi_b, phi_a, X_beta + Y_beta, k);
        rme_sum_dir_c += (Ckala * Ckbeb) * Rkabcd;
      }

      // exchange part (X):
      const auto l_min_X =
          std::max(std::abs(tjalpha - tjbeta), std::abs(tja - tjb)) / 2;
      const auto l_max_X = std::min((tjalpha + tjbeta), (tja + tjb)) / 2;

      for (int l = l_min_X; l <= l_max_X; ++l) {
        const auto sixj = m_6j.get_6j_mutable(tja, tjalpha, tjbeta, tjb, k, l);
        if (sixj == 0)
          continue;
        const auto m1kpl = Wigner::evenQ(k + l) ? 1 : -1;
        const auto Ckba = Wigner::Ck_kk(l, phi_a.k, phi_b.k);
        const auto Ckalbe = Wigner::Ck_kk(l, phi_alpha.k, X_beta.k);
        if (Ckba == 0 || Ckalbe == 0)
          continue;
        const auto Rk =
            Coulomb::Rk_abcd_rhs(phi_alpha, phi_a, X_beta, phi_b, l);
        rme_sum_exc_c += (m1kpl * Ckba * Ckalbe * sixj) * Rk;
      }

      // exchange part (Y):
      const auto l_min_Y =
          std::max(std::abs(tjalpha - tjb), std::abs(tja - tjbeta)) / 2;
      const auto l_max_Y = std::min((tjalpha + tjb), (tja + tjbeta)) / 2;

      auto s = Wigner::evenQ_2(tjalpha + tja) ? -1 : 1; //??? XXX GUESS??

      for (int l = l_min_Y; l <= l_max_Y; ++l) {
        const auto sixj = m_6j.get_6j(tja, tjalpha, tjb, tjbeta, k, l);
        if (sixj == 0)
          continue;
        const auto m1kpl = Wigner::evenQ(k + l) ? 1 : -1;
        const auto Ckbea = Wigner::Ck_kk(l, phi_a.k, Y_beta.k);
        const auto Ckbal = Wigner::Ck_kk(l, phi_alpha.k, phi_b.k);
        if (Ckbea == 0 || Ckbal == 0)
          continue;
        const auto Rk =
            Coulomb::Rk_abcd_rhs(phi_alpha, phi_a, phi_b, Y_beta, l);
        rme_sum_exc_c += (s * m1kpl * Ckbea * Ckbal * sixj) * Rk;
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

  const std::size_t nspl = 40;
  const std::size_t kspl = 4;
  const double rmin = 1.0e-5; //?
  const double rmax = 40.0;   //?

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

  auto tmp_X = m_X;
  auto tmp_Y = m_Y;
#pragma omp parallel for
  for (std::size_t ic = 0; ic < p_core->size(); ic++) {
    const auto &phic = (*p_core)[ic];
    auto &dPsiX = tmp_X[ic];
    auto &dPsiY = tmp_Y[ic];
    auto de = 0.0; // m_h->reducedME(phic, phic); // no dV ?
    for (auto ibeta = 0ul; ibeta < dPsiX.size(); ++ibeta) {
      auto &Xx = dPsiX[ibeta];
      auto &Yx = dPsiY[ibeta];
      auto ki = Angular::indexFromKappa(Xx.k);
      const auto &basis = basis_kappa[ki];

      LinAlg::Vector bi_X((int)basis.size());
      LinAlg::Vector bi_Y((int)basis.size());
      for (int i = 0; i < (int)basis.size(); ++i) {
        const auto &xi = basis[i];
        // fill LHS vector, b
        const auto hi = m_h->reducedME(xi, phic);
        const auto hi_dag = m_h->reducedME(xi, phic); //??? XXX
        const auto dV = first ? 0.0 : dV_ab(xi, phic);
        const auto dV_dag = first ? 0.0 : dV_ab(xi, phic, true);
        const auto s = m_h->imaginaryQ() ? -1 : 1;
        // auto deS = (xi.k == phic.k) ? de * (xi * phic) : 0.0; // reduced? OK?
        bi_X[i] = -hi - dV;               // + deS;
        bi_Y[i] = -s * (hi_dag + dV_dag); // + deS;
      }
      const auto [Hij, Sij] = SplineBasis::fill_Hamiltonian_matrix(basis, wf);
      auto Aij_X = Hij - (phic.en - m_omega) * Sij;
      auto Aij_Y = Hij - (phic.en + m_omega) * Sij;

      auto c_X = LinAlg::solve_Axeqb(Aij_X, bi_X);
      auto c_Y = LinAlg::solve_Axeqb(Aij_Y, bi_Y);

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
