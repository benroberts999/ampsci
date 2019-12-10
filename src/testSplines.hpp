#include "Dirac/DiracSpinor.hpp"
#include "Dirac/Wavefunction.hpp"
#include "IO/ChronoTimer.hpp"
#include "Maths/BSplines.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/AtomData.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>

inline auto form_spline_basis(const int kappa, const std::size_t n_states,
                              const std::size_t k_spl, const double r0_spl,
                              const double rmax_spl, const Grid &rgrid,
                              const double alpha) {
  //
  const auto imin = static_cast<std::size_t>(std::abs(kappa));
  const auto n_spl = n_states + imin + 1;
  const auto imax = n_spl - 1;

  // uses sepperate B-splines for each partial wave! OK?
  BSplines bspl(n_spl, k_spl, rgrid, r0_spl, rmax_spl);
  bspl.derivitate();

  std::pair<std::vector<DiracSpinor>, std::vector<DiracSpinor>> basis_both;
  auto &[basis, d_basis_g] = basis_both;

  // auto n_count = 1;
  for (auto i = imin; i < imax; i++) {
    basis.emplace_back(0, kappa, rgrid);
    auto &phi = basis.back();

    auto Bi = bspl.get_spline(i);
    auto dBi = bspl.get_spline_deriv(i);
    phi.f = Bi;
    auto gtmp = NumCalc::mult_vectors(rgrid.inverse_r(), Bi);
    NumCalc::scaleVec(gtmp, double(kappa));
    phi.g = NumCalc::add_vectors(dBi, gtmp);
    NumCalc::scaleVec(phi.g, 0.5 * alpha);

    auto [p0, pinf] = bspl.get_ends(i);
    phi.pinf = pinf;
    phi.p0 = p0;
  }

  for (auto i = imin; i < imax; i++) {
    basis.emplace_back(0, kappa, rgrid);
    auto &phi = basis.back();

    auto Bi = bspl.get_spline(i);
    auto dBi = bspl.get_spline_deriv(i);
    phi.g = Bi;
    auto ftmp = NumCalc::mult_vectors(rgrid.inverse_r(), Bi);
    NumCalc::scaleVec(ftmp, double(-kappa));
    phi.f = NumCalc::add_vectors(dBi, ftmp);
    NumCalc::scaleVec(phi.f, 0.5 * alpha);

    auto [p0, pinf] = bspl.get_ends(i);
    phi.pinf = pinf;
    phi.p0 = p0;
  }

  // n_count = 1;
  for (auto i = imin; i < imax; i++) {
    d_basis_g.emplace_back(0, kappa, rgrid);
    auto &phi = d_basis_g.back();

    auto Bi = bspl.get_spline(i);
    auto dBi = bspl.get_spline_deriv(i);
    auto d2Bi = bspl.get_spline_deriv2(i);
    auto dBior = NumCalc::mult_vectors(rgrid.inverse_r(), dBi);
    auto tmp = NumCalc::mult_vectors(rgrid.inverse_r(), Bi);
    auto Bor2 = NumCalc::mult_vectors(rgrid.inverse_r(), tmp);
    NumCalc::scaleVec(dBior, double(kappa));
    NumCalc::scaleVec(Bor2, double(-kappa));
    phi.f = NumCalc::add_vectors(d2Bi, dBior, Bor2);
    NumCalc::scaleVec(phi.f, 0.5 * alpha);

    auto [p0, pinf] = bspl.get_ends(i);
    phi.pinf = pinf;
    phi.p0 = p0;
  }

  for (auto i = imin; i < imax; i++) {
    d_basis_g.emplace_back(0, kappa, rgrid);
    auto &phi = d_basis_g.back();

    auto dBi = bspl.get_spline_deriv(i);
    phi.f = dBi;

    auto [p0, pinf] = bspl.get_ends(i);
    phi.pinf = pinf;
    phi.p0 = p0;
  }

  return basis_both;
}

inline auto form_basis(const int n_spl, const std::string &states_str,
                       const std::size_t k_spl, const double r0_spl,
                       const double rmax_spl, const Wavefunction &wf) {
  //
  std::vector<DiracSpinor> basis;
  std::vector<DiracSpinor> basis_positron;

  auto nklst = AtomData::listOfMaxn_k(states_str);

  auto Hd = DirectHamiltonian(wf.vnuc, wf.vdir, wf.get_alpha());

  for (const auto &nk : nklst) {
    auto kappa = nk.k;
    auto l = AtomData::l_k(kappa);
    auto max_n = nk.n;
    auto min_n = l + 1;
    // auto n_states = max_n - min_n + 1;
    // std::cout << max_n << " " << min_n << " " << n_states << " " << k_spl
    //           << "\n";
    auto [b_basis, d_basis] = form_spline_basis(
        kappa, n_spl, k_spl, r0_spl, rmax_spl, wf.rgrid, wf.get_alpha());

    LinAlg::SqMatrix Aij((int)b_basis.size());
    LinAlg::SqMatrix Sij((int)b_basis.size());
#pragma omp parallel for
    for (auto i = 0; i < (int)b_basis.size(); i++) {
      const auto &si = b_basis[i];
      auto VexPsi_i = HartreeFock::vex_psia_any(si, wf.core_orbitals);
      for (auto j = 0; j < (int)b_basis.size(); j++) {
        const auto &sj = b_basis[j];
        auto VexPsi_j = HartreeFock::vex_psia_any(sj, wf.core_orbitals);

        // auto aij = Hd.matrixEl_noD1(si, sj) +
        //            0.5 * ((si * VexPsi_j) + (sj * VexPsi_i));
        // aij -= (si * d_basis[j] + d_basis[i] * sj) / wf.get_alpha();

        auto aij =
            Hd.matrixEl(si, sj) + 0.5 * ((si * VexPsi_j) + (sj * VexPsi_i));

        Aij[i][j] = aij;
        Sij[i][j] = si * sj;
      }
    }

    auto [e_values, e_vectors] = LinAlg::realSymmetricEigensystem(Aij, Sij);

    const auto negmc2 = -1.0 / (wf.get_alpha() * wf.get_alpha());
    auto pqn = min_n - 1;
    auto pqn_pstrn = -min_n + 1;
    for (int i = 0; i < e_values.n; i++) {
      // std::cout << i << " " << pqn << " " << pqn_pstrn << "\n";
      const auto &en = e_values[i];
      const auto &pvec = e_vectors[i];
      if (en > negmc2)
        pqn++;
      else
        pqn_pstrn--;
      if (en > negmc2 && pqn > max_n) {
        break; // XXX positron states too!
      }
      auto &phi = (en > negmc2)
                      ? basis.emplace_back(pqn, kappa, wf.rgrid)
                      : basis_positron.emplace_back(pqn_pstrn, kappa, wf.rgrid);
      phi.en = en;
      for (std::size_t ib = 0; ib < b_basis.size(); ++ib) {
        // XXX This isn't working!
        phi += pvec[ib] * b_basis[ib];
      }
    }

  } // end loop over kappa

  return basis;
}
