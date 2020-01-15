#include "BSplineBasis.hpp"
#include "Dirac/DiracSpinor.hpp"
#include "Dirac/Hamiltonian.hpp"
// #include "Dirac/Operators.hpp"
#include "Dirac/Wavefunction.hpp"
#include "Maths/BSplines.hpp"
#include "Maths/Grid.hpp"
#include "Maths/LinAlg_MatrixVector.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/AtomData.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <utility>

//******************************************************************************
std::vector<DiracSpinor> form_basis(const std::string &states_str,
                                    const std::size_t n_spl,
                                    const std::size_t k_spl,
                                    const double r0_spl, const double rmax_spl,
                                    const Wavefunction &wf) {
  std::vector<DiracSpinor> basis;
  std::vector<DiracSpinor> basis_positron;

  const auto nklst = AtomData::listOfMaxn_k(states_str);

  std::cout << "\nConstructing B-spline basis with N=" << n_spl
            << ", k=" << k_spl << ". Storing: " << states_str << "\n";

  for (const auto &nk : nklst) {
    const auto max_n = nk.n;
    const auto kappa = nk.k;
    const auto spl_basis = form_spline_basis(
        kappa, n_spl, k_spl, r0_spl, rmax_spl, wf.rgrid, wf.get_alpha());

    auto [Aij, Sij] = fill_Hamiltonian_matrix(spl_basis, wf);

    const auto [e_values, e_vectors] =
        LinAlg::realSymmetricEigensystem(Aij, Sij);

    expand_basis_orbitals(&basis, &basis_positron, spl_basis, kappa, max_n,
                          e_values, e_vectors, wf);
  }

  if (!basis.empty()) {
    printf("Spline cavity: (%7.1e,%5.1f)aB.\n", basis.front().r0(),
           basis.front().rinf());
  }

  return basis;
  // basis_positron is ignored for now, but it is calculated
}

//******************************************************************************
std::vector<DiracSpinor>
form_spline_basis(const int kappa, const std::size_t n_states,
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

  std::vector<DiracSpinor> basis;

  for (auto i = imin; i < imax; i++) {
    auto &phi = basis.emplace_back(0, kappa, rgrid);

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
    auto &phi = basis.emplace_back(0, kappa, rgrid);

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

  return basis;
}

//******************************************************************************
std::pair<LinAlg::SqMatrix, LinAlg::SqMatrix>
fill_Hamiltonian_matrix(const std::vector<DiracSpinor> &spl_basis,
                        const Wavefunction &wf) {
  auto n_spl = (int)spl_basis.size();

  std::pair<LinAlg::SqMatrix, LinAlg::SqMatrix> A_and_S =
      std::make_pair(n_spl, n_spl);
  auto &[Aij, Sij] = A_and_S;

  auto excl_exch = wf.exclude_exchangeQ();

  // XXX Move this into wf ??
  // auto Hd = DirectHamiltonian(wf.vnuc, wf.vdir, wf.get_alpha());
  auto Hd = RadialHamiltonian(wf.rgrid, wf.get_alpha());
  Hd.set_v(-1, wf.vnuc, wf.vdir); // same each kappa
  Hd.set_v_mag(wf.Hse_mag);       // same each kappa

#pragma omp parallel for
  for (auto i = 0; i < (int)spl_basis.size(); i++) {
    const auto &si = spl_basis[i];
    auto VexPsi_i = HartreeFock::vex_psia_any(si, wf.core_orbitals);
    for (auto j = 0; j < (int)spl_basis.size(); j++) {
      const auto &sj = spl_basis[j];

      auto aij = Hd.matrixEl(si, sj);
      if (!excl_exch)
        aij += (sj * VexPsi_i);

      Aij[i][j] = aij;
      Sij[i][j] = si * sj;
    }
  }
  Aij.make_symmetric(); // very slight asymmetry from exchange pot.
  return A_and_S;
}

//******************************************************************************
void expand_basis_orbitals(std::vector<DiracSpinor> *basis,
                           std::vector<DiracSpinor> *basis_positron,
                           const std::vector<DiracSpinor> &spl_basis,
                           const int kappa, const int max_n,
                           const LinAlg::Vector &e_values,
                           const LinAlg::SqMatrix &e_vectors,
                           const Wavefunction &wf) {
  auto l = AtomData::l_k(kappa);
  auto min_n = l + 1;

  const auto neg_mc2 = -1.0 / (wf.get_alpha() * wf.get_alpha());
  auto pqn = min_n - 1;
  auto pqn_pstrn = -min_n + 1;
  for (int i = 0; i < e_values.n; i++) {
    const auto &en = e_values[i];
    const auto &pvec = e_vectors[i];
    if (en > neg_mc2) {
      if (++pqn > max_n)
        continue;
    } else {
      if (--pqn_pstrn < -max_n)
        continue;
    }
    auto &phi = (en > neg_mc2)
                    ? basis->emplace_back(pqn, kappa, wf.rgrid)
                    : basis_positron->emplace_back(pqn_pstrn, kappa, wf.rgrid);
    phi.en = en;
    phi.p0 = spl_basis[0].p0;
    phi.pinf = spl_basis[0].pinf;
    for (std::size_t ib = 0; ib < spl_basis.size(); ++ib) {
      phi += pvec[ib] * spl_basis[ib];
    }
    phi.normalise();
  }
}

//******************************************************************************

//******************************************************************************

//******************************************************************************

//******************************************************************************
