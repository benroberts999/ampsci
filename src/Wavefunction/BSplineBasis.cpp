#include "BSplineBasis.hpp"
#include "Maths/BSplines.hpp"
#include "Maths/Grid.hpp"
#include "Maths/LinAlg_MatrixVector.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Hamiltonian.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <utility>

// size_t vs int for vectors + LinAlg..currently a mess
#pragma GCC diagnostic ignored "-Wsign-conversion"

namespace SplineBasis {

//******************************************************************************
std::vector<DiracSpinor>
form_basis(const std::string &states_str, const std::size_t n_spl,
           const std::size_t k_spl, const double r0_spl, const double r0_eps,
           const double rmax_spl, const Wavefunction &wf, const bool positronQ,
           const bool correlationsQ)
// Forms the pseudo-spectrum basis by diagonalising Hamiltonian over B-splines
{
  std::vector<DiracSpinor> basis;
  std::vector<DiracSpinor> basis_positron;

  const auto nklst = AtomData::listOfMaxn_k(states_str);

  std::cout << "\nConstructing B-spline basis with N=" << n_spl
            << ", k=" << k_spl << ". Storing: " << states_str << "\n";

  for (const auto &nk : nklst) {
    const auto max_n = nk.n;
    const auto kappa = nk.k;
    const auto kmin = std::size_t(AtomData::l_k(kappa) + 3);
    const auto k = k_spl; // < kmin ? kmin : k_spl;
    if (k_spl < kmin) {
      std::cout << "Warning: Spline order k=" << k
                << " may be small for kappa=" << kappa << " (kmin=" << kmin
                << ")\n";
    }

    // Chose larger r0 depending on core density:
    const auto l = AtomData::l_k(kappa);
    const auto l_tmp = std::min(l, wf.maxCore_l());

    const auto [rmin_l, rmax_l] = wf.lminmax_core_range(l_tmp, r0_eps);
    (void)rmax_l;
    const auto r0_eff = std::max(rmin_l, r0_spl);

    const auto spl_basis = form_spline_basis(kappa, n_spl, k, r0_eff, rmax_spl,
                                             wf.rgrid, wf.get_alpha());

    auto [Aij, Sij] = fill_Hamiltonian_matrix(spl_basis, wf, correlationsQ);
    const auto [e_values, e_vectors] =
        LinAlg::realSymmetricEigensystem(Aij, Sij);

    expand_basis_orbitals(&basis, &basis_positron, spl_basis, kappa, max_n,
                          e_values, e_vectors, wf);
    if (!basis.empty() && kappa < 0) {
      printf("Spline cavity l=%i %1s: (%7.1e,%5.1f)aB.\n", l,
             AtomData::l_symbol(l).c_str(), basis.back().r0(),
             basis.back().rinf());
    }
  }

  // Optionally add positron basis to end of basis. Prob not best way?
  if (positronQ) {
    for (const auto &Fp : basis_positron)
      basis.push_back(Fp);
  }

  for (auto &Fp : basis) {
    auto l = Fp.l();
    if (l == 0)
      continue;
    auto comp = [=](const auto &Fa) { return Fa.l() == l - 1; };
    auto prev = std::find_if(basis.begin(), basis.end(), comp);
    if (prev == basis.end())
      continue;
    auto nmc2 = -1.0 / (wf.get_alpha() * wf.get_alpha());
    if (Fp.en < prev->en && Fp.en > nmc2) {
      // XXX Better: count nodes? ['Spurious node at large r?']
      std::cout << "WARNING: "
                << "Spurious state?? " << Fp.symbol() << " " << Fp.en << "\n";
      // Fp *= 0.0;
    }
  }

  return basis;
}

//******************************************************************************
std::vector<DiracSpinor>
form_spline_basis(const int kappa, const std::size_t n_states,
                  const std::size_t k_spl, const double r0_spl,
                  const double rmax_spl, const Grid &rgrid, const double alpha)
// Forms the "base" basis of B-splines (DKB/Reno Method)
{
  //
  const auto imin = static_cast<std::size_t>(std::abs(kappa));
  const auto n_spl = n_states + imin;
  const auto imax = n_spl - 1;

  // uses sepperate B-splines for each partial wave! OK?
  BSplines bspl(n_spl, k_spl, rgrid, r0_spl, rmax_spl);
  bspl.derivitate();

  std::vector<DiracSpinor> basis;

  for (auto i = imin; i < imax; i++) {
    auto &phi = basis.emplace_back(0, kappa, rgrid);

    const auto &Bi = bspl.get_spline(i);
    const auto &dBi = bspl.get_spline_deriv(i);
    phi.f = Bi;
    auto gtmp = NumCalc::mult_vectors(rgrid.rpow(-1), Bi);
    NumCalc::scaleVec(gtmp, double(kappa));
    phi.g = NumCalc::add_vectors(dBi, gtmp);
    NumCalc::scaleVec(phi.g, 0.5 * alpha);

    auto [p0, pinf] = bspl.get_ends(i);
    phi.pinf = pinf;
    phi.p0 = p0;
  }

  for (auto i = imin; i < imax; i++) {
    auto &phi = basis.emplace_back(0, kappa, rgrid);

    const auto &Bi = bspl.get_spline(i);
    const auto &dBi = bspl.get_spline_deriv(i);
    phi.g = Bi;
    auto ftmp = NumCalc::mult_vectors(rgrid.rpow(-1), Bi);
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
                        const Wavefunction &wf, const bool correlationsQ) {
  auto n_spl = (int)spl_basis.size();

  std::pair<LinAlg::SqMatrix, LinAlg::SqMatrix> A_and_S =
      std::make_pair(n_spl, n_spl);
  // auto &[Aij, Sij] = A_and_S; //XXX error w/ clang? Why?
  auto &Aij = A_and_S.first;
  auto &Sij = A_and_S.second;

  const auto excl_exch = wf.exclude_exchangeQ();
  const auto sigmaQ = wf.m_Sigma != nullptr && correlationsQ;

  // Move this into wf ??
  auto Hd = RadialHamiltonian(wf.rgrid, wf.get_alpha());
  Hd.set_v(-1, wf.get_Vlocal(0)); // same each kappa //XXX
  Hd.set_v_mag(wf.get_Hmag(0));   // Magnetic QED form-factor [usually empty]

#pragma omp parallel for
  for (auto i = 0; i < (int)spl_basis.size(); i++) {
    const auto &si = spl_basis[i];
    const auto VexSi =
        excl_exch ? 0.0 * si : HF::vex_psia_any(si, wf.core_orbitals);
    const auto SigmaSi = sigmaQ ? (*wf.m_Sigma)(si) : 0.0 * si;

    for (auto j = 0; j < (int)spl_basis.size(); j++) {
      const auto &sj = spl_basis[j];

      auto aij = Hd.matrixEl(si, sj); // + si * VexSi + si * SigmaSi;
      if (!excl_exch)
        aij += (sj * VexSi);
      if (sigmaQ)
        aij += (sj * SigmaSi);

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
                           const Wavefunction &wf)
// Expands the pseudo-spectrum basis in terms of B-spline basis and expansion
// coeficient found from diagonalising the Hamiltonian over Bsplns
{
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
    // In some cases, first point causes issues? Unlcear
    phi.f[phi.p0] = 0.0;
    phi.g[phi.p0] = 0.0;
    phi.normalise();
  }
}

} // namespace SplineBasis
