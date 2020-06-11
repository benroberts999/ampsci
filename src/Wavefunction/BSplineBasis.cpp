#include "BSplineBasis.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/UserInput.hpp"
#include "IO/safeProfiler.hpp"
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

namespace SplineBasis {

// plot
// for [i=2:6] "Cs_core.txt" u 1:i every :::0::0 w l lc i lw 2 t
// columnheader(i), for [i=2:6] "Cs_basis.txt" u 1:i every :::0::0 w l dt 0 lc i
// notitle

static constexpr bool ND_type = false;

//******************************************************************************
Parameters::Parameters(IO::UserInputBlock input)
    : states(input.get<std::string>("states", "")),
      n(input.get("number", 0ul)),
      k(input.get("order", 0ul)),
      r0(input.get("r0", 0.0)),
      reps(input.get("r0_eps", 0.0)),
      rmax(input.get("rmax", 0.0)),
      positronQ(input.get("positron", false)) {}

Parameters::Parameters(std::string istates, std::size_t in, std::size_t ik,
                       double ir0, double ireps, double irmax, bool ipositronQ)
    : states(istates),
      n(in),
      k(ik),
      r0(ir0),
      reps(ireps),
      rmax(irmax),
      positronQ(ipositronQ) {}

//******************************************************************************
std::vector<DiracSpinor> form_basis(const Parameters &params,
                                    const Wavefunction &wf,
                                    const bool correlationsQ)
// Forms the pseudo-spectrum basis by diagonalising Hamiltonian over B-splines
{
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  const auto &[states_str, n_spl, k_spl, r0_spl, r0_eps, rmax_spl, positronQ] =
      params;
  std::vector<DiracSpinor> basis;
  std::vector<DiracSpinor> basis_positron;

  const auto nklst = AtomData::listOfStates_singlen(states_str);

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

    [[maybe_unused]] // don't warn on unused rmax_l
    const auto [rmin_l, rmax_l] = wf.lminmax_core_range(l_tmp, r0_eps);
    auto r0_eff = std::max(rmin_l, r0_spl);
    if (l_tmp < l) {
      // For l's that arn't in core, make r0 20% larger (per delta l) ?? XXX
      r0_eff *= 1.0 + 0.20 * (l - l_tmp);
      const auto r0_min = l <= 1 ? 1.0e-4 : l <= 3 ? 1.0e-3 : 1.0e-2; // ?
      r0_eff = std::max(r0_eff, r0_min);
    }

    if (ND_type)
      r0_eff = l <= 1 ? 1.0e-4 : l <= 3 ? 1.0e-3 : 1.0e-2; // Notre-Dame XXX

    const auto spl_basis = form_spline_basis(kappa, n_spl, k, r0_eff, rmax_spl,
                                             wf.rgrid, wf.alpha);

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
    auto nmc2 = -1.0 / (wf.alpha * wf.alpha);
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
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  //
  auto imin = static_cast<std::size_t>(std::abs(kappa));
  auto n_spl = n_states + imin + 1; // subtract l (n_min)?
  auto imax = n_spl - 1;

  if (ND_type) {
    imin = 1; // XXX
    n_spl = n_states + imin;
    imax = n_spl;
  }

  // uses sepperate B-splines for each partial wave! OK?
  BSplines bspl(n_spl, k_spl, rgrid, r0_spl, rmax_spl);
  bspl.derivitate();

  std::vector<DiracSpinor> basis;
  basis.reserve(2 * imax);

  for (auto i = imin; i < imax; i++) {
    auto &phi = basis.emplace_back(0, kappa, rgrid);

    const auto &Bi = bspl.get_spline(i);
    const auto &dBi = bspl.get_spline_deriv(i);
    phi.f = Bi;
    if (!ND_type) {
      auto gtmp = NumCalc::mult_vectors(rgrid.rpow(-1), Bi);
      NumCalc::scaleVec(gtmp, double(kappa));
      phi.g = NumCalc::add_vectors(dBi, gtmp);
      NumCalc::scaleVec(phi.g, 0.5 * alpha);
    }

    // phi.df = dBi;
    // auto bor2 = NumCalc::mult_vectors(rgrid.rpow(-2), Bi);
    // auto dbor = NumCalc::mult_vectors(rgrid.rpow(-1), dBi);
    // NumCalc::scaleVec(bor2, -double(kappa));
    // NumCalc::scaleVec(dbor, double(kappa));
    // const auto &d2Bi = bspl.get_spline_deriv2(i);
    // phi.dg = NumCalc::add_vectors(bor2, dbor, d2Bi);
    // NumCalc::scaleVec(phi.dg, 0.5 * alpha);

    auto [p0, pinf] = bspl.get_ends(i);
    phi.pinf = pinf;
    phi.p0 = p0;
    phi.normalise(1.0);
  }

  for (auto i = imin; i < imax; i++) {
    auto &phi = basis.emplace_back(0, kappa, rgrid);

    const auto &Bi = bspl.get_spline(i);
    const auto &dBi = bspl.get_spline_deriv(i);
    phi.g = Bi;
    if (!ND_type) {
      auto ftmp = NumCalc::mult_vectors(rgrid.rpow(-1), Bi);
      NumCalc::scaleVec(ftmp, double(-kappa));
      phi.f = NumCalc::add_vectors(dBi, ftmp);
      NumCalc::scaleVec(phi.f, 0.5 * alpha);
    }

    // phi.dg = dBi;
    // auto bor2 = NumCalc::mult_vectors(rgrid.rpow(-2), Bi);
    // auto dbor = NumCalc::mult_vectors(rgrid.rpow(-1), dBi);
    // NumCalc::scaleVec(bor2, double(kappa));
    // NumCalc::scaleVec(dbor, -double(kappa));
    // const auto &d2Bi = bspl.get_spline_deriv2(i);
    // phi.df = NumCalc::add_vectors(bor2, dbor, d2Bi);
    // NumCalc::scaleVec(phi.df, 0.5 * alpha);

    auto [p0, pinf] = bspl.get_ends(i);
    phi.pinf = pinf;
    phi.p0 = p0;
    phi.normalise(1.0);
  }

  return basis;
}

//******************************************************************************
std::pair<LinAlg::SqMatrix, LinAlg::SqMatrix>
fill_Hamiltonian_matrix(const std::vector<DiracSpinor> &spl_basis,
                        const Wavefunction &wf, const bool correlationsQ) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  auto n_spl = (int)spl_basis.size();

  std::pair<LinAlg::SqMatrix, LinAlg::SqMatrix> A_and_S =
      std::make_pair(n_spl, n_spl);
  // auto &[Aij, Sij] = A_and_S; //XXX error w/ clang? Why?
  auto &Aij = A_and_S.first;
  auto &Sij = A_and_S.second;

  const auto excl_exch = wf.exclude_exchangeQ();
  const auto sigmaQ = wf.getSigma() != nullptr && correlationsQ;
  const auto VBr = wf.getHF() ? wf.getHF()->get_Breit() : nullptr;

  // Move this into wf ??
  auto Hd = RadialHamiltonian(wf.rgrid, wf.alpha);
  Hd.set_v(-1, wf.get_Vlocal(0)); // same each kappa //XXX
  Hd.set_v_mag(wf.get_Hmag(0));   // Magnetic QED form-factor [usually empty]

#pragma omp parallel for
  for (auto i = 0ul; i < Aij.n; i++) {
    const auto &si = spl_basis[i];
    const auto VexSi = excl_exch ? 0.0 * si : HF::vexFa(si, wf.core);
    const auto SigmaSi = sigmaQ ? (*wf.getSigma())(si) : 0.0 * si;
    const auto BreitSi = VBr ? (*VBr)(si) : 0.0 * si;

    for (auto j = 0ul; j <= i; j++) {
      // for (auto j = 0ul; j < Aij.n; j++) {
      const auto &sj = spl_basis[j];

      auto aij = Hd.matrixEl(sj, si);
      if (!excl_exch)
        aij += (sj * VexSi);
      if (sigmaQ)
        aij += (sj * SigmaSi);
      if (VBr)
        aij += sj * BreitSi;

      Aij[i][j] = aij;
      Sij[i][j] = sj * si;
    }
  }
  // Fill second - half of symmetric matrix
  for (auto i = 0ul; i < Aij.n; i++) {
    for (auto j = i + 1; j < Aij.n; j++) {
      Aij[i][j] = Aij[j][i];
      Sij[i][j] = Sij[j][i];
    }
  }
  // Aij.enforce_symmetric();
  // Note: This is work-around, since Breit seems not to be 100% symmetric!
  if (ND_type)
    add_NotreDameBoundary(&Aij, spl_basis.front().k, wf.alpha);

  return A_and_S;
}

//******************************************************************************
void add_NotreDameBoundary(LinAlg::SqMatrix *pAij, const int kappa,
                           const double alpha) {
  auto &Aij = *pAij;
  const auto n2 = Aij.n;
  const auto n1 = n2 / 2;
  const auto c = 1.0 / alpha;
  // r=0 boundary conds
  Aij[0][0] += (kappa < 0) ? c * c : 2.0 * c * c * c;
  Aij[0][n1] += 0.5 * c;
  Aij[n1][0] -= 0.5 * c;
  // f(rmax)=g(rmax)
  Aij[n1 - 1][n1 - 1] += 0.5 * c;
  Aij[n1 - 1][n2 - 1] -= 0.5;
  Aij[n2 - 1][n1 - 1] += 0.5;
  // Aij[n2 - 1][n2 - 1] -= 0.5 * alpha; //?
  Aij[n2 - 1][n2 - 1] -= 0.5 * c;
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
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  auto l = AtomData::l_k(kappa);
  auto min_n = l + 1;

  const auto neg_mc2 = -1.0 / (wf.alpha * wf.alpha);
  auto pqn = min_n - 1;
  auto pqn_pstrn = -min_n + 1;
  for (auto i = 0ul; i < e_values.n; i++) {
    const auto &en = e_values[i];
    const auto &pvec = e_vectors[i];
    const auto positive_energy = en > neg_mc2;
    positive_energy ? ++pqn : --pqn_pstrn;

    if ((positive_energy && pqn > max_n) ||
        (!positive_energy && pqn_pstrn < -max_n))
      continue;

    auto &phi = (positive_energy)
                    ? basis->emplace_back(pqn, kappa, wf.rgrid)
                    : basis_positron->emplace_back(pqn_pstrn, kappa, wf.rgrid);
    phi.en = en;
    phi.p0 = spl_basis[0].pinf; // yes, backwards (updated below)
    phi.pinf = spl_basis[0].p0;
    auto sign = pvec[0] > 0 ? 1 : -1; // mostly, but not completely, works
    for (std::size_t ib = 0; ib < spl_basis.size(); ++ib) {
      phi += sign * pvec[ib] * spl_basis[ib];
    }
    // Note: they are not even roughly normalised...I think they should be??
    phi.normalise();
  }

  // ensure correct sign (doesn't seem to matter.., but nicer)
  const auto ir = wf.rgrid.getIndex(0.005);
  for (auto &Fb : *basis) {
    if (Fb.f[ir] < 0)
      Fb *= -1;
  }
}

} // namespace SplineBasis
