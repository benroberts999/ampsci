#include "BSplineBasis.hpp"
#include "DiracOperator/include.hpp" //for Drake-Gordon
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/InputBlock.hpp"
#include "LinAlg/include.hpp"
#include "Maths/BSpline.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <utility>

namespace SplineBasis {

//==============================================================================
Parameters::Parameters(IO::InputBlock input)
  : states(input.get<std::string>("states", "")),
    n(input.get("number", 30ul)),
    k(input.get("order", 7ul)),
    r0(input.get("r0", 1.0e-4)),
    reps(input.get("r0_eps", 0.0)),
    rmax(input.get("rmax", 40.0)),
    positron(input.get("positron", std::string{""})),
    type(parseSplineType(input.get<std::string>("type", "Derevianko"))),
    orthogonalise(input.get("orthogonalise", false)),
    verbose(input.get("verbose", true)) {}

Parameters::Parameters(const std::string &istates, std::size_t in,
                       std::size_t ik, double ir0, double ireps, double irmax,
                       const std::string &ipositron, SplineType itype,
                       bool iorthogonalise, bool iverbose)
  : states(istates),
    n(in),
    k(ik),
    r0(ir0),
    reps(ireps),
    rmax(irmax),
    positron(ipositron),
    type(itype),
    orthogonalise(iorthogonalise),
    verbose(iverbose) {}

//==============================================================================
std::vector<DiracSpinor> form_basis(const Parameters &params,
                                    const Wavefunction &wf,
                                    const bool correlationsQ) {
  // Forms the pseudo-spectrum basis by diagonalising Hamiltonian over B-splines

  const auto &[states_str, n_spl, k_spl, r0_spl, r0_eps, rmax_spl,
               tmp_positron_str, basis_type, ortho, verbose] = params;

  // Not required, but helpful to be back compatible
  const std::string positron_str =
    tmp_positron_str == "false" ? "" :
    tmp_positron_str == "true"  ? states_str :
                                  tmp_positron_str;

  const auto nklst = AtomData::listOfStates_singlen(states_str);
  const auto nklst_positron = AtomData::listOfStates_singlen(positron_str);

  std::vector<DiracSpinor> basis;
  std::vector<DiracSpinor> basis_positron;
  const auto positronQ = !positron_str.empty();

  if (verbose) {
    std::cout << "\nConstructing B-spline basis with N=" << n_spl
              << ", k=" << k_spl
              << ".\n"
                 "Storing: "
              << states_str << " electron states\n";
    if (!positron_str.empty()) {
      std::cout << "and    : " << positron_str << " -ve energy states\n";
    }
    std::cout << "Using "
              << (basis_type == SplineType::Derevianko ?
                    "Derevianko (Duel Kinetic Balance)" :
                    "Johnson")
              << " type splines.\n";
  }

  for (const auto &nk : nklst) {

    const auto max_n = nk.n;
    const auto kappa = nk.k;
    const auto k = k_spl;

    // Find the corresponding entry in nklst2 with same kappa
    auto it = std::find_if(nklst_positron.begin(), nklst_positron.end(),
                           [kappa](const auto &nk2) { return nk2.k == kappa; });
    const auto max_n_positron = (it != nklst_positron.end()) ? it->n : 0;

    // Choose larger r0 depending on core density:
    const auto l = Angular::l_k(kappa);
    const auto l_tmp = std::min(l, DiracSpinor::max_l(wf.core()));

    const auto [rmin_l, rmax_l] = wf.lminmax_core_range(l_tmp, r0_eps);
    (void)rmax_l; // don't warn on unused rmax_l
    auto r0_eff = std::max(rmin_l, r0_spl);

    if (r0_eff == 0.0)
      r0_eff = wf.grid().r(0);

    const auto [spl_basis, d_basis] =
      form_spline_basis(kappa, n_spl, k, r0_eff, rmax_spl, wf.grid_sptr(),
                        wf.alpha(), basis_type);

    auto [Aij, Sij] = fill_Hamiltonian_matrix(spl_basis, d_basis, wf,
                                              correlationsQ, basis_type);
    const auto [e_values, e_vectors] = LinAlg::symmhEigensystem(Aij, Sij);

    expand_basis_orbitals(&basis, &basis_positron, spl_basis, kappa, max_n,
                          max_n_positron, e_values, e_vectors, wf);
    if (!basis.empty() && kappa < 0 && verbose) {
      printf("Spline cavity l=%i %1s: (%7.1e,%5.1f)aB.\n", l,
             AtomData::l_symbol(l).c_str(), r0_eff, rmax_spl);
    }
  }

  // Optionally add positron basis to end of basis. Prob not best way?
  if (positronQ) {
    for (const auto &Fp : basis_positron)
      basis.push_back(Fp);
  }

  return basis;
}

//==============================================================================
double check(const std::vector<DiracSpinor> &basis,
             const std::vector<DiracSpinor> &orbs, bool print_warning) {
  // Compare basis states to core/valence states
  // Check normality, orthogonality and energy differences
  // for (const auto *orbs : {&core, &valence})

  if (orbs.empty())
    return 0.0;

  std::string wrong_sign_list = "";
  double worst_dN = 0.0;
  double worst_dE = 0.0;
  std::string wFN, wFE;
  for (const auto &Fc : orbs) {
    // find corresponding basis state:
    auto pFbc = std::find(basis.cbegin(), basis.cend(), Fc);
    if (pFbc == basis.cend())
      continue;
    const auto FcFb = Fc * (*pFbc);
    if (Fc == *pFbc && FcFb < 0.0) {
      // basis state has wrong sign!
      wrong_sign_list += Fc.shortSymbol() + ",";
    }
    const auto dN = std::abs(std::abs(FcFb) - 1.0);
    const auto dE = std::abs((pFbc->en() - Fc.en()) / pFbc->en());
    if (dN > worst_dN) {
      worst_dN = dN;
      wFN = Fc.shortSymbol();
    }
    if (dE > worst_dE) {
      worst_dE = dE;
      wFE = Fc.shortSymbol();
    }
  }

  printf(" |<%3s|%3s>-1| = %.1e", wFN.c_str(), wFN.c_str(), worst_dN);
  if (worst_dN > 1.0e-3 && print_warning) {
    std::cout << "  ** OK?";
  }
  printf("\n dE/E(%3s)     = %.1e", wFE.c_str(), worst_dE);
  if (worst_dE > 1.0e-3 && print_warning) {
    std::cout << "  ** OK?";
  }
  const auto [eps, str] = DiracSpinor::check_ortho(orbs, basis);
  printf("\n %-10s    = %.1e", str.c_str(), eps);
  if (eps > 1.0e-3 && print_warning) {
    std::cout << "  ** OK?";
  }
  std::cout << "\n";
  if (wrong_sign_list != "") {
    std::cout << "Warning: Some basis states have opposite sign (e.g.): "
              << wrong_sign_list << "\n";
  }
  return std::max({worst_dN, std::abs(worst_dE), eps});
}

//==============================================================================
std::pair<std::vector<DiracSpinor>, std::vector<DiracSpinor>>
form_spline_basis(const int kappa, const std::size_t n_states,
                  const std::size_t k_spl, const double r0_spl,
                  const double rmax_spl, std::shared_ptr<const Grid> rgrid,
                  const double alpha, SplineType type) {
  // Forms the "base" basis of B-splines (DKB/Reno Method)

  auto imin = static_cast<std::size_t>(std::abs(kappa));
  auto n_spl = n_states + imin + 1; // subtract l (n_min)?
  auto imax = n_spl - 1;
  auto lambda_DKB = 1.0; // include kinetic balance term

  if (type == SplineType::Johnson) {
    imin = 0;
    n_spl = n_states + imin;
    imax = n_spl;
    lambda_DKB = 0.0;
  }

  // uses sepperate B-splines for each partial wave! OK?

  BSpline bspl(n_spl, k_spl, r0_spl, rmax_spl, BSpline::KnotDistro::loglinear);

  std::pair<std::vector<DiracSpinor>, std::vector<DiracSpinor>> out;
  auto &basis = out.first;
  auto &d_basis = out.second;
  basis.resize(2 * n_states, {0, kappa, rgrid});
  d_basis.resize(2 * n_states, {0, kappa, rgrid});

  for (auto ir = 0ul; ir < rgrid->num_points(); ++ir) {
    const auto r = rgrid->r(ir);
    auto [i0, bij] = bspl.get_nonzero(r, 2);
    for (auto i = 0ul; i < bspl.K(); ++i) {

      const auto nB = i + i0; // which basis spline

      // Only include
      if (nB < imin || nB >= imax)
        continue;

      // Set the splines
      auto &Sn1 = basis.at(nB - imin);
      // First "f-like" set
      Sn1.f(ir) = bij[i][0];
      Sn1.g(ir) =
        lambda_DKB * 0.5 * alpha * (bij[i][1] + (kappa / r) * bij[i][0]);
      // second "g-like"
      auto &Sn2 = basis.at(nB + n_states - imin);
      Sn2.f(ir) =
        lambda_DKB * 0.5 * alpha * (bij[i][1] - (kappa / r) * bij[i][0]);
      Sn2.g(ir) = bij[i][0];

      // Set the spline derivatives
      auto &dSn1 = d_basis.at(nB - imin);
      // First "f-like" set
      dSn1.f(ir) = bij[i][1];
      dSn1.g(ir) =
        lambda_DKB * 0.5 * alpha *
        (bij[i][2] + (kappa / r) * bij[i][1] - (kappa / r / r) * bij[i][0]);
      // second "g-like"
      auto &dSn2 = d_basis.at(nB + n_states - imin);
      dSn2.f(ir) =
        lambda_DKB * 0.5 * alpha *
        (bij[i][2] - (kappa / r) * bij[i][1] + (kappa / r / r) * bij[i][0]);
      dSn2.g(ir) = bij[i][1];
    }
  }

  return out;
}

//==============================================================================
std::pair<LinAlg::Matrix<double>, LinAlg::Matrix<double>>
fill_Hamiltonian_matrix(const std::vector<DiracSpinor> &spl_basis,
                        const std::vector<DiracSpinor> &d_basis,
                        const Wavefunction &wf, const bool correlationsQ,
                        SplineType type) {
  const auto n_spl = spl_basis.size();

  auto A_and_S = std::make_pair(LinAlg::Matrix<double>{n_spl, n_spl},
                                LinAlg::Matrix<double>{n_spl, n_spl});
  // auto &[Aij, Sij] = A_and_S; //XXX error w/ clang? Why?
  auto &Aij = A_and_S.first;
  auto &Sij = A_and_S.second;

  const auto excl_exch = wf.vHF() ? wf.vHF()->excludeExchangeQ() : true;
  const auto sigmaQ = wf.Sigma() != nullptr && correlationsQ;
  const auto VBr = wf.vHF() ? wf.vHF()->vBreit() : nullptr;

#pragma omp parallel for
  for (auto i = 0ul; i < Aij.rows(); i++) {
    const auto &si = spl_basis[i];
    const auto &dsi = d_basis[i];
    const auto VexSi = excl_exch ? 0.0 * si : HF::vexFa(si, wf.core());
    const auto SigmaSi = sigmaQ ? (*wf.Sigma())(si) : 0.0 * si;
    const auto BreitSi = VBr ? VBr->VbrFa(si, wf.core()) : 0.0 * si;

    for (auto j = 0ul; j <= i; j++) {
      const auto &sj = spl_basis[j];
      const auto &dsj = d_basis[j];

      auto aij = wf.H0ab(sj, dsj, si, dsi);
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
  for (auto i = 0ul; i < Aij.rows(); i++) {
    for (auto j = i + 1; j < Aij.cols(); j++) {
      Aij[i][j] = Aij[j][i];
      Sij[i][j] = Sij[j][i];
    }
  }
  // Note: This is work-around, since Breit seems not to be 100% symmetric!
  if (type == SplineType::Johnson)
    add_NotreDameBoundary(&Aij, spl_basis.front().kappa(), wf.alpha());

  return A_and_S;
}

//==============================================================================
void add_NotreDameBoundary(LinAlg::Matrix<double> *pAij, const int kappa,
                           const double alpha) {
  // XX These should be re-checked..
  // W. R. Johnson, S. A. Blundell, J. Sapirstein, Phys. Rev. A 37, 307 (1988)
  auto &Aij = *pAij;
  const auto n2 = Aij.rows();
  const auto n1 = n2 / 2;
  const auto c = 1.0 / alpha;
  // r=0 boundary conds
  // Add (c/2)*f(0)^2 + (c/2)*f(0)*g(0)
  Aij[0][0] += (kappa < 0) ? c * c : 2.0 * c * c * c;
  Aij[0][n1] += 0.5 * c;
  Aij[n1][0] += 0.5 * c;
  // f(rmax)=g(rmax)
  Aij[n1 - 1][n1 - 1] += 0.5 * c * c; // f(R)^2
  Aij[n2 - 1][n2 - 1] -= 0.5 * c * c; // g(R)^2
}

//==============================================================================
void expand_basis_orbitals(std::vector<DiracSpinor> *basis,
                           std::vector<DiracSpinor> *basis_positron,
                           const std::vector<DiracSpinor> &spl_basis,
                           const int kappa, const int max_n, int max_n_positron,
                           const LinAlg::Vector<double> &e_values,
                           const LinAlg::Matrix<double> &e_vectors,
                           const Wavefunction &wf)
// Expands the pseudo-spectrum basis in terms of B-spline basis and expansion
// coeficient found from diagonalising the Hamiltonian over Bsplns
{
  auto l = Angular::l_k(kappa);
  auto min_n = l + 1;

  const auto neg_mc2 = -1.0 / (wf.alpha() * wf.alpha());
  auto pqn = min_n - 1;
  auto pqn_pstrn = -min_n + 1;
  for (auto i = 0ul; i < e_values.rows(); i++) {
    const auto &en = e_values[i];
    const auto &pvec = e_vectors[i];
    const auto positive_energy = en > neg_mc2;
    positive_energy ? ++pqn : --pqn_pstrn;

    if ((positive_energy && pqn > max_n) ||
        (!positive_energy && pqn_pstrn < -max_n_positron))
      continue;

    auto &Fi = (positive_energy) ?
                 basis->emplace_back(pqn, kappa, wf.grid_sptr()) :
                 basis_positron->emplace_back(pqn_pstrn, kappa, wf.grid_sptr());
    Fi.en() = en;
    for (std::size_t ib = 0; ib < spl_basis.size(); ++ib) {
      Fi += pvec[ib] * spl_basis[ib];
    }
    // Check positive at low r:
    {
      const auto low_r = 0.3 / wf.Znuc();
      const auto low_i = wf.grid().getIndex(low_r);
      double fsum = Fi.f(low_i);
      if (fsum < 0.0) {
        Fi *= -1.0;
      }
    }
    const auto Fnorm = Fi * Fi;
    if (std::abs(Fnorm - 1) > 1.0e-4) {
      std::cout << "Warning: Possible spurious state: " << Fi.shortSymbol()
                << " e=" << Fi.en() << ", Norm=" << Fnorm << " (" << Fnorm - 1.0
                << ")\n";
    }
    // Note: they are not even roughly normalised...I think they should be??
    Fi.normalise();

    // find first non-zero point
    {
      auto p0 = 0ul;
      for (auto ir = 0ul; ir < Fi.grid().num_points(); ++ir) {
        if (Fi.f(ir) != 0 || Fi.g(ir) != 0) {
          p0 = ir;
          break;
        }
      }
      auto pinf = Fi.grid().num_points();
      for (auto ir = Fi.grid().num_points(); ir != 0; --ir) {
        if (Fi.f(ir - 1) != 0 || Fi.g(ir - 1) != 0) {
          pinf = ir;
          break;
        }
      }
      Fi.min_pt() = p0;
      Fi.max_pt() = pinf;
    }
  }
}

//==============================================================================
std::vector<double> sumrule_TKR(const std::vector<DiracSpinor> &basis,
                                const std::vector<double> &r, bool print) {
  std::vector<double> result;
  if (basis.empty())
    return result;

  const auto &Fa = basis.front();
  const auto max_l =
    std::max_element(cbegin(basis), cend(basis), DiracSpinor::comp_l)->l();
  result.reserve(std::size_t(max_l + 1));

  for (int l = 0; l <= max_l; l++) {
    double sum_el = 0.0;
    double sum_p = 0.0;
    for (const auto &Fn : basis) {
      if (Fn == Fa)
        continue;
      const auto f = (Fn.kappa() == l) ? l : (Fn.kappa() == -l - 1) ? l + 1 : 0;
      if (f == 0)
        continue;
      const auto Ran = Fa * (r * Fn);
      const auto term = f * (Fn.en() - Fa.en()) * Ran * Ran / (2 * l + 1);
      if (Fn.n() > 0)
        sum_el += term;
      else
        sum_p += term;
    }
    result.push_back(sum_el + sum_p);
    if (print)
      printf("l=%1i, sum = %10.6f%+10.6f = %8.1e\n", l, sum_el, sum_p,
             sum_el + sum_p);
  }
  return result;
}

//==============================================================================
std::vector<double> sumrule_DG(int nDG, const std::vector<DiracSpinor> &basis,
                               const Grid &gr, double alpha, bool print) {
  // Drake-Goldman sum rules: w^n |<a|r|b>|^2  (n=0,1,2)
  // Only up to lmax-1, since need to have states with l'=l+1

  assert(nDG < 3); // better way?

  auto rhat = DiracOperator::E1(gr);          // vector E1
  auto r2hat = DiracOperator::RadialF(gr, 2); // scalar r^2

  std::vector<double> result;
  if (basis.empty())
    return result;

  const auto max_ki =
    std::max_element(cbegin(basis), cend(basis), DiracSpinor::comp_ki)
      ->k_index();
  const auto max_l =
    std::max_element(cbegin(basis), cend(basis), DiracSpinor::comp_l)->l();

  for (int ki = 0; ki <= max_ki; ki++) {
    const auto kappa = Angular::kappaFromIndex(ki);
    auto find_ka = [=](const auto &Fn) { return Fn.kappa() == kappa; };
    const auto &Fa = *std::find_if(basis.begin(), basis.end(), find_ka);
    // need to have l_n = la+1 terms, or sum doesn't work:
    if (Fa.l() == max_l)
      continue;
    // if (print)
    //   std::cout << "kappa: " << kappa << " (" << Fa.symbol() << ")\n";
    auto sum = 0.0;
    for (const auto &Fn : basis) {
      const auto w = Fn.en() - Fa.en();
      const auto Ran = rhat.reducedME(Fa, Fn);
      const double c = 1.0 / (2 * std::abs(Fa.kappa()));
      const auto term = std::pow(w, nDG) * Ran * Ran * c;
      sum += term;
    }
    if (nDG == 2) {
      sum *= alpha * alpha / 3;
    }
    const auto s0 = (nDG == 0) ? r2hat.radialIntegral(Fa, Fa) :
                    (nDG == 1) ? 0.0 :
                                 1.0;
    if (print)
      printf("%i| k=%2i: sum=%10.5f, exact=%+10.5f, diff = %8.1e\n", nDG, kappa,
             sum, s0, sum - s0);
    result.push_back(sum - s0);
  }
  return result;
}

std::pair<double, double> r_completeness(const DiracSpinor &Fa,
                                         const std::vector<DiracSpinor> &basis,
                                         const Grid &gr, bool print) {
  // <a|a> = sum_n <a|r|n><n|1/r|a> = 1
  // <a|r^2|a> = sum_n <a|r|n><n|r|a>

  const auto r = gr.r();
  const auto r2 = gr.rpow(2);
  const auto rinv = gr.rpow(-1);

  double sum1 = 0.0;
  double sumr2 = 0.0;
  const auto sumr2_expect = Fa * (r2 * Fa);
  for (const auto &Fn : basis) {
    if (Fn.kappa() != Fa.kappa())
      continue;
    const auto arn = Fa * (r * Fn);
    const auto nra = arn; // symmetric
    const auto nrinva = Fn * (rinv * Fa);
    sum1 += arn * nrinva;
    sumr2 += arn * nra;
    // std::cout << Fn.shortSymbol() << " " << sum1 - 1.0 << " .. "
    //           << (sumr2 - sumr2_expect) / sumr2_expect << "\n";
  }
  const auto eps1 = std::abs((sum1 - 1.0) / 1.0);
  const auto epsr2 = std::abs((sumr2 - sumr2_expect) / sumr2_expect);
  if (print) {
    printf("%3s : %.0f [%.0f] %.1e | %9.4f [%9.4f] %.1e\n",
           Fa.shortSymbol().c_str(), sum1, 1.0, eps1, sumr2, sumr2_expect,
           epsr2);
  }
  return {eps1, epsr2};
}

} // namespace SplineBasis
