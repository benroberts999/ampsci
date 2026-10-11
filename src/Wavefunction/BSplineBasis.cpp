#include "BSplineBasis.hpp"
#include "DiracOperator/include.hpp" //for Drake-Goldman
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
#include <limits>
#include <string>
#include <utility>

namespace SplineBasis {

//==============================================================================
Parameters::Parameters(IO::InputBlock input)
  : states(input.get<std::string>("states", "")),
    n(input.get("number", 30ul)),
    k(input.get("order", 7ul)),
    r0(input.get("r0", std::vector<double>{})),
    rmax(input.get("rmax", 40.0)),
    positron(input.get("positron", std::string{""})),
    type(parseSplineType(input.get<std::string>("type", "Derevianko"))),
    orthogonalise(input.get("orthogonalise", false)),
    verbose(input.get("verbose", true)) {
  // The same options are used for both the Basis and Spectrum blocks
  input.check(
    {{"", "Options for 'Basis' and 'Spectrum'. Spectrum is the same as "
          "Basis, but includes correlations. Spectrum is used for "
          "sum-over-states (while Basis is used for MBPT)."},
     {"number", "Number of splines used in expansion [30]"},
     {"order", "order of splines ~7-9 [7]"},
     {"r0", "minimum cavity radius (first internal knot), for each l: "
            "comma-separated list, last value used for all higher l "
            "[1.0e-4; Johnson and Fischer: 3.0e-4,1.0e-4,1.0e-3,1.0e-3,"
            "1.0e-2]"},
     {"rmax", "maximum cavity radius [40.0]"},
     {"states", "states to keep (e.g., 30spdf20ghi)"},
     {"orthogonalise",
      "Force orthogonal to core (Basis) or valence (Spectrum) [false]"},
     {"print", "Print all spline energies (for testing) [false]"},
     {"positron",
      "Basis string for negative energy states (same format as states). []"},
     {"type", "Derevianko (DKB), Shabaev (DKB), Johnson, Fischer, or Hybrid "
              "(DKB) [Derevianko]"},
     {"verbose", "Print construction details [true]"}});
}

Parameters::Parameters(const std::string &istates, std::size_t in,
                       std::size_t ik, std::vector<double> ir0, double irmax,
                       const std::string &ipositron, SplineType itype,
                       bool iorthogonalise, bool iverbose)
  : states(istates),
    n(in),
    k(ik),
    r0(std::move(ir0)),
    rmax(irmax),
    positron(ipositron),
    type(itype),
    orthogonalise(iorthogonalise),
    verbose(iverbose) {}

//==============================================================================
std::vector<double> default_r0(SplineType type) {
  switch (type) {
  case SplineType::Johnson:
  case SplineType::Fischer:
    return {3.0e-4, 1.0e-4, 1.0e-3, 1.0e-3, 1.0e-2};
  case SplineType::Derevianko:
  case SplineType::Shabaev:
  case SplineType::Hybrid:
    return {1.0e-4};
  }
  assert(false && "Unreachable: missing spline type in default_r0");
  return {1.0e-4};
}

//==============================================================================
std::vector<DiracSpinor> form_basis(const Parameters &params,
                                    const Wavefunction &wf,
                                    const bool correlationsQ) {
  // Forms the pseudo-spectrum basis by diagonalising Hamiltonian over B-splines

  const auto &[states_str, n_spl, k_spl, r0_in, rmax_spl, tmp_positron_str,
               basis_type, ortho, verbose] = params;

  // First internal knot for each l; last entry used for all higher l
  const auto r0_list = r0_in.empty() ? default_r0(basis_type) : r0_in;

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
    std::cout << "Using " << parseSplineType(basis_type) << " type splines.\n";
  }

  // Negative-energy states are only formed for kappas that are also in states
  for (const auto &nk_positron : nklst_positron) {
    const auto kappa_in_states =
      std::any_of(nklst.cbegin(), nklst.cend(), [&nk_positron](const auto &nk) {
        return nk.k == nk_positron.k;
      });
    if (!kappa_in_states) {
      std::cout << "Warning: No -ve energy states formed for "
                << DiracSpinor::shortSymbol(0, nk_positron.k)
                << ", since not included in electron states\n";
    }
  }

  for (const auto &nk : nklst) {

    const auto max_n = nk.n;
    const auto kappa = nk.k;
    const auto k = k_spl;

    // Find the corresponding entry in nklst2 with same kappa
    auto it = std::find_if(nklst_positron.begin(), nklst_positron.end(),
                           [kappa](const auto &nk2) { return nk2.k == kappa; });
    const auto max_n_positron = (it != nklst_positron.end()) ? it->n : 0;

    // First knot for this l; no closer than the first grid point (r0 > 0)
    const auto l = Angular::l_k(kappa);
    const auto il = std::min(std::size_t(l), r0_list.size() - 1);
    const auto r0_eff = std::max(r0_list.at(il), wf.grid().r0());

    const auto splines =
      form_spline_basis(kappa, n_spl, k, r0_eff, rmax_spl, wf.grid_sptr(),
                        wf.alpha(), basis_type);
    const auto &spl_basis = splines.basis;

    auto [Aij, Sij] =
      fill_Hamiltonian_matrix(splines, wf, correlationsQ, basis_type);
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

  // nk_index keys the Coulomb/MBPT tables, and must fit in DiracSpinor::Index.
  // Large n arises for negative-energy states, labelled after the electrons
  const auto max_nk_index = std::numeric_limits<DiracSpinor::Index>::max();
  for (const auto &Fn : basis) {
    if (Angular::nk_to_index(Fn.n(), Fn.kappa()) > max_nk_index) {
      std::cout << "Warning: n too large for nk_index (" << Fn.shortSymbol()
                << "): reduce number of splines or -ve energy states\n";
      break;
    }
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
    if (FcFb < 0.0) {
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
SplineParams spline_params(SplineType type, int kappa, std::size_t N_states) {
  const auto lk = static_cast<std::size_t>(std::abs(kappa));
  const auto N = N_states;

  // {n_spl, imin_f, imax_f, imin_g, imax_g, lambda_DKB}
  switch (type) {
  case SplineType::Johnson:
    return {N, 0, N, 0, N, 0.0};
  case SplineType::Fischer: {
    // f = ~r^(kappa+1) for kappa > 0: the second spline (linear at r = 0)
    // is also dropped from f, which removes the spurious kappa = 1 state
    // (Fischer and Parpia)
    const auto imin_f = kappa > 0 ? 2ul : 1ul;
    return {N + 2, imin_f, N + 1, 1, N + 1, 0.0};
  }
  case SplineType::Hybrid:
    return {N + lk, lk, N + lk, lk, N + lk, 1.0};
  case SplineType::Shabaev:
    return {N + 2, 1, N + 1, 1, N + 1, 1.0};
  case SplineType::Derevianko:
    return {N + lk + 1, lk, N + lk, lk, N + lk, 1.0};
  }
  assert(false && "Unreachable: missing spline type in spline_params");
  return {};
}

//==============================================================================
SplineSpinorValues spline_spinor_values(double B, double dB, double d2B,
                                        int kappa, double r, double kb) {
  // f-like: f = B, g = kb (B' + kappa B/r)
  // g-like: g = B, f = kb (B' - kappa B/r)
  const auto k_r = kappa / r;
  SplineSpinorValues v{};
  v.f_f = B;
  v.g_f = kb * (dB + k_r * B);
  v.df_f = dB;
  v.dg_f = kb * (d2B + k_r * dB - k_r / r * B);
  v.f_g = kb * (dB - k_r * B);
  v.g_g = B;
  v.df_g = kb * (d2B - k_r * dB + k_r / r * B);
  v.dg_g = dB;
  return v;
}

//==============================================================================
SplineSpinors form_spline_basis(const int kappa, const std::size_t N_states,
                                const std::size_t k_spl, const double r0_spl,
                                const double rmax_spl,
                                std::shared_ptr<const Grid> rgrid,
                                const double alpha, SplineType type) {
  // Forms the underlying spline spinors for the given SplineType

  const auto params = spline_params(type, kappa, N_states);
  const auto kb = params.lambda_DKB * 0.5 * alpha;
  const auto num_f = params.num_f();
  const auto num_spinors = num_f + params.num_g();

  // uses sepperate B-splines for each partial wave! OK?

  BSpline bspl(params.n_spl, k_spl, r0_spl, rmax_spl,
               BSpline::KnotDistro::loglinear);

  SplineSpinors out;
  out.basis.resize(num_spinors, {0, kappa, rgrid});
  out.d_basis.resize(num_spinors, {0, kappa, rgrid});
  out.f_rmax.resize(num_spinors, 0.0);
  out.g_rmax.resize(num_spinors, 0.0);

  for (auto ir = 0ul; ir < rgrid->num_points(); ++ir) {
    const auto r = rgrid->r(ir);
    auto [i0, bij] = bspl.get_nonzero(r, 2);
    for (auto i = 0ul; i < bspl.K(); ++i) {

      const auto nB = i + i0; // which basis spline
      const auto v =
        spline_spinor_values(bij[i][0], bij[i][1], bij[i][2], kappa, r, kb);

      // "f-like" set: B_i in the large component
      if (nB >= params.imin_f && nB < params.imax_f) {
        auto &Sn = out.basis.at(nB - params.imin_f);
        auto &dSn = out.d_basis.at(nB - params.imin_f);
        Sn.f(ir) = v.f_f;
        Sn.g(ir) = v.g_f;
        dSn.f(ir) = v.df_f;
        dSn.g(ir) = v.dg_f;
      }

      // "g-like" set: B_i in the small component
      if (nB >= params.imin_g && nB < params.imax_g) {
        auto &Sn = out.basis.at(num_f + nB - params.imin_g);
        auto &dSn = out.d_basis.at(num_f + nB - params.imin_g);
        Sn.f(ir) = v.f_g;
        Sn.g(ir) = v.g_g;
        dSn.f(ir) = v.df_g;
        dSn.g(ir) = v.dg_g;
      }
    }
  }

  // Values exactly at rmax (generally between grid points), for the rmax
  // boundary terms
  const auto [i0, bij] = bspl.get_nonzero(rmax_spl, 2);
  for (auto i = 0ul; i < bspl.K(); ++i) {
    const auto nB = i + i0;
    const auto v = spline_spinor_values(bij[i][0], bij[i][1], bij[i][2], kappa,
                                        rmax_spl, kb);
    if (nB >= params.imin_f && nB < params.imax_f) {
      out.f_rmax.at(nB - params.imin_f) = v.f_f;
      out.g_rmax.at(nB - params.imin_f) = v.g_f;
    }
    if (nB >= params.imin_g && nB < params.imax_g) {
      out.f_rmax.at(num_f + nB - params.imin_g) = v.f_g;
      out.g_rmax.at(num_f + nB - params.imin_g) = v.g_g;
    }
  }

  return out;
}

//==============================================================================
std::pair<LinAlg::Matrix<double>, LinAlg::Matrix<double>>
fill_Hamiltonian_matrix(const SplineSpinors &splines, const Wavefunction &wf,
                        const bool correlationsQ, SplineType type) {
  const auto &spl_basis = splines.basis;
  const auto &d_basis = splines.d_basis;
  const auto n_spl = spl_basis.size();

  auto A_and_S = std::make_pair(LinAlg::Matrix<double>{n_spl, n_spl},
                                LinAlg::Matrix<double>{n_spl, n_spl});
  // auto &[Aij, Sij] = A_and_S; //XXX error w/ clang? Why?
  auto &Aij = A_and_S.first;
  auto &Sij = A_and_S.second;

  const auto excl_exch = wf.vHF() ? wf.vHF()->is_localQ() : true;
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
    add_JohnsonBoundary(&Aij, spl_basis.front().kappa(), wf.alpha(),
                        splines.f_rmax, splines.g_rmax);
  else if (type == SplineType::Hybrid)
    add_rmaxBoundary(&Aij, wf.alpha(), splines.f_rmax, splines.g_rmax);

  return A_and_S;
}

//==============================================================================
void add_JohnsonBoundary(LinAlg::Matrix<double> *pAij, const int kappa,
                         const double alpha, const std::vector<double> &f_rmax,
                         const std::vector<double> &g_rmax) {
  // W. R. Johnson, S. A. Blundell, J. Sapirstein, Phys. Rev. A 37, 307 (1988)
  // Second derivatives of the boundary term S', plus the surface term dropped
  // by H0ab (see header); f(0) is the coefficient of the first f-like spline.
  // The f(0)g(0) terms of the two cancel, so there is no r=0 cross term.
  auto &Aij = *pAij;
  const auto c = 1.0 / alpha;
  // r=0: (c/2) f(0)^2 [kappa<0], c^2 f(0)^2 [kappa>0]
  Aij[0][0] += (kappa < 0) ? c : 2.0 * c * c;
  add_rmaxBoundary(pAij, alpha, f_rmax, g_rmax);
}

//==============================================================================
void add_rmaxBoundary(LinAlg::Matrix<double> *pAij, const double alpha,
                      const std::vector<double> &f_rmax,
                      const std::vector<double> &g_rmax) {
  // r=R: (c/4)[f(R)^2 - g(R)^2] + (c/2) f(R) g(R), where f(R) = sum_i p_i
  // f_i(R) (likewise g) over the spline spinors: with kinetic balance, several
  // of the last splines contribute to both. The cross term restores the
  // surface term dropped by H0ab, so that the natural boundary condition is
  // f(R) = -g(R)
  auto &Aij = *pAij;
  assert(f_rmax.size() == Aij.rows() && g_rmax.size() == Aij.rows());
  const auto c = 1.0 / alpha;
  std::vector<std::size_t> nonzero_at_rmax;
  for (std::size_t i = 0; i < f_rmax.size(); ++i) {
    if (f_rmax[i] != 0.0 || g_rmax[i] != 0.0)
      nonzero_at_rmax.push_back(i);
  }
  for (const auto i : nonzero_at_rmax) {
    for (const auto j : nonzero_at_rmax) {
      const auto fi = f_rmax[i], gi = g_rmax[i];
      const auto fj = f_rmax[j], gj = g_rmax[j];
      Aij[i][j] += 0.5 * c * (fi * fj - gi * gj + fi * gj + gi * fj);
    }
  }
}

//==============================================================================
void expand_basis_orbitals(std::vector<DiracSpinor> *basis,
                           std::vector<DiracSpinor> *basis_positron,
                           const std::vector<DiracSpinor> &spl_basis,
                           const int kappa, const int max_n, int max_n_positron,
                           const LinAlg::Vector<double> &e_values,
                           const LinAlg::Matrix<double> &e_vectors,
                           const Wavefunction &wf) {
  // Expands the pseudo-spectrum basis in terms of B-spline basis and expansion
  // coeficient found from diagonalising the Hamiltonian over Bsplns

  const auto l = Angular::l_k(kappa);
  const auto min_n = l + 1;

  // Eigenvalues are in ascending order. The threshold -mc^2 lies inside the
  // gap (energies exclude the rest energy), so the negative-energy solutions
  // occupy the first num_negative slots.
  const auto neg_mc2 = -1.0 / (wf.alpha() * wf.alpha());
  std::size_t num_negative = 0;
  while (num_negative < e_values.rows() && e_values[num_negative] < neg_mc2) {
    ++num_negative;
  }

  // Electron states are taken from the bottom of the positive-energy branch
  // upwards, labelled n = l+1, l+2, ..., l+num_positive. Negative-energy
  // states are taken from just below the gap downwards (i.e., in reverse
  // order), so that on both branches the states nearest the gap are kept and
  // those far from it discarded. Their labels continue on from the
  // positive-energy branch, n = l+num_positive+1, ..., so {n, kappa} is unique
  // across the basis (there is a gap in n if not all electron states are
  // kept). max_n_positron is the number of negative-energy states kept (for
  // each kappa).
  const auto num_positive = static_cast<int>(e_values.rows() - num_negative);
  auto pqn_electron = min_n - 1;
  auto count_negative = 0;
  for (std::size_t j = 0; j < e_values.rows(); j++) {
    const auto negative_energy = j < num_negative;
    const auto i = negative_energy ? num_negative - 1 - j : j;
    const auto &en = e_values[i];
    const auto &pvec = e_vectors[i];

    if (negative_energy && ++count_negative > max_n_positron)
      continue;
    if (!negative_energy && ++pqn_electron > max_n)
      continue;
    const auto pqn =
      negative_energy ? l + num_positive + count_negative : pqn_electron;

    auto &Fi =
      negative_energy ?
        basis_positron->emplace_back(pqn, kappa, wf.grid_sptr(), true) :
        basis->emplace_back(pqn, kappa, wf.grid_sptr());
    Fi.en() = en;
    for (std::size_t ib = 0; ib < spl_basis.size(); ++ib) {
      Fi += pvec[ib] * spl_basis[ib];
    }
    // Enforce f > 0 at small r: find first significant amplitude
    {
      const auto &fv = Fi.f();
      const double fmax =
        *std::max_element(fv.cbegin(), fv.cend(), [](double a, double b) {
          return std::abs(a) < std::abs(b);
        });
      for (const auto fi : fv) {
        if (std::abs(fi) > 1.0e-4 * std::abs(fmax)) {
          if (fi < 0.0) {
            Fi *= -1.0;
          }
          break;
        }
      }
    }
    const auto Fnorm = Fi * Fi;
    Fi.normalise();
    if (std::abs(Fnorm - 1) > 1.0e-4) {
      std::cout << "Warning: Possible spurious state: " << Fi.shortSymbol()
                << " e=" << Fi.en() << ", Norm=" << Fnorm << " (" << Fnorm - 1.0
                << "). Killed\n";
      Fi *= 0.0;
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
      if (Fn.negativeEnergyStateQ())
        sum_p += term;
      else
        sum_el += term;
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

  assert(nDG >= 0 && nDG < 3); // better way?

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

  for (auto ki = 0ul; ki <= max_ki; ki++) {
    const auto kappa = Angular::kindex_to_kappa(ki);
    auto find_ka = [=](const auto &Fn) { return Fn.kappa() == kappa; };
    const auto pFa = std::find_if(basis.begin(), basis.end(), find_ka);
    // kappas missing from the basis (e.g., "20s5d" has no p states) skipped
    if (pFa == basis.end())
      continue;
    const auto &Fa = *pFa;
    // need to have l_n = la+1 terms, or sum doesn't work:
    if (Fa.l() == max_l)
      continue;
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
    printf("%-4s  %7.5f [%7.5f] %6.0e | %10.5f [%10.5f] %6.1e\n",
           Fa.shortSymbol().c_str(), sum1, 1.0, eps1, sumr2, sumr2_expect,
           epsr2);
  }
  return {eps1, epsr2};
}

} // namespace SplineBasis
