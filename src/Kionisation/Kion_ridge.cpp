#include "Kionisation/Kion_ridge.hpp"
#include "Angular/Wigner369j.hpp"
#include "DiracODE/FreeDirac.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/ChronoTimer.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/UnitConv_conversions.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/color.hpp"
#include "fmt/format.hpp"
#include "qip/Widgets.hpp"
#include "qip/omp.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace Kion {

//==============================================================================
MomentumOrbital momentum_orbital(const DiracSpinor &Fa, double en, double p_min,
                                 double p_max, std::size_t points_per_decade,
                                 double density_cut) {
  assert(p_min > 0.0 && p_max > p_min && points_per_decade > 0);
  MomentumOrbital orb;
  orb.kappa = Fa.kappa();
  orb.en = en;
  orb.num_electrons = (Fa.twoj() + 1) * Fa.occ_frac();

  const auto &gr = Fa.grid();
  const auto l = Fa.l();
  const auto l_tilde = Angular::l_tilde_k(Fa.kappa());
  // Radial range: the extent of the orbital
  const auto n_r = Fa.max_pt();
  const std::vector<double> r(gr.r().begin(), gr.r().begin() + long(n_r));

  // Transform in blocks of p (parallel), stopping once the density has been
  // negligible for a whole block beyond its peak. A node of the momentum
  // density (n > l+1 orbitals) is a single point, never a block.
  const std::size_t block = 32;
  double peak = 0.0;
  for (std::size_t ip0 = 0;; ip0 += block) {
    std::vector<double> p_block(block), F_block(block), G_block(block);
    for (std::size_t j = 0; j < block; ++j) {
      p_block[j] =
        p_min * std::pow(10.0, double(ip0 + j) / double(points_per_decade));
    }
    if (p_block.front() > p_max)
      break;
#pragma omp parallel for
    for (std::size_t j = 0; j < block; ++j) {
      const auto p = p_block[j];
      const auto jl = SphericalBessel::fillBesselVec_kr(l, p, r);
      const auto jl_tilde = SphericalBessel::fillBesselVec_kr(l_tilde, p, r);
      // f~(p) = 4 pi int f(r) j_l(pr) r dr: the Fourier transform of the
      // orbital, normalised to (2 pi)^3
      F_block[j] =
        4.0 * M_PI *
        NumCalc::integrate(gr.du(), 0, n_r, Fa.f(), jl, gr.r(), gr.drdu());
      G_block[j] = 4.0 * M_PI *
                   NumCalc::integrate(gr.du(), 0, n_r, Fa.g(), jl_tilde, gr.r(),
                                      gr.drdu());
    }
    bool all_negligible = true;
    for (std::size_t j = 0; j < block; ++j) {
      if (p_block[j] > p_max)
        break;
      orb.p.push_back(p_block[j]);
      orb.F.push_back(F_block[j]);
      orb.G.push_back(G_block[j]);
      const auto density = F_block[j] * F_block[j] + G_block[j] * G_block[j];
      peak = std::max(peak, density);
      if (density >= density_cut * peak) {
        all_negligible = false;
      }
    }
    if (all_negligible)
      break;
  }
  return orb;
}

//==============================================================================
std::array<double, 11> planewave_traces(double p, double F, double G, int kappa,
                                        double pf, double q, double ef,
                                        double alpha) {
  // Dirac representation, z along q, p in the xz plane (the returned
  // combinations are independent of the azimuth of p). Closed forms of
  //   Tr[(E_f + c alpha.p_f + beta m c^2) Gamma rho Gamma'^dag]
  // with rho = [[F^2, s F G sigma.n], [s F G sigma.n, G^2]] (2x2 blocks):
  // the momentum-space density matrix of the shell without its N_a/(8 pi)
  // prefactor; n = p-hat, p_f = p + q, s = sign(kappa) from the phases of
  // the momentum-space spinor. A = E_f + m c^2, B = E_f - m c^2 = e_f.
  const double c = 1.0 / alpha;
  const double mc2 = c * c;
  const double Ef = mc2 + ef;
  const double A = Ef + mc2;
  const double B = ef;

  // Energy conservation fixes the angle between p and q
  auto cos_t = (pf * pf - p * p - q * q) / (2.0 * p * q);
  cos_t = std::clamp(cos_t, -1.0, 1.0);
  const double sin_t = std::sqrt(1.0 - cos_t * cos_t);
  const double nz = cos_t;
  const double uz = p * cos_t + q;
  // n.p_f, and the components n_z u_z, n_x u_x
  const double n_u = p + q * cos_t;
  const double nzuz = nz * uz;
  const double nxux = p * sin_t * sin_t;

  const double s = kappa > 0 ? 1.0 : -1.0;
  const double FG = s * F * G;
  const double F2 = F * F;
  const double G2 = G * G;

  std::array<double, 11> T;
  // vector: temporal, longitudinal, transverse (E + M), temporal-longitudinal
  T[0] = 2.0 * A * F2 + 2.0 * B * G2 + 4.0 * c * FG * n_u;
  T[1] = 2.0 * A * G2 + 2.0 * B * F2 + 4.0 * c * FG * (nzuz - nxux);
  T[2] = 4.0 * A * G2 + 4.0 * B * F2 - 8.0 * c * FG * nzuz;
  T[3] = 2.0 * c * uz * (F2 + G2) + 4.0 * Ef * nz * FG;
  // axial: the same with the roles of F and G exchanged in the squares
  T[4] = 2.0 * A * G2 + 2.0 * B * F2 + 4.0 * c * FG * n_u;
  T[5] = 2.0 * A * F2 + 2.0 * B * G2 + 4.0 * c * FG * (nzuz - nxux);
  T[6] = 4.0 * A * F2 + 4.0 * B * G2 - 8.0 * c * FG * nzuz;
  T[7] = T[3];
  // vector-axial transverse interference: Im R^{12} (Z is twice this: the
  // antisymmetric part, Im(R^{12} - R^{21}))
  T[8] = 2.0 * c * uz * (F2 + G2) - 4.0 * Ef * nz * FG;
  // scalar, pseudoscalar
  T[9] = 2.0 * A * F2 + 2.0 * B * G2 - 4.0 * c * FG * n_u;
  T[10] = 2.0 * A * G2 + 2.0 * B * F2 - 4.0 * c * FG * n_u;
  return T;
}

//==============================================================================
std::array<double, 13> planewave_formFactors(const MomentumOrbital &orb,
                                             double E, double q, double alpha) {
  std::array<double, 13> K{};
  const auto ef = E + orb.en;
  if (ef <= 0.0 || q <= 0.0 || orb.p.size() < 2)
    return K;

  const double c = 1.0 / alpha;
  const double pf = std::sqrt(ef * (2.0 + alpha * alpha * ef));
  // Integration limits: |p_f - q| <= p <= p_f + q; nothing below the first
  // grid point (O(p^3)), or beyond the last (density negligible)
  const double p_lo = std::max(std::abs(pf - q), orb.p.front());
  const double p_hi = std::min(pf + q, orb.p.back());
  if (p_hi <= p_lo)
    return K;

  // F, G at an arbitrary p by linear interpolation (for the end points)
  const auto interpolate = [&](double p, std::size_t i1) {
    // i1: first grid index with p[i1] >= p
    const auto i0 = i1 - 1;
    const auto t = (p - orb.p[i0]) / (orb.p[i1] - orb.p[i0]);
    return std::pair{(1.0 - t) * orb.F[i0] + t * orb.F[i1],
                     (1.0 - t) * orb.G[i0] + t * orb.G[i1]};
  };

  // Points of the quadrature: the grid points inside [p_lo, p_hi], plus the
  // two limits themselves (interpolated) where they fall between grid points
  const auto it_lo = std::lower_bound(orb.p.begin(), orb.p.end(), p_lo);
  const auto it_hi = std::upper_bound(orb.p.begin(), orb.p.end(), p_hi);
  const auto i_lo = std::size_t(std::distance(orb.p.begin(), it_lo));
  const auto i_hi = std::size_t(std::distance(orb.p.begin(), it_hi));

  std::vector<double> ps, Fs, Gs;
  ps.reserve(i_hi - i_lo + 2);
  Fs.reserve(ps.capacity());
  Gs.reserve(ps.capacity());
  if (i_lo > 0 && orb.p[i_lo] > p_lo) {
    const auto [F, G] = interpolate(p_lo, i_lo);
    ps.push_back(p_lo);
    Fs.push_back(F);
    Gs.push_back(G);
  }
  for (auto i = i_lo; i < i_hi; ++i) {
    ps.push_back(orb.p[i]);
    Fs.push_back(orb.F[i]);
    Gs.push_back(orb.G[i]);
  }
  if (i_hi < orb.p.size() && orb.p[i_hi - 1] < p_hi) {
    const auto [F, G] = interpolate(p_hi, i_hi);
    ps.push_back(p_hi);
    Fs.push_back(F);
    Gs.push_back(G);
  }
  if (ps.size() < 2)
    return K;

  // int p dp T(p) = int p^2 T d(ln p): Simpson on the (uniform in ln p)
  // grid points, trapezoid on the partial end segments
  std::vector<std::array<double, 11>> Ts;
  Ts.reserve(ps.size());
  for (std::size_t i = 0; i < ps.size(); ++i) {
    Ts.push_back(
      planewave_traces(ps[i], Fs[i], Gs[i], orb.kappa, pf, q, ef, alpha));
  }
  std::array<double, 11> sum{};
  const auto add_trapezoid = [&](std::size_t i0, std::size_t i1) {
    const auto dx = std::log(ps[i1] / ps[i0]);
    for (std::size_t j = 0; j < sum.size(); ++j) {
      sum[j] +=
        0.5 * dx * (ps[i0] * ps[i0] * Ts[i0][j] + ps[i1] * ps[i1] * Ts[i1][j]);
    }
  };
  // ps: [p_lo if added] + grid points i_lo..i_hi-1 + [p_hi if added]
  const std::size_t n_grid = i_hi - i_lo;
  if (n_grid == 0) {
    // Both limits inside one grid interval: a single trapezoid
    add_trapezoid(0, 1);
  } else {
    const bool lo_end = ps.front() < orb.p[i_lo];
    const bool hi_end = ps.back() > orb.p[i_hi - 1];
    const std::size_t g0 = lo_end ? 1 : 0;
    const std::size_t g1 = g0 + n_grid - 1;
    if (lo_end)
      add_trapezoid(0, g0);
    if (hi_end)
      add_trapezoid(g1, g1 + 1);
    // Simpson over the grid points [g0, g1]: pairs of equal ln-p intervals;
    // an odd interval count leaves one trapezoid
    auto i = g0;
    for (; i + 2 <= g1; i += 2) {
      const auto dx = std::log(ps[i + 1] / ps[i]);
      for (std::size_t j = 0; j < sum.size(); ++j) {
        sum[j] += (dx / 3.0) * (ps[i] * ps[i] * Ts[i][j] +
                                4.0 * ps[i + 1] * ps[i + 1] * Ts[i + 1][j] +
                                ps[i + 2] * ps[i + 2] * Ts[i + 2][j]);
      }
    }
    if (i < g1)
      add_trapezoid(i, g1);
  }

  // R^{mu nu} = 1/(8 pi^2 c^2 q) int p dp Tr[...], with the density matrix
  // of the shell rho_a = (N_a / 8 pi) x (the matrix of planewave_traces)
  const double density_norm = orb.num_electrons / (8.0 * M_PI);
  const double pref = density_norm / (8.0 * M_PI * M_PI * q * c * c);
  // FormFactorSet order: V_T, V_E, V_M, V_L, X, A_T, A_E, A_M, A_L, Y, Z, S, P
  // The transverse response is not split between E and M: all in E
  K[0] = pref * sum[0];
  K[1] = pref * sum[2];
  K[2] = 0.0;
  K[3] = pref * sum[1];
  K[4] = -pref * sum[3];
  K[5] = pref * sum[4];
  K[6] = pref * sum[6];
  K[7] = 0.0;
  K[8] = pref * sum[5];
  K[9] = pref * sum[7];
  K[10] = 2.0 * pref * sum[8];
  K[11] = pref * sum[9];
  K[12] = pref * sum[10];
  return K;
}

//==============================================================================
double captured_fraction(const DiracSpinor &Fa, int Kmax, std::size_t iq,
                         const SphericalBessel::JL_table &jK_tab) {
  const auto &gr = Fa.grid();
  const auto n_r = Fa.max_pt();
  // Bare f^2 + g^2 (DiracSpinor::rho() carries the occupation factor, zero
  // for a model state)
  std::vector<double> density(n_r);
  for (std::size_t i = 0; i < n_r; ++i) {
    density[i] = Fa.f(i) * Fa.f(i) + Fa.g(i) * Fa.g(i);
  }
  double sum = 0.0;
  for (int K = 0; K <= Kmax; ++K) {
    const auto &jK = jK_tab.at(std::size_t(K), iq);
    sum += (2.0 * K + 1.0) *
           NumCalc::integrate(gr.du(), 0, n_r, density, jK, jK, gr.drdu());
  }
  return sum;
}

//==============================================================================
std::vector<FormFactorSet> calculate_ridge_correction(
  const HF::HartreeFock *vHF, const std::vector<DiracSpinor> &bound_states,
  double ec_min, double ec_max, const std::vector<double> &Egrid,
  const std::vector<double> &qgrid, const SphericalBessel::JL_table &jK_tab,
  int Kmax, bool vectorQ, bool axialQ, bool scalarQ, bool pseudoscalarQ,
  bool spatialQ, double ridge_eps) {
  IO::ChronoTimer timer("Ridge correction");

  const bool print = false;

  assert(vHF != nullptr);
  const auto &core = vHF->core();
  const auto n_core = core.size();
  assert(bound_states.size() == n_core);
  const auto E_steps = Egrid.size();
  const auto q_steps = qgrid.size();
  const auto alpha = vHF->alpha();
  const auto c = 1.0 / alpha;

  std::vector<FormFactorSet> dK(
    n_core, allocate_formFactors(E_steps, q_steps, vectorQ, axialQ, scalarQ,
                                 pseudoscalarQ, spatialQ));

  // Bound-bound leakage of the plane waves is confined by the triangle
  // rule on j, K <= j_a + j_b (= l_a + l_b + 1 for the operators with
  // C^K(kappa_b, -kappa_a)); below Kmax = 2 j_max the correction is not
  // defined: none is applied
  const auto twoj_max_core = DiracSpinor::max_tj(core);
  if (Kmax < twoj_max_core) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("ridge correction requires Kmax >= 2 j_max(core) = {} (have "
               "{}); no correction applied\n",
               twoj_max_core, Kmax);
    return dK;
  }

  // Where the correction is active: the q at which the multipoles above
  // Kmax carry more than ridge_eps of the norm. Their fraction grows with
  // q, so this is a threshold. The operators take qc as their frequency.
  std::vector<std::vector<std::pair<std::size_t, double>>> active_q(n_core);
  std::cout << "Ridge correction (plane-wave completion of K > " << Kmax
            << "):\n";
  for (std::size_t ia = 0; ia < n_core; ++ia) {
    for (std::size_t iq = 0; iq < q_steps; ++iq) {
      const auto captured =
        captured_fraction(bound_states[ia], Kmax, iq, jK_tab);
      if (1.0 - captured > ridge_eps) {
        active_q[ia].emplace_back(iq, qgrid[iq] * c);
      }
    }
    if (print) {
      if (active_q[ia].empty()) {
        fmt::print("  {:4s}: not needed (K <= {} complete to {:.0e})\n",
                   core[ia].shortSymbol(), Kmax, ridge_eps);
      } else {
        fmt::print("  {:4s}: active for q >= {:.3e} eV ({} of {} q points)\n",
                   core[ia].shortSymbol(),
                   qgrid[active_q[ia].front().first] *
                     UnitConv::Momentum_au_to_eV,
                   active_q[ia].size(), q_steps);
      }
    }
  }

  // Momentum-space orbitals (one at a time; the transform is parallel)
  // IO::ChronoTimer timer("");
  std::vector<MomentumOrbital> p_orbitals;
  p_orbitals.reserve(n_core);
  for (std::size_t ia = 0; ia < n_core; ++ia) {
    p_orbitals.push_back(active_q[ia].empty() ?
                           MomentumOrbital{} :
                           momentum_orbital(bound_states[ia], core[ia].en()));
  }

  // fmt::print("  momentum-space orbitals: {:.1f} s\n",
  //            timer.lap_reading_ms() / 1000.0);
  // timer.start();

  // Work list: every (orbital, E) with an active q and the ejected electron
  // energy within the limits
  std::vector<std::pair<std::size_t, std::size_t>> orbital_E_pairs;
  for (std::size_t ia = 0; ia < n_core; ++ia) {
    if (active_q[ia].empty())
      continue;
    for (std::size_t iE = 0; iE < E_steps; ++iE) {
      const auto ec = Egrid[iE] + core[ia].en();
      if (ec > ec_min && ec <= ec_max) {
        orbital_E_pairs.emplace_back(ia, iE);
      }
    }
  }
  if (orbital_E_pairs.empty()) {
    return dK;
  }

  // One operator set per thread (rank and frequency are set in place)
  const auto multipoles =
    multipole_operators(vHF->grid(), false, &jK_tab, vectorQ, axialQ, scalarQ,
                        pseudoscalarQ, spatialQ);
  const auto max_threads = std::size_t(omp_get_max_threads());
  std::vector<std::array<std::unique_ptr<DiracOperator::TensorOperator>, 10>>
    thread_multipoles(max_threads);
  for (auto &clones : thread_multipoles) {
    for (std::size_t i = 0; i < clones.size(); ++i) {
      clones[i] = multipoles[i] ? multipoles[i]->clone() : nullptr;
    }
  }

  // Squared factors (sums of squares): clamped at zero after the subtraction
  const std::array<bool, 13> squared{true,  true, true, true, false,
                                     true,  true, true, true, false,
                                     false, true, true};

  qip::ProgressBar bar(orbital_E_pairs.size());
#pragma omp parallel for schedule(dynamic)
  for (std::size_t it = 0; it < orbital_E_pairs.size(); ++it) {
    const auto [ia, iE] = orbital_E_pairs[it];
    const auto &Fa_hf = core[ia];
    const auto &Fa = bound_states[ia];
    const auto E = Egrid[iE];
    const auto ec = E + Fa_hf.en();

    // Far off the ridge in E, where |p_f - q| lies beyond the momentum grid
    // of the orbital, the plane-wave response is zero by construction (and
    // its partial sum negligible): nothing to complete
    const auto pf = std::sqrt(ec * (2.0 + alpha * alpha * ec));
    std::vector<std::pair<std::size_t, double>> q_columns;
    for (const auto &column : active_q[ia]) {
      if (std::abs(pf - qgrid[column.first]) <= p_orbitals[ia].p.back()) {
        q_columns.push_back(column);
      }
    }
    if (q_columns.empty()) {
      bar.update();
      continue;
    }

    // Plane-wave multipole sum, K = 0..Kmax, full continuum l range
    const auto [lc_min, lc_max] = continuum_l_range(Fa_hf, Kmax, std::nullopt);
    const auto free_waves =
      DiracODE::freeDirac(ec, lc_min, lc_max, Fa.grid_sptr(), alpha);
    auto PW = allocate_formFactors(1, q_steps, vectorQ, axialQ, scalarQ,
                                   pseudoscalarQ, spatialQ);
    auto &own_multipoles = thread_multipoles[std::size_t(omp_get_thread_num())];
    accumulate_multipole_sum(&PW, 0, Fa, Fa.occ_frac(), free_waves,
                             own_multipoles, 0, Kmax, q_columns);

    // All-K plane-wave response minus the plane-wave partial sum
    for (const auto &[iq, qc] : q_columns) {
      const auto IA =
        planewave_formFactors(p_orbitals[ia], E, qgrid[iq], alpha);
      for (std::size_t i = 0; i < 13; ++i) {
        if (dK[ia][i].empty())
          continue;
        auto value = IA[i] - PW[i](0, iq);
        // The transverse correction goes into E; M is subtracted there too
        if (i == 1 || i == 6) {
          value -= PW[i + 1](0, iq);
        } else if (i == 2 || i == 7) {
          value = 0.0;
        }
        if (squared[i]) {
          value = std::max(value, 0.0);
        }
        dK[ia][i](iE, iq) = value;
      }
    }
    bar.update();
  }
  // fmt::print("  plane-wave sums: {:.1f} s\n", timer.lap_reading_ms() / 1000.0);

  return dK;
}

} // namespace Kion
