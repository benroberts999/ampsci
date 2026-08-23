#include "Feynman.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracODE/include.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "MBPT/RadialMatrix.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "qip/omp.hpp"
#include <algorithm>
#include <cassert>
#include <memory>
#include <utility>
#include <vector>

namespace MBPT {

//==============================================================================
Feynman::Feynman(const HF::HartreeFock *vHF, std::size_t i0, std::size_t stride,
                 std::size_t size, const FeynmanOptions &options,
                 int n_min_core, bool include_G, bool verbose,
                 const std::string &ident, bool form_qpq)
  : m_HF(vHF),
    m_grid(vHF->grid_sptr()),
    m_i0(i0),
    m_stride(stride),
    m_subgrid_points(size),
    m_max_ki_core(std::size_t(2 * DiracSpinor::max_l(m_HF->core()))),
    m_max_ki(std::max(std::size_t(2 * options.max_l_internal), m_max_ki_core)),
    m_max_k(int(std::max(m_max_ki_core + 1, m_max_ki + 1))),
    m_min_core_n(n_min_core),
    m_include_G(include_G || vHF->vBreit() != nullptr),
    m_omre(options.omre),
    m_w0(options.w0),
    m_wratio(options.w_ratio),
    m_hole_particle(options.hole_particle == HoleParticle::include ||
                    options.hole_particle == HoleParticle::include_k0),
    m_include_higher_order_hp(options.hole_particle !=
                              HoleParticle::include_k0),
    m_screen_Coulomb(options.screening == Screening::include),
    m_Complex_green_method(options.complex_green) {

  form_w_quadrature(m_w0, m_wratio);

  if (verbose) {
    std::cout << "\nFeynman diagrams:\n";
    fmt::print("lmax = {} (internal lines)\n", Angular::kindex_to_l(m_max_ki));
    fmt::print("Including n ≥ {} in polarisation loops\n", m_min_core_n);
    if (m_screen_Coulomb) {
      std::cout << "Including all-orders Coulomb screening\n";
    }
    if (m_hole_particle) {
      std::cout << "Including hole-particle interaction: "
                << (m_include_higher_order_hp ? "(all k in hp)\n" :
                                                "(only k=0)\n");
    }
    std::cout << "Re(w) = " << m_omre << "\n";
    fmt::print("Im(w) : 0 + log grid [{}, {:.1f}], ratio={}, N={} "
               "(log-Simpson + tail)\n",
               m_w0, m_wgrid_points.back(), m_wratio, m_wgrid_points.size());
    std::cout << "Complex Green method: "
              << (m_Complex_green_method ? "direct (solve at complex energy)" :
                                           "Dyson (solve at Re(e), extend)")
              << "\n";
    if (m_include_G) {
      std::cout << "Inlcuding lower (g) components throughout\n";
    }
  }

  // Construct qk, and dri/drj
  form_qk();
  // Construct pa, Vx
  form_pa();
  form_vx();

  // Construct Q*Pi*Q (polarisation loop): the expensive step
  if (form_qpq) {
    form_qpiq(ident);
  }
}

//==============================================================================
void Feynman::form_qk() {
  // q_ij includes right-hand integration measure!
  // q_ij := (r_< / r_>)^k/r_> * dr_j

  // Numerically stable method (avoid calling pow)
  // r_> = rg = max(r_i, r_j)
  // r_< = rl = min(r_i, r_j)
  // q = rl^k / rg^{k+1}
  //   = (rl/rg)^k / rg
  // => q^k = q^{k-1} * (rl/rg)

  std::size_t size = std::size_t(m_max_k + 1);

  // radial coulomb operator without spin indices
  m_qk = std::vector<ComplexRMatrix>(
    size, {m_i0, m_stride, m_subgrid_points, m_HF->grid_sptr()});

  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto fi = m_qk.front().index_to_fullgrid(i);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto fj = m_qk.front().index_to_fullgrid(j);
      const auto drj = m_grid->drdu(fj) * m_grid->du() * double(m_stride);
      const auto rless = m_grid->r(std::min(fi, fj));
      const auto rgreater = m_grid->r(std::max(fi, fj));
      const auto ratio = rless / rgreater; // = r_< / r_>
      m_qk[0](i, j) = (1.0 / rgreater) * drj;
      for (int k = 1; k <= m_max_k; ++k) {
        const auto sk = std::size_t(k);
        m_qk[sk](i, j) = m_qk[sk - 1](i, j) * ratio;
      }
    }
  }
}

//==============================================================================
void Feynman::form_pa() {
  // Fill core |a><a|
  const auto &core = m_HF->core();
  constexpr auto one = std::complex<double>{1.0, 0.0};

  m_pa.resize(core.size(),
              {m_i0, m_stride, m_subgrid_points, m_include_G, m_grid});

  // Summed core projector, P = sum_a |a><a|, for each kappa
  m_Pcore.resize(m_max_ki + 1,
                 {m_i0, m_stride, m_subgrid_points, m_include_G, m_grid});

  for (auto ia = 0ul; ia < core.size(); ia++) {
    m_pa[ia] = green_single(core[ia], core[ia], one);
    // Projector must be idempotent on the sub-grid: with the large component
    // only, <a|a> = 1 - int(g^2) is short by up to a few % for the deepest
    // shells (and the deep poles would be incompletely removed), so
    // normalise by the discrete ff norm. With g included, <a|a> = 1 already
    auto norm = 1.0;
    if (!m_include_G) {
      norm = 0.0;
      for (auto i = 0ul; i < m_subgrid_points; ++i) {
        const auto si = m_pa[ia].index_to_fullgrid(i);
        norm += core[ia].f(si) * core[ia].f(si) * m_pa[ia].dr(i);
      }
    }
    m_Pcore[Angular::kappa_to_kindex(core[ia].kappa())] +=
      (1.0 / norm) * m_pa[ia];
  }
}

//==============================================================================
ComplexGMatrix Feynman::green_single(const DiracSpinor &ket,
                                     const DiracSpinor &bra,
                                     const std::complex<double> f) const {
  ComplexGMatrix Gmat(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);
  Gmat.add(ket, bra, f);
  return Gmat;
}

//==============================================================================
// new test way to form the radial exchange potential coordinate matrix using change of basis formula -- the complete basis that is used is the hydrogenic wave functions
void Feynman::form_vx() {

  const auto localQ = m_HF->is_localQ();
  if (localQ)
    return;

  m_Vx_kappa = std::vector<GMatrix>(
    m_max_ki + 1, {m_i0, m_stride, m_subgrid_points, m_include_G, m_grid});

  // Use basic H-like basis to express matrices in coordinate/orbital space
  // initialises the hydrogen wave function object
  Wavefunction wfH(m_grid, {"1", 0, "Ball"});
  wfH.set_HF();
  wfH.solve_core(false);

  const auto r0 = std::max(1.0e-4, 3.0 * m_grid->r0());
  const auto rmax = std::min(90.0, 0.9 * m_grid->rmax());
  const int number = 95;
  const auto order = 9;
  const auto max_l =
    std::max(Angular::kindex_to_l(m_max_ki), DiracSpinor::max_l(m_HF->core()));
  const auto l_string =
    AtomData::spectroscopic_notation.substr(0, std::size_t(max_l));
  const auto basis_string = l_string;

  wfH.formBasis(SplineBasis::Parameters(
    basis_string, number, order, r0, 0.0, rmax, basis_string,
    SplineBasis::SplineType::Derevianko, false, false));

  // constructs the exchange matrix as:
  // V(r1,r2) = \sum_n [Vx*F_n](r1) F_n^†(r2)
  for (const auto &Fn : wfH.basis()) {

    auto VxFn = m_HF->vexFa(Fn);
    if (m_HF->vBreit()) {
      VxFn += m_HF->VBr(Fn);
    }

    m_Vx_kappa[Fn.k_index()].add(VxFn, Fn, 1.0);
  }

  // includes both integration measures
  for (auto &Vx : m_Vx_kappa) {
    Vx.dri_in_place();
    Vx.drj_in_place();
  }
}

//==============================================================================
// // this function constructs the radial exchange coordinate matrix using the typical way that ampsci does it -- know this works since this is what the feynman method used by default
// void Feynman::form_vx_old() {
//   // this function forms the matrix form of the radial exchange matrix
//   // Vx = -|a>Q<a|

//   m_Vx_kappa =
//       std::vector<GMatrix>(std::size_t(m_max_ki + 1),
//                            {m_i0, m_stride, m_subgrid_points, true, m_grid});

//   for (int kapi = 0; kapi <= m_max_ki; ++kapi) {

//     const auto kappa = Angular::kindex_to_kappa(kapi);
//     const auto twojp1 = Angular::twoj_k(kappa) + 1;

//     for (int k = 0; k <= m_max_k; ++k) {
//       // GMatrix V_k(m_i0, m_stride, m_subgrid_points, false, m_grid);
//       ComplexGMatrix V_k(m_i0, m_stride, m_subgrid_points, true, m_grid);

//       const auto &core = m_HF->core();
//       for (auto ia = 0ul; ia < core.size(); ia++) {

//         const auto ck = Angular::Ck_kk(k, kappa, core.at(ia).kappa());
//         if (ck == 0.0)
//           continue;
//         const auto c_ang = -1.0 * ck * ck / double(twojp1);
//         // V_k.add(core[ia], core[ia], c_ang);
//         // V_k += c_ang * m_pa[ia].real();
//         V_k += c_ang * m_pa[ia];
//       }

//       const auto &qk = get_qk(k);

//       // This needs updating: Spinor matrix times radial matrix?
//       m_Vx_kappa[std::size_t(kapi)] += V_k.mult_elements_by(qk.dri()).real();
//     }
//   }
// }

//==============================================================================
ComplexGMatrix Feynman::green(int kappa, std::complex<double> en,
                              GreenStates states) const {
  if (states == GreenStates::core) {
    return green_core(kappa, en);
  } else if (states == GreenStates::excited) {
    return green_excited(kappa, en);
  }
  if (m_Complex_green_method) {
    // Hybrid: solve directly at complex energy where the kernel ridge
    // (width ~ 1/sqrt(2*Im(e))) is resolvable on the sub-grid; at larger
    // Im(e) the Dyson (resolvent) construction is exactly consistent with
    // the discrete radial sums, so use it there. Results are insensitive
    // to the switch point over Im(e) ~ 1 - 30.
    constexpr double u_switch = 10.0;
    if (std::abs(en.imag()) > u_switch) {
      return green_hf(kappa, en);
    }
    return green_hf_complex_dirac(kappa, en);
  }
  return green_hf(kappa, en);
}

//==============================================================================
ComplexGMatrix Feynman::green_complex_dirac(int kappa,
                                            std::complex<double> en) const {
  return green_hf_complex_dirac(kappa, en);
}

//==============================================================================
ComplexGMatrix Feynman::green_dyson(int kappa, std::complex<double> en) const {
  return green_hf(kappa, en);
}

//==============================================================================
ComplexGMatrix Feynman::green_hf(int kappa, std::complex<double> en,
                                 const DiracSpinor *Fc_hp) const {
  // Solve DE (no exchange), regular at 0, infinity ("pinf")
  DiracSpinor x0(0, kappa, m_grid);
  DiracSpinor xI(0, kappa, m_grid);

  const auto alpha = m_HF->alpha();
  const auto &Hmag = m_HF->Hmag(x0.l());

  using namespace qip::overloads;
  const auto vl = m_HF->vlocal(x0.l());

  DiracODE::regularAtOrigin(x0, en.real(), vl, Hmag, alpha);
  DiracODE::regularAtInfinity(xI, en.real(), vl, Hmag, alpha);

  // Evaluate Wronskian at ~65% of the way to pinf. Should be independent of r
  const auto pp = std::size_t(0.65 * double(xI.max_pt()));
  const auto w = (x0.f(pp) * xI.g(pp) - xI.f(pp) * x0.g(pp)) / alpha;

  // Get G0 (Green's function, without exchange):
  const auto g0 = construct_green_g0(x0, xI, w);

  // Don't include exchange if local!
  const auto localQ = m_HF->is_localQ();
  if (localQ) {
    if (en.imag() == 0.0)
      return g0.complex();
    std::complex<double> iw{0.0, en.imag()};
    return g0.complex() *
           ((iw * g0.complex().dri_in_place() + 1.0).invert_in_place());
  }

  // Include exchange (optionally, with hole-particle correction)
  auto Vx = get_Vx_kappa(kappa);
  if (Fc_hp != nullptr) {
    // Include hole-particle interaction (w/ [1-P]V[1-P]):
    Vx += calculate_Vhp(kappa, *Fc_hp);
  }

  // Include exchange, and imaginary energy part:
  if (en.imag() == 0.0) {
    // G = [1 - G0*Vx]^{-1} * G0 = -[G0*Vx-1]^{-1} * G0
    // nb: much faster to invert _before_ make complex!
    // (but, only if imag. part is zero)
    return -1.0 * g0.complex() *
           (((Vx * g0).complex() - 1.0).invert_in_place());
  }

  // G0 := G0(re{e}) - no exchange, only real part
  // G(e) = [1 + i*Im{e}*G0 - G0*Vx]^{-1} * G0
  // Note: differential dr is included in Vx (via Q)
  std::complex<double> iw{0.0, en.imag()};
  return g0.complex() *
         ((iw * g0.complex().dri_in_place() + 1.0 - (Vx * g0).complex())
            .invert_in_place());
}

//==============================================================================
ComplexGMatrix Feynman::green_hf_complex_dirac(int kappa,
                                               std::complex<double> en,
                                               const DiracSpinor *Fc_hp) const {
  if (en.imag() == 0.0) {
    // If en is real, don't need to use complex version:
    return green_hf(kappa, en, Fc_hp);
  }

  // Solve DE (no exchange), regular at 0, infinity ("pinf")
  DiracSpinor x0(0, kappa, m_grid);
  DiracSpinor xI(0, kappa, m_grid);
  DiracSpinor Ix0(0, kappa, m_grid);
  DiracSpinor IxI(0, kappa, m_grid);

  const auto alpha = m_HF->alpha();
  const auto &Hmag = m_HF->Hmag(x0.l());

  using namespace qip::overloads;
  const auto vl = m_HF->vlocal(x0.l());

  DiracODE::regularAtOrigin_C(x0, Ix0, en, vl, Hmag, alpha);
  DiracODE::regularAtInfinity_C(xI, IxI, en, vl, Hmag, alpha);

  // Evaluate Wronskian at ~65% of the way to pinf. Should be independent of r
  const auto pp = std::size_t(0.65 * double(xI.max_pt()));
  const auto I = std::complex{0.0, 1.0};
  // f -> f_r + I*f_i and same for g. (No conjugation: complex symmetric)
  const auto w = ((x0.f(pp) + I * Ix0.f(pp)) * (xI.g(pp) + I * IxI.g(pp)) -
                  (xI.f(pp) + I * IxI.f(pp)) * (x0.g(pp) + I * Ix0.g(pp))) /
                 alpha;

  // Get G0 (Green's function, without exchange):
  const auto g0 = construct_green_g0(x0, Ix0, xI, IxI, w);

  auto G = g0;

  // Include exchange (optionally, with hole-particle correction),
  // unless local (no exchange)
  if (!m_HF->is_localQ()) {
    auto Vx = get_Vx_kappa(kappa);
    if (Fc_hp != nullptr) {
      // Include hole-particle interaction (w/ [1-P]V[1-P]):
      Vx += calculate_Vhp(kappa, *Fc_hp);
    }
    // Include exchange using Dyson:
    G = -1.0 * ((g0 * Vx.complex() - 1.0).invert_in_place() * g0);
  }

  // Cell-average the continuum part of the sampled kernel. The core-pole
  // part is smooth and separable (not a ridge), so it is excluded from the
  // averaging and restored afterwards; averaging it distorts the delicate
  // monopole (k=0) cancellation, where the same-kappa pole sits only
  // |omre| from the contour.
  const auto Gcore = green_core(kappa, en);
  G -= Gcore;
  cell_average(&G, x0, Ix0, xI, IxI);
  G += Gcore;
  return G;
}

//==============================================================================
void Feynman::cell_average(ComplexGMatrix *G, const DiracSpinor &x0,
                           const DiracSpinor &Ix0, const DiracSpinor &xI,
                           const DiracSpinor &IxI) const {
  // Kernel:
  //   G(r1,r2) = chi0(r<) * chiI(r>) / W
  // Within one sub-grid cell (width h) each solution is locally exponential,
  // with complex rate p (real part: growth/decay, imag part: oscillation):
  //   chi0(r_i + t) ~ chi0(r_i) * exp(+p*t)
  //   chiI(r_i + t) ~ chiI(r_i) * exp(-p*t)
  // so G falls off like exp(-p*|r1-r2|) across the diagonal, with width 1/p
  // often below the sub-grid spacing h. Downstream double integrations
  // sample each dr*dr cell at a single point, which overestimates these
  // unresolved near-diagonal cells.
  //
  // So, we replace the sampled value of each cell by the average of the local
  // exponential over the cell. With x = p*h for the cell, the average of
  // exp(p*t) over t in [-h/2, h/2], relative to its central value, is
  //   s(x) = (1/h) * Int_{-h/2}^{+h/2} exp(p*t) dt = sinh(x/2) / (x/2)
  // Off-diagonal cells (i,j) factorise into one average per solution:
  //   G(i,j) *= s(x_i) * s(x_j)
  // On the diagonal cell (i,i), r< and r> cross inside the cell, so average
  // the ridge profile over the 2D cell instead:
  //   f(x) = (1/h^2) * Int Int exp(-p*|t1 - t2|) dt1 dt2
  //        = 2*(x - 1 + exp(-x)) / x^2
  //   G(i,i) *= f(x_i)
  //
  // The rate is measured from the solutions' own ratio across one cell:
  //   x = [ ln(chi0_{i+1}/chi0_i) - ln(chiI_{i+1}/chiI_i) ] / 2
  // The complex log keeps sign and phase: Im(x) averages the local
  // oscillation; Re(x) is clamped >= 0 so the model always decays away
  // from the diagonal. All factors -> 1 as x -> 0, so cells that resolve
  // the kernel are unchanged.
  const auto I = std::complex{0.0, 1.0};
  std::vector<std::complex<double>> sfac(m_subgrid_points, {1.0, 0.0});
  std::vector<std::complex<double>> ffac(m_subgrid_points, {1.0, 0.0});
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = G->index_to_fullgrid(i);
    const bool last = (i + 1 == m_subgrid_points);
    const auto sj = G->index_to_fullgrid(last ? i - 1 : i + 1);
    const auto c0_i = x0.f(si) + I * Ix0.f(si);
    const auto c0_j = x0.f(sj) + I * Ix0.f(sj);
    const auto cI_i = xI.f(si) + I * IxI.f(si);
    const auto cI_j = xI.f(sj) + I * IxI.f(sj);
    // beyond the practical infinity the solutions are zero: leave cell as-is
    if (std::abs(c0_i) == 0.0 || std::abs(c0_j) == 0.0 ||
        std::abs(cI_i) == 0.0 || std::abs(cI_j) == 0.0) {
      continue;
    }
    auto xx = 0.5 * (std::log(c0_j / c0_i) - std::log(cI_j / cI_i));
    if (last) {
      xx = -xx;
    }
    if (xx.real() < 0.0) {
      xx = std::complex{-xx.real(), xx.imag()};
    }
    if (std::abs(xx) < 1.0e-3) {
      sfac[i] = 1.0 + xx * xx / 24.0;
      ffac[i] = 1.0 - xx / 3.0 + xx * xx / 12.0;
    } else {
      sfac[i] = std::sinh(0.5 * xx) / (0.5 * xx);
      ffac[i] = 2.0 * (xx - 1.0 + std::exp(-xx)) / (xx * xx);
    }
  }
  const bool include_g = G->includes_g();
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto cell = (i == j) ? ffac[i] : sfac[i] * sfac[j];
      G->ff(i, j) *= cell;
      if (include_g) {
        G->fg(i, j) *= cell;
        G->gf(i, j) *= cell;
        G->gg(i, j) *= cell;
      }
    }
  }
}

//==============================================================================
ComplexGMatrix
Feynman::green_basis(int kappa, std::complex<double> en,
                     const std::vector<DiracSpinor> &basis) const {
  ComplexGMatrix Gmat(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);
  for (const auto &n : basis) {
    if (n.kappa() != kappa)
      continue;
    const auto inv_de = 1.0 / (en - std::complex<double>{n.en()});
    Gmat.add(n, n, inv_de);
  }
  return Gmat;
}

//==============================================================================
ComplexGMatrix Feynman::green_core(int kappa, std::complex<double> en) const {
  // G_core = \sum_a |a><a|/(e_r + i*e_i-ea), for all a with a.kappa() = k
  ComplexGMatrix Gcore(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

  // loop over HF core, not Sigma core (used in subtraction to get
  // G^excited)
  const auto &core = m_HF->core();
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    if (a.kappa() != kappa)
      continue;
    const auto inv_de = 1.0 / (en - std::complex<double>{a.en()});
    Gcore += inv_de * m_pa[ia]; // Pa = |a><a|
  }
  return Gcore;
}

//==============================================================================
ComplexGMatrix Feynman::green_excited(int kappa, std::complex<double> en,
                                      const DiracSpinor *Fc_hp) const {
  ComplexGMatrix Gk(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

  // Subtract core states, by forceing Gk to be orthogonal to core:
  // Gk -> Gk - \sum_a|a><a|G

  const auto g0 = m_Complex_green_method ?
                    green_hf_complex_dirac(kappa, en, Fc_hp) :
                    green_hf(kappa, en, Fc_hp);
  // nb: can also subtract of core part, but doesn't seem to make a difference
  // return orthogonalise_wrt_core(g0 - green_core(kappa, en), kappa);
  return orthogonalise_wrt_core(g0, kappa);
}

//==============================================================================
ComplexGMatrix Feynman::orthogonalise_wrt_core(const ComplexGMatrix &g_in,
                                               int kappa) const {
  // Project the core states out of G from both sides:
  //   G -> (1 - P) G (1 - P),  P = sum_a |a><a| dr  (sub-grid brakets).
  // One-sided projection (1-P)G leaves the right-side pole residual
  // G|a><a|, which survives discretisation and dominates the monopole
  // (k=0) channel, where the sandwich Fa^dag G Fa picks it up linearly.
  // The symmetric projection removes both sides, and remains valid for
  // the hole-particle-dressed G: core states are still eigenstates under
  // the [1-P]V[1-P] dressing, but the numerical G carries core-excited
  // cross terms that one-sided subtraction leaves behind.
  const auto &P = get_Pcore(kappa);
  return (-1.0 * P.drj() + 1.0) * g_in * (-1.0 * P.dri() + 1.0);
}

//==============================================================================
GMatrix Feynman::calculate_Vhp(int kappa, const DiracSpinor &Fc) const {
  // Hole-particle interaction, extra potential

  GMatrix V0(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

  const auto y0cc = Coulomb::yk_ab(0, Fc, Fc);
  // Local potential: same on both spinor components
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = V0.index_to_fullgrid(i);
    V0.ff(i, i) = -y0cc[si];
    if (m_include_G) {
      V0.gg(i, i) = -y0cc[si];
    }
  }

  // Higher-k hole-particle
  if (m_include_higher_order_hp) {
    const auto tjcp1 = Fc.twojp1();
    for (int k = 2; k <= m_max_k; k += 2) {
      const auto ykcc = Coulomb::yk_ab(k, Fc, Fc);
      const auto ck = Angular::Ck_kk(k, Fc.kappa(), Fc.kappa());
      const auto coef = (ck / tjcp1) * (ck / tjcp1);
      for (auto i = 0ul; i < m_subgrid_points; ++i) {
        const auto si = V0.index_to_fullgrid(i);
        V0.ff(i, i) += coef * ykcc[si];
        if (m_include_G) {
          V0.gg(i, i) += coef * ykcc[si];
        }
      }
    }
  }

  V0.drj_in_place();

  // (1-P), with P = sum_a |a><a| drj
  const auto OneNegPc = -1.0 * get_Pcore(kappa).real().drj_in_place() + 1.0;
  return OneNegPc * V0 * OneNegPc;
}

//==============================================================================
ComplexGMatrix Feynman::green_to_complex(const ComplexGMatrix &Gr,
                                         double om_imag) const {
  // Given G(wr) and wi, returns G(wr+i*wi)
  // G(w) =  G(re(w)+im(w)) ;  Gr = G(re(w)), G = G(w),   im(w) = wi
  // G = [1 + i*wi*Gr]^-1 * Gr
  const std::complex<double> iw{0.0, om_imag};
  return ((iw * Gr).drj_in_place() + 1.0).invert_in_place() * Gr;
}

//==============================================================================
GMatrix Feynman::construct_green_g0(const DiracSpinor &x0,
                                    const DiracSpinor &xI,
                                    const double w) const {
  // Takes sub-grid into account; ket,bra are on full grid, G on sub-grid
  // G(r1,r2) = x0(rmin)*xI(imax)/w

  // G(ri,rj) = xI(ri) x0^†(rj) when ri > rj
  // G(ri,rj) = x0(ri) xI^†(rj) when ri < rj
  // ==>, when rj < ri
  // G(rj,ri) = x0(rj) xI^†(ri)
  //          = (xI(ri) x0^†(rj))^†
  //          = G(ri,rj)^†
  // ^† refers to spinor space

  GMatrix g0I(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

  const auto winv = 1.0 / w;

  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = g0I.index_to_fullgrid(i);
    for (auto j = 0ul; j <= i; ++j) {
      // i >= j
      const auto sj = g0I.index_to_fullgrid(j);

      g0I.ff(i, j) = xI.f(si) * x0.f(sj) * winv;
      g0I.ff(j, i) = g0I.ff(i, j);

      if (!m_include_G)
        continue;

      g0I.fg(i, j) = xI.f(si) * x0.g(sj) * winv;
      g0I.gf(j, i) = g0I.fg(i, j);

      g0I.gf(i, j) = xI.g(si) * x0.f(sj) * winv;
      g0I.fg(j, i) = g0I.gf(i, j);

      g0I.gg(i, j) = xI.g(si) * x0.g(sj) * winv;
      g0I.gg(j, i) = g0I.gg(i, j);
    }
  }

  return g0I;
}

//==============================================================================
ComplexGMatrix Feynman::construct_green_g0(const DiracSpinor &x0,
                                           const DiracSpinor &Ix0,
                                           const DiracSpinor &xI,
                                           const DiracSpinor &IxI,
                                           const std::complex<double> w) const {
  // Takes sub-grid into account; ket,bra are on full grid, G on sub-grid
  // G(r1,r2) = x0(rmin)*xI(imax)/w
  // Same structure as real version: for i >= j, G(ri,rj) = xI(ri) x0^T(rj).
  // At complex energy G is complex-symmetric (transpose), NOT Hermitian:
  // no complex conjugation anywhere.
  ComplexGMatrix g0I(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

  const auto winv = 1.0 / w;

  const auto I = std::complex{0.0, 1.0};

  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = g0I.index_to_fullgrid(i);
    for (auto j = 0ul; j <= i; ++j) {
      const auto sj = g0I.index_to_fullgrid(j);

      // j <= i
      const auto x0f = x0.f(sj) + I * Ix0.f(sj);
      const auto x0g = x0.g(sj) + I * Ix0.g(sj);
      const auto xIf = xI.f(si) + I * IxI.f(si);
      const auto xIg = xI.g(si) + I * IxI.g(si);

      g0I.ff(i, j) = xIf * x0f * winv;
      g0I.ff(j, i) = g0I.ff(i, j);

      if (!m_include_G)
        continue;

      g0I.fg(i, j) = xIf * x0g * winv;
      g0I.gf(j, i) = g0I.fg(i, j);

      g0I.gf(i, j) = xIg * x0f * winv;
      g0I.fg(j, i) = g0I.gf(i, j);

      g0I.gg(i, j) = xIg * x0g * winv;
      g0I.gg(j, i) = g0I.gg(i, j);
    }
  }

  return g0I;
}

//==============================================================================
std::vector<ComplexRMatrix>
Feynman::polarisation_each_k(std::complex<double> omega,
                             bool hole_particle) const {

  // polarisation operator is ~ Fa^dag * [Gex(ea + w) + Gex(ea - w)] * Fa
  // The Green's functions, and the sandwich Fa^dag * Gex * Fa, are independent
  // of the multipolarity k: form them once per (a, kappa_n), and accumulate
  // into each pi_k with its angular factor.

  std::vector<ComplexRMatrix> pi_k(
    std::size_t(m_max_k + 1),
    ComplexRMatrix{m_i0, m_stride, m_subgrid_points, m_grid});

  const auto Iunit = std::complex<double>{0.0, 1.0};
  const auto &core = m_HF->core();
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &Fa = core[ia];
    if (Fa.n() < m_min_core_n)
      continue;

    const auto ea_minus_w = std::complex<double>{Fa.en()} - omega;
    const auto ea_plus_w = std::complex<double>{Fa.en()} + omega;

    // not m_hole_particle, as need both for "screen only"
    const auto *Fa_hp = hole_particle ? &Fa : nullptr;

    for (auto in = 0ul; in <= m_max_ki; ++in) {
      const auto kn = Angular::kindex_to_kappa(in);

      // Angular factor for each multipole k; skip Green's fns if all zero
      std::vector<std::pair<std::size_t, double>> k_cang;
      for (int k = 0; k <= m_max_k; ++k) {
        const auto ck_an = Angular::Ck_kk(k, Fa.kappa(), kn);
        if (ck_an != 0.0) {
          k_cang.emplace_back(std::size_t(k),
                              ck_an * ck_an / double(2 * k + 1));
        }
      }
      if (k_cang.empty())
        continue;

      const ComplexGMatrix Gx_pm = green_excited(kn, ea_minus_w, Fa_hp) +
                                   green_excited(kn, ea_plus_w, Fa_hp);

      // sandwich: Sa ~ Fa^dag(r1)[Gex(r1,r2,ea-w) + Gex(r1,r2,ea+w)]Fa(r2)
      // pi is symmetric in r1,r2: fill lower half, then mirror
      ComplexRMatrix Sa(m_i0, m_stride, m_subgrid_points, m_grid);
      const bool include_g = Gx_pm.includes_g();
      for (auto i = 0ul; i < m_subgrid_points; ++i) {
        const auto si = Gx_pm.index_to_fullgrid(i);
        for (auto j = 0ul; j <= i; ++j) {
          const auto sj = Gx_pm.index_to_fullgrid(j);
          Sa(i, j) = Fa.f(si) * Gx_pm.ff(i, j) * Fa.f(sj);
          if (include_g) {
            Sa(i, j) += Fa.g(si) * Gx_pm.gf(i, j) * Fa.f(sj) +
                        Fa.f(si) * Gx_pm.fg(i, j) * Fa.g(sj) +
                        Fa.g(si) * Gx_pm.gg(i, j) * Fa.g(sj);
          }
          Sa(j, i) = Sa(i, j);
        }
      }

      for (const auto &[k, c_ang] : k_cang) {
        pi_k[k] += std::complex<double>{c_ang} * Sa;
      }
    }
  }

  for (auto &pik : pi_k) {
    pik *= Iunit;
  }
  return pi_k;
}

//==============================================================================
void Feynman::form_w_quadrature(double w0, double wratio) {
  // Quadrature for int_0^infty F(u) du, u = Im(w), in linear measure:
  //   int_0^infty F(u) du =~ sum_i W_i F(u_i).
  // The integrand is finite at u = 0 (generically its maximum, and flat:
  // F(u) - F(0) is O(u^2)), so u = 0 is included as an explicit point:
  //  - [0, w0] panel: trapezoid, (w0/2)[F(0) + F(w0)]
  //  - [w0, wmax]: composite Simpson's rule in t = ln(u) (uniform in t,
  //    jacobian du/dt = u; endpoint coefficients 1/3, requires odd N_log)
  //  - tail: F ~ 1/u^3 at large u [g ~ 1/u, qpiq ~ 1/u^2],
  //    so int_wmax^infty =~ (wmax/2) F(wmax)

  // Maximum Im(w): the polarisation loop is O(1) for u up to ~|e_core|;
  // beyond a few |e_core| the u^-3 tail correction handles the remainder
  // (with the tail included, the result is insensitive to the exact wmax)
  auto wmax_core = 30.0;
  const auto &core = m_HF->core();
  for (const auto &Fc : core) {
    if (Fc.n() < m_min_core_n)
      continue;
    if (std::abs(Fc.en()) > wmax_core)
      wmax_core = std::abs(Fc.en());
  }
  const auto wmax_t = 2.0 * wratio * wmax_core;

  // Solve wmax < w0 * ratio^{N-1} for N; odd number of log points
  std::size_t wsteps =
    std::size_t(std::log(wratio * wmax_t / w0) / std::log(wratio)) + 1;
  if (wsteps % 2 == 0)
    ++wsteps;

  const auto dt = std::log(wratio);

  m_wgrid_points.clear();
  m_wgrid_weights.clear();
  m_wgrid_points.reserve(wsteps + 1);
  m_wgrid_weights.reserve(wsteps + 1);

  // u = 0 point: half of the [0, w0] trapezoid panel
  m_wgrid_points.push_back(0.0);
  m_wgrid_weights.push_back(0.5 * w0);

  for (std::size_t i = 0; i < wsteps; ++i) {
    const auto u = w0 * std::pow(wratio, double(i));
    const auto simpson = (i == 0 || i == wsteps - 1) ? 1.0 / 3.0 :
                         (i % 2 == 1)                ? 4.0 / 3.0 :
                                                       2.0 / 3.0;
    auto W = simpson * dt * u;
    if (i == 0) {
      // other half of the [0, w0] trapezoid panel
      W += 0.5 * w0;
    }
    if (i == wsteps - 1) {
      // u^-3 tail beyond wmax
      W += 0.5 * u;
    }
    m_wgrid_points.push_back(u);
    m_wgrid_weights.push_back(W);
  }
}

//==============================================================================
double best_omre(const std::vector<DiracSpinor> &core,
                 const std::vector<DiracSpinor> &valence, bool print) {
  // Real part of frequency Re{omega} = omre:
  // should sit as far as possible from the poles.
  // The poles:
  //  - Green's function, g(e_v + omre): bound-state poles. The core
  //    ones require omre <= e_core_max - e_v; the excited ones require
  //    omre >= 0, except the few excited states BELOW e_v, which give
  //    omre = e_m - e_v (small, near zero; only for non-lowest valence).
  //  - Polarisation loop, gex(e_a +/- omre): true (excited) poles all lie
  //    at omre <= -Delta. The numerical core subtraction in gex is
  //    imperfect, leaving weak fictitious poles at the core-core
  //    differences, omre = -|e_b - e_a|.
  // So with Delta = e_lowest_excited - e_core_max, the window (-Delta, 0)
  // is free of true poles for the lowest valence state; inside it are only
  // the fictitious core-difference poles and the small e_m - e_v
  // differences.
  // Returns the midpoint of the largest pole-free gap. One omre serves all
  // valence states (the window is set by the lowest one; the others only
  // add near-zero poles, which the gap search avoids).

  constexpr double omre_default = -0.3;
  // Assumed lowest-valence energy, if no valence states are given
  constexpr double e_valence_typical = -0.1;

  if (core.empty())
    return omre_default;

  // Search window (-Delta, 0), Delta = lowest excited - highest core
  const auto e_core_max = DiracSpinor::max_En(core);
  const auto e_v_min =
    valence.empty() ? e_valence_typical : DiracSpinor::min_En(valence);
  const auto Delta = e_v_min - e_core_max;
  if (Delta <= 0.0)
    return omre_default;

  // Keep only poles inside the window; sort; drop duplicates
  const auto tidy = [=](std::vector<double> &list) {
    list.erase(std::remove_if(list.begin(), list.end(),
                              [=](double w) { return w < -Delta || w > 0.0; }),
               list.end());
    std::sort(list.begin(), list.end());
    list.erase(
      std::unique(list.begin(), list.end(),
                  [](double a, double b) { return std::abs(a - b) < 1.0e-9; }),
      list.end());
  };

  // True poles: the window ends, and the Green's-function (valence-line)
  // poles omre = e_m - e_v (excited states below e_v)
  std::vector<double> true_poles{-Delta, 0.0};
  for (const auto &Fv : valence) {
    for (const auto &Fm : valence) {
      true_poles.push_back(Fm.en() - Fv.en());
    }
  }
  tidy(true_poles);

  // Fictitious poles (imperfect core subtraction in the loop gex) at
  // (minus) core-core differences
  std::vector<double> fict_poles;
  for (const auto &Fa : core) {
    for (const auto &Fb : core) {
      // each unordered pair once; skip self-pairs
      if (Fb >= Fa)
        continue;
      fict_poles.push_back(-std::abs(Fa.en() - Fb.en()));
    }
  }
  tidy(fict_poles);

  // Combined list for the gap search
  auto poles = true_poles;
  poles.insert(poles.end(), fict_poles.cbegin(), fict_poles.cend());
  tidy(poles);

  // Best omre: midpoint of the largest gap between consecutive poles
  auto best = -0.5 * Delta;
  auto best_gap = 0.0;
  for (auto i = 1ul; i < poles.size(); ++i) {
    const auto gap = poles[i] - poles[i - 1];
    if (gap > best_gap) {
      best_gap = gap;
      best = 0.5 * (poles[i] + poles[i - 1]);
    }
  }

  bool print_all_true_poles = false;
  if (print) {
    fmt::print("\nFeynman contour:\n");
    fmt::print("Delta = {:.4f}\n", Delta);
    if (print_all_true_poles) {
      fmt::print("True poles in window (valence-valence):\n  ");
      for (const auto w : true_poles) {
        fmt::print("{:.4f} ", w);
      }
    } else {
      if (true_poles.size() >= 2) {
        // True poles includes one *at* Delta, rest between 0.0 and next
        fmt::print("True poles in window (valence-valence):\n  ");
        fmt::print("{} poles between {:.4f} and 0.0", true_poles.size() - 1,
                   true_poles.at(1));
      }
    }
    fmt::print("\nFictitious poles in window (core-core):\n  ");
    for (const auto w : fict_poles) {
      fmt::print("{:.4f} ", w);
    }
    fmt::print("\nBest omre = {:.4f}\n", best);
    fmt::print("Distance to nearest pole: {:.4f}\n", 0.5 * best_gap);
  }

  return best;
}

//==============================================================================
std::string Feynman::qpiq_filename(const std::string &ident) const {
  const auto prefix = ident.substr(0, ident.find('.'));
  if (prefix == "" || prefix == "false")
    return "";
  return prefix + ".qpq" + (m_hole_particle ? "h" : "") +
         (m_screen_Coulomb ? "s" : "") +
         (m_HF->vBreit() == nullptr ? "" : "b") + (m_include_G ? "g" : "") +
         std::to_string(m_min_core_n) + ".abf";
}

//==============================================================================
bool Feynman::read_qpiq(const std::string &ident) {
  const auto fname = qpiq_filename(ident);
  return fname.empty() ? false : readwrite_qpiq(IO::FRW::read, fname);
}

//==============================================================================
bool Feynman::write_qpiq(const std::string &ident) {
  const auto fname = qpiq_filename(ident);
  return fname.empty() ? false : readwrite_qpiq(IO::FRW::write, fname);
}

//==============================================================================
bool Feynman::readwrite_qpiq(IO::FRW::RoW rw, const std::string &fname) {

  const auto readQ = rw == IO::FRW::read;

  if (!readQ && fname == "")
    return false;

  if (readQ && !IO::FRW::file_exists(fname))
    return false;

  // For comparing floats:
  constexpr double eps = 1.0e-10;
  auto fequal = [](double a, double b) { return std::abs(a - b) <= eps; };

  std::fstream iofs;
  IO::FRW::open_binary(iofs, fname, rw);

  // Check screening / hole-particle (should be different filename)
  bool t_hp{m_hole_particle}, t_sc{m_screen_Coulomb},
    t_hohp{m_include_higher_order_hp}, t_cgm{m_Complex_green_method},
    t_iG{m_include_G};
  rw_binary(iofs, rw, t_hp, t_sc, t_hohp, t_cgm, t_iG);
  if (t_hp != m_hole_particle || t_sc != m_screen_Coulomb ||
      t_hohp != m_include_higher_order_hp || t_cgm != m_Complex_green_method ||
      t_iG != m_include_G)
    return false;

  // Other parameters that make a difference
  std::size_t t_max_ki{m_max_ki}, t_max_ki_core{m_max_ki_core};
  int t_min_core_n{m_min_core_n};
  rw_binary(iofs, rw, t_max_ki, t_min_core_n, t_max_ki_core);
  if (t_max_ki != m_max_ki || t_min_core_n != m_min_core_n ||
      t_max_ki_core != m_max_ki_core)
    return false;

  // Check HartreeFock things
  const bool BreitQ = m_HF->vBreit() != nullptr;
  bool t_Breit{BreitQ};
  double t_alpha{m_HF->alpha()};
  int t_num_core_e{m_HF->num_core_electrons()};
  rw_binary(iofs, rw, t_Breit, t_alpha, t_num_core_e);
  if (t_Breit != BreitQ || !fequal(t_alpha, m_HF->alpha()) ||
      t_num_core_e != m_HF->num_core_electrons())
    return false;

  // Check subgrid:
  std::size_t t_i0{m_i0}, t_stride{m_stride},
    t_subgrid_points{m_subgrid_points}, t_size(m_grid->size());
  rw_binary(iofs, rw, t_i0, t_stride, t_subgrid_points);
  if (t_i0 != m_i0 || t_stride != m_stride ||
      t_subgrid_points != m_subgrid_points || t_size != m_grid->size())
    return false;

  // Check regular grid:
  double t_r0{m_grid->r0()}, t_rmax{m_grid->rmax()},
    t_loglin{m_grid->loglin_b()}, t_du{m_grid->du()};
  rw_binary(iofs, rw, t_r0, t_rmax, t_loglin, t_du);
  if (!fequal(t_r0, m_grid->r0()) || !fequal(t_rmax, m_grid->rmax()) ||
      !fequal(t_loglin, m_grid->loglin_b()) || !fequal(t_du, m_grid->du()))
    return false;

  // Check frequency quadrature (points define the stored qpiq)
  std::size_t twsize{m_wgrid_points.size()};
  rw_binary(iofs, rw, twsize);
  if (twsize != m_wgrid_points.size())
    return false;
  for (const auto &u : m_wgrid_points) {
    double t_u{u};
    rw_binary(iofs, rw, t_u);
    if (!fequal(t_u, u))
      return false;
  }

  int t_max_k{m_max_k};
  double t_omre{m_omre};
  rw_binary(iofs, rw, t_max_k, t_omre);
  if (t_max_k != m_max_k || !fequal(t_omre, m_omre))
    return false;

  // Now, do actual read/write of data:

  const auto num_ks = std::size_t(m_max_k + 1);
  const auto num_ws = m_wgrid_points.size();
  if (readQ) {
    m_qpiq_wk.resize(num_ws, num_ks,
                     ComplexRMatrix{m_i0, m_stride, m_subgrid_points, m_grid});
  }

  for (auto iw = 0ul; iw < num_ws; ++iw) {
    for (auto k = 0ul; k < num_ks; ++k) {
      for (std::size_t i = 0; i < m_subgrid_points; ++i) {
        for (std::size_t j = 0; j < m_subgrid_points; ++j) {
          double re = m_qpiq_wk[iw][k](i, j).real();
          double im = m_qpiq_wk[iw][k](i, j).imag();
          rw_binary(iofs, rw, re, im);
          if (readQ)
            m_qpiq_wk[iw][k](i, j) = {re, im};
        }
      }
    }
  }

  const auto rw_str = !readQ ? "Written QPQ to " : "Read QPQ from ";
  std::cout << rw_str << "file: " << fname << "\n";
  return true;
}

//==============================================================================
void Feynman::form_qpiq(const std::string &ident) {
  if (read_qpiq(ident))
    return;
  form_qpiq(polarisation_wk());
  write_qpiq(ident);
}

//==============================================================================
std::vector<std::vector<ComplexRMatrix>> Feynman::polarisation_wk() const {
  std::cout << "Forming Pi(w,k)" << (m_hole_particle ? " (w/ hp)" : "")
            << " .. " << std::flush;

  const auto num_ws = m_wgrid_points.size();

  // BLAS must run single-threaded inside the omp region below
  const qip::SingleThreadBlas single_thread_blas{};

  // Parallel over w only: Green's fns are computed once per w, and re-used
  // across all k
  std::vector<std::vector<ComplexRMatrix>> pi_wk(num_ws);
#pragma omp parallel for schedule(dynamic)
  for (auto iw = 0ul; iw < num_ws; ++iw) {
    const auto omega = std::complex<double>{m_omre, m_wgrid_points[iw]};
    pi_wk[iw] = polarisation_each_k(omega, m_hole_particle);
  }

  std::cout << " done\n" << std::flush;
  return pi_wk;
}

//==============================================================================
void Feynman::form_qpiq(const std::vector<std::vector<ComplexRMatrix>> &pi_wk) {

  const auto num_ks = std::size_t(m_max_k + 1);
  const auto num_ws = m_wgrid_points.size();

  assert(pi_wk.size() == num_ws && "pi_wk must be on the same frequency grid");
  assert(pi_wk.front().size() == num_ks && "pi_wk must have k = 0..max_k");
  assert(pi_wk.front().front().size() == m_subgrid_points &&
         "pi_wk must be on the same sub-grid");

  m_qpiq_wk.resize(num_ws, num_ks,
                   ComplexRMatrix{m_i0, m_stride, m_subgrid_points, m_grid});

  // BLAS must run single-threaded inside the omp region below
  const qip::SingleThreadBlas single_thread_blas{};

  // q*pi*q products (+ screening): parallel over all (w, k); near instant
#pragma omp parallel for collapse(2)
  for (auto iw = 0ul; iw < num_ws; ++iw) {
    for (auto k = 0ul; k < num_ks; ++k) {
      const auto &q = get_qk(int(k)); // has drj
      const auto qdri = q.dri();      // has drj, and dri
      const auto &pi = pi_wk[iw][k];

      if (m_screen_Coulomb) {
        const auto X = X_screen(pi, qdri);
        m_qpiq_wk[iw][k] = q * pi * X * qdri;
      } else {
        m_qpiq_wk[iw][k] = q * pi * qdri;
      }
    }
  }
}

//==============================================================================
ComplexRMatrix Feynman::X_screen(const ComplexRMatrix &pik,
                                 const ComplexRMatrix &qdri) const {
  // X = [1 + i q*pi(w)]^-1
  constexpr auto Iunit = std::complex<double>{0.0, 1.0};
  return (Iunit * qdri * pik + 1.0).inverse();
}

//==============================================================================
GMatrix Feynman::Sigma_direct(int kv, double env,
                              std::optional<int> in_k) const {
  // If in_k is set, only calculate for single k
  // Used both for testing, and for calculating f_k factors

  auto Sigma_k = Sigma_direct_each_k(kv, env);

  if (in_k) {
    if (*in_k >= 0 && *in_k <= m_max_k) {
      return Sigma_k[std::size_t(*in_k)];
    }
    // k out of range: zero matrix
    return GMatrix{m_i0, m_stride, m_subgrid_points, m_include_G, m_grid};
  }

  GMatrix Sigma(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);
  for (const auto &Sk : Sigma_k) {
    Sigma += Sk;
  }
  return Sigma;
}

//==============================================================================
double L_exchange(int k, int l, int kappa_v, int kappa_alpha, int kappa_beta,
                  int kappa_gamma, const Angular::SixJTable &sixj) {
  const auto Ck_va = Angular::tildeCk_kk(k, kappa_v, kappa_alpha);
  const auto Ck_bg = Angular::tildeCk_kk(k, kappa_beta, kappa_gamma);
  const auto Cl_vg = Angular::tildeCk_kk(l, kappa_v, kappa_gamma);
  const auto Cl_ba = Angular::tildeCk_kk(l, kappa_beta, kappa_alpha);
  if (Ck_va == 0.0 || Ck_bg == 0.0 || Cl_vg == 0.0 || Cl_ba == 0.0)
    return 0.0;
  const auto tjv = Angular::twoj_k(kappa_v);
  const auto sj = sixj.get_2(tjv, Angular::twoj_k(kappa_alpha), 2 * k,
                             Angular::twoj_k(kappa_beta),
                             Angular::twoj_k(kappa_gamma), 2 * l);
  return Angular::neg1pow(k + l) * Ck_va * Ck_bg * Cl_vg * Cl_ba * sj /
         double(tjv + 1);
}

//==============================================================================
std::vector<std::vector<double>>
Feynman::subgrid_components(const DiracSpinor &Fa) const {
  const auto num_sp = m_include_G ? 2ul : 1ul;
  std::vector<std::vector<double>> F(num_sp,
                                     std::vector<double>(m_subgrid_points));
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    F[0][i] = Fa.f(m_i0 + i * m_stride);
    if (m_include_G) {
      F[1][i] = Fa.g(m_i0 + i * m_stride);
    }
  }
  return F;
}

//==============================================================================
GMatrix Feynman::Sigma_exchange(int kv, double env) const {
  // Exchange correlation potential, Methods Eq. (RadialSigmaExch):
  //   Sigma_12 = Int dw1/2pi Int dw2/2pi sum_{k l} sum_{alpha beta gamma}
  //     L^{kl}_{v beta alpha gamma}
  //     g^alpha_1i(e+w1) q^k_1j g^beta_ij(e+w1+w2) q^l_i2 g^gamma_j2(e+w2)
  // alpha, beta, gamma: partial waves of the lines 1i, ij, j2; k, l:
  // multipoles of the Coulomb lines 1j, i2; L includes (-1)^(k+l)/[j_v].
  // Spinor indices (f, g) run along the electron line 1i, ij, j2 (q is a
  // spinor scalar): Sigma^{mu nu}_12 = ... g^{mu s}_1i q_1j g^{st}_ij q_i2 g^{t nu}_j2
  //
  // The w1 integral is done analytically, by closing the contour on the
  // core poles [Methods Eq. (ExchInt1)]:
  //   Int dw1/2pi g^alpha_1i(e+w1) g^beta_ij(e+w1+w) = i Gamma^{alpha beta}_1iij(w)
  //   Gamma^{alpha beta}_1iij(w) = sum_{a in alpha} p^a_1i gex^beta_ij(e_a+w)
  //                              + sum_{a in beta} gex^alpha_1i(e_a-w) p^a_ij
  // p^a_ij = F_a(r_i) F_a(r_j)^T is the core projector [p^{st}_ij =
  // F_a^s(r_i) F_a^t(r_j)], gex is the Green's function of the excited
  // states only, 'a in alpha' runs over the core states of partial wave
  // alpha.
  //
  // The w2 = w = omre + iu integral is done numerically:
  //   Sigma_12 = -(1/pi) Re Int_0^inf du sum_{k l} sum_{alpha beta gamma}
  //     L^{kl}_{v beta alpha gamma} Gamma^{alpha beta}_1iij(w) q^k_1j q^l_i2 g^gamma_j2(e+w)
  // -1/pi = (1/2pi) x i (dw = i du) x i (from Gamma) x 2 (the integrand at
  // -u is the conjugate of that at +u). The integrand is finite at u = 0
  // and falls as u^-2.
  //
  // Gamma is not formed. Since p^a is rank one, each of its two terms, with
  // q^l_i2 and L attached and summed over the internal partial wave and l
  // (exchange_Gamma_q), is a product of matrices:
  //   sum_{beta l} L^{kl}_{v beta a gamma} [p^a_1i gex^beta_ij q^l_i2]^{mu t} = F_a^mu(r_1) [pa_gex]^t_j2
  //   sum_{alpha l} L^{kl}_{v a alpha gamma} [gex^alpha_1i p^a_ij q^l_i2]^{mu t} = [gex_pa]^mu_12 F_a^t(r_j)
  // which leaves, for each core state a, and each (k, gamma):
  //   Sigma^{mu nu}_12 += F_a^mu(r_1) sum_t [q^k (g^{t nu} o [pa_gex]^t)]_12
  //                     + [gex_pa]^mu_12 [q^k sum_t F_a^t g^{t nu}]_12
  // ('o': element-wise product). Integration measures: q^k_1j carries dr_j,
  // q^l_i2 carries dr_i dr_2 (so Sigma carries dr_2, as for the direct term).

  const auto num_ks = std::size_t(m_max_k + 1);
  const auto num_kappas = m_max_ki + 1;
  const auto num_sp = m_include_G ? 2ul : 1ul;
  const Angular::SixJTable sixj(2 * m_max_k);

  // Quadrature weights: the u^-3 tail correction of the w grid (direct term)
  // becomes u_max F(u_max) for the u^-2 exchange integrand
  auto weights = m_wgrid_weights;
  weights.back() += 0.5 * m_wgrid_points.back();

  // BLAS must run single-threaded inside the omp region below
  const qip::SingleThreadBlas single_thread_blas{};

  const ComplexRMatrix zero(m_i0, m_stride, m_subgrid_points, m_grid);
  const GMatrix Sigma_zero(m_i0, m_stride, m_subgrid_points, m_include_G,
                           m_grid);
  std::vector<GMatrix> Sigma_ts(std::size_t(omp_get_max_threads()), Sigma_zero);

#pragma omp parallel for schedule(dynamic)
  for (auto iw = 0ul; iw < m_wgrid_points.size(); ++iw) {
    auto &Sigma_t = Sigma_ts[std::size_t(omp_get_thread_num())];

    const auto w = std::complex<double>{m_omre, m_wgrid_points[iw]};
    const auto du = weights[iw];

    // Valence-line Green's function, g^gamma(e + w), for each partial wave
    std::vector<ComplexGMatrix> g_gamma;
    for (auto ig = 0ul; ig < num_kappas; ++ig) {
      g_gamma.push_back(green(Angular::kindex_to_kappa(ig), env + w));
    }

    for (const auto &Fa : m_HF->core()) {
      if (Fa.n() < m_min_core_n)
        continue;
      const auto ka = Fa.kappa();
      const auto F = subgrid_components(Fa);
      const auto Gamma = exchange_Gamma_q(kv, Fa, w, sixj);

      for (auto k = 0ul; k < num_ks; ++k) {
        const auto &qk = get_qk(int(k));
        // pa_gex requires C^k_{v a}; gex_pa requires C^k_{a gamma}
        const bool pa_gex_nonzero = Angular::Ck_kk_SR(int(k), kv, ka);

        for (auto ig = 0ul; ig < num_kappas; ++ig) {
          const auto kg = Angular::kindex_to_kappa(ig);
          const auto &g = g_gamma[ig];
          const bool gex_pa_nonzero = Angular::Ck_kk_SR(int(k), ka, kg);

          for (auto nu = 0ul; nu < num_sp; ++nu) {

            if (pa_gex_nonzero) {
              // F_a^mu(r_1) sum_t [q^k (g^{t nu} o [pa_gex]^t)]_12
              ComplexRMatrix qk_g_Gamma(zero);
              for (auto t = 0ul; t < num_sp; ++t) {
                qk_g_Gamma +=
                  qk * mult_elements(g.radial(t, nu), Gamma.pa_gex[t](k, ig));
              }
              for (auto mu = 0ul; mu < num_sp; ++mu) {
                Sigma_t.sp(mu, nu) +=
                  du * mult_rows(qk_g_Gamma, F[mu]).real().Rmatrix();
              }
            }

            if (gex_pa_nonzero) {
              // [gex_pa]^mu_12 [q^k sum_t F_a^t g^{t nu}]_12
              ComplexRMatrix Fa_g(zero);
              for (auto t = 0ul; t < num_sp; ++t) {
                Fa_g += mult_rows(g.radial(t, nu), F[t]);
              }
              const auto qk_Fa_g = qk * Fa_g;
              for (auto mu = 0ul; mu < num_sp; ++mu) {
                Sigma_t.sp(mu, nu) +=
                  du * mult_elements(Gamma.gex_pa[mu](k, ig), qk_Fa_g)
                         .real()
                         .Rmatrix();
              }
            }
          }
        }
      }
    }
  }

  GMatrix Sigma(Sigma_zero);
  for (const auto &Sigma_t : Sigma_ts) {
    Sigma += Sigma_t;
  }
  Sigma *= (-1.0 / M_PI);
  return Sigma;
}

//==============================================================================
Feynman::GammaQ
Feynman::exchange_Gamma_q(int kv, const DiracSpinor &Fa, std::complex<double> w,
                          const Angular::SixJTable &sixj) const {
  // The two terms of Gamma (see Sigma_exchange) for core state a, with the
  // Coulomb line q^l_i2 and the angular factor L attached, summed over the
  // internal partial wave (beta, alpha) and l. With p^{st}_ij = F_a^s(r_i)
  // F_a^t(r_j), the factor of p^a outside the i sum comes out:
  //   sum_{beta l} L^{kl}_{v beta a gamma} [p^a_1i gex^beta_ij(e_a+w) q^l_i2]^{mu t} = F_a^mu(r_1) [pa_gex]^t_j2
  //   sum_{alpha l} L^{kl}_{v a alpha gamma} [gex^alpha_1i(e_a-w) p^a_ij q^l_i2]^{mu t} = [gex_pa]^mu_12 F_a^t(r_j)
  // so that
  //   [pa_gex]^t = sum_{beta l} L^{kl}_{v beta a gamma} [gex^beta(e_a+w) F_a q^l]^t
  //   [gex_pa]^mu = sum_{alpha l} L^{kl}_{v a alpha gamma} [gex^alpha(e_a-w) F_a q^l]^mu
  // with [gex F_a]^t_ji = sum_s gex^{ts}_ji F_a^s(r_i) (using gex_ij = gex_ji^T).
  // Returned for each (k, gamma), and each spinor index t (mu)

  const auto num_ks = std::size_t(m_max_k + 1);
  const auto num_kappas = m_max_ki + 1;
  const auto num_sp = m_include_G ? 2ul : 1ul;

  const auto ka = Fa.kappa();
  const auto ea = std::complex<double>{Fa.en()};
  const auto F = subgrid_components(Fa);

  const ComplexRMatrix zero(m_i0, m_stride, m_subgrid_points, m_grid);
  const LinAlg::Matrix<ComplexRMatrix> zeros(num_ks, num_kappas, zero);
  GammaQ Gamma{std::vector(num_sp, zeros), std::vector(num_sp, zeros)};

  for (auto ik = 0ul; ik < num_kappas; ++ik) {
    const auto kappa = Angular::kindex_to_kappa(ik);

    const auto gex_plus = green_excited(kappa, ea + w);
    const auto gex_minus = green_excited(kappa, ea - w);

    // [gex F_a]^t = sum_s gex^{ts} F_a^s
    std::vector<ComplexRMatrix> gexFa_plus(num_sp, zero),
      gexFa_minus(num_sp, zero);
    for (auto t = 0ul; t < num_sp; ++t) {
      for (auto s = 0ul; s < num_sp; ++s) {
        gexFa_plus[t] += mult_cols(gex_plus.radial(t, s), F[s]);
        gexFa_minus[t] += mult_cols(gex_minus.radial(t, s), F[s]);
      }
    }

    for (auto l = 0; l <= m_max_k; ++l) {
      // q^l joins core state a to the internal line
      if (!Angular::Ck_kk_SR(l, ka, kappa))
        continue;
      const auto ql = get_qk(l).dri();

      // [gex(e_a +/- w) F_a q^l]^t
      std::vector<ComplexRMatrix> gexFaq_plus, gexFaq_minus;
      for (auto t = 0ul; t < num_sp; ++t) {
        gexFaq_plus.push_back(gexFa_plus[t] * ql);
        gexFaq_minus.push_back(gexFa_minus[t] * ql);
      }

      for (auto k = 0ul; k < num_ks; ++k) {
        for (auto ig = 0ul; ig < num_kappas; ++ig) {
          const auto kg = Angular::kindex_to_kappa(ig);
          const auto L_pa_gex = L_exchange(int(k), l, kv, ka, kappa, kg, sixj);
          const auto L_gex_pa = L_exchange(int(k), l, kv, kappa, ka, kg, sixj);
          for (auto t = 0ul; t < num_sp; ++t) {
            if (L_pa_gex != 0.0) {
              Gamma.pa_gex[t](k, ig) += L_pa_gex * gexFaq_plus[t];
            }
            if (L_gex_pa != 0.0) {
              Gamma.gex_pa[t](k, ig) += L_gex_pa * gexFaq_minus[t];
            }
          }
        }
      }
    }
  }
  return Gamma;
}

//==============================================================================
std::vector<GMatrix> Feynman::Sigma_direct_each_k(int kv, double env) const {
  // Direct Sigma, for each multipole k separately: Sigma_d = sum_k Sigma_d^k.
  // Green's functions are shared by all k, so calculating every k at once
  // costs the same as a single k (used for the effective screening factors
  // fk, which need the ratio of each k term separately).

  assert(has_qpiq() && "Q*Pi*Q must be formed (form_qpiq) before Sigma_direct");

  const auto num_ks = std::size_t(m_max_k + 1);
  const GMatrix zero(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);
  std::vector<GMatrix> Sigma_k(num_ks, zero);

  constexpr std::complex<double> I{0.0, 1.0};
  const auto num_kappas = m_max_ki + 1;

  // BLAS must run single-threaded inside the omp region below
  const qip::SingleThreadBlas single_thread_blas{};

  std::vector<std::vector<GMatrix>> Sigma_ts(std::size_t(omp_get_max_threads()),
                                             Sigma_k);

#pragma omp parallel for collapse(2)
  for (auto iw = 0ul; iw < m_wgrid_points.size(); iw++) {
    for (auto iB = 0ul; iB < num_kappas; ++iB) {

      auto &Sigma_t = Sigma_ts[std::size_t(omp_get_thread_num())];

      const auto omega = std::complex{m_omre, m_wgrid_points[iw]};

      // I, since dw is on imag. grid
      const auto dw = I * m_wgrid_weights[iw];

      const auto kB = Angular::kindex_to_kappa(iB);

      const auto gB = green(kB, env + omega);

      for (auto k = 0ul; k < num_ks; k++) {

        const auto ck_vB = Angular::Ck_kk(int(k), kv, kB);
        if (ck_vB == 0.0)
          continue;

        const auto &qpq_dw = m_qpiq_wk[iw][k];

        const auto c_ang_dw =
          dw * ck_vB * ck_vB / double(Angular::twoj_k(kv) + 1);

        Sigma_t[k] += (c_ang_dw * mult_elements(gB, qpq_dw)).real();
      }
    }
  }

  for (const auto &Sigma_t : Sigma_ts) {
    for (auto k = 0ul; k < num_ks; k++) {
      Sigma_k[k] += Sigma_t[k];
    }
  }

  // 1/pi = (1/2pi)*2; the 2 from the symmetric +/- w contour halves
  // (integrand at -w is conjugate of that at +w, and we take Re part)
  for (auto &Sk : Sigma_k) {
    Sk *= (1.0 / M_PI);
  }

  return Sigma_k;
}

} // namespace MBPT