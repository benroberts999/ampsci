#include "Feynman.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracODE/DiracODE.hpp"
#include "HF/HartreeFock.hpp"
#include "MBPT/RDMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/color.hpp"
#include <cassert>
#include <memory>
#include <vector>

namespace MBPT {

//==============================================================================
Feynman::Feynman(const HF::HartreeFock *vHF, MBPT::rgrid_params subgrid,
                 double omre, double w_0, double w_ratio, Screening scr_option,
                 HoleParticle hp_option, int max_l, int n_min_core,
                 bool verbose)
    : m_HF(vHF),
      m_grid(vHF->grid_sptr()),
      m_i0(vHF->grid().getIndex(subgrid.r0)),
      m_imax(vHF->grid().getIndex(subgrid.rmax)),
      m_stride(subgrid.stride), // check lines up? XXX
      m_subgrid_points((m_imax - m_i0) / m_stride + 1),
      m_max_ki_core(2 * DiracSpinor::max_l(m_HF->core())),
      m_max_ki(std::max(2 * max_l, m_max_ki_core)),
      m_max_k(std::max(m_max_ki_core + 1, m_max_ki + 1)),
      m_min_core_n(n_min_core),
      m_omre(omre),
      m_wgrid(form_w_grid(w_0, w_ratio)),
      m_verbose(verbose),
      m_hole_particle(hp_option == HoleParticle::include ||
                      hp_option == HoleParticle::include_k0),
      m_include_higher_order_hp(hp_option != HoleParticle::include_k0),
      m_screen_Coulomb(scr_option == Screening::include),
      m_only_screen(scr_option == Screening::only),
      // creats empty arrays, these are filled below
      m_dri(m_i0, m_stride, m_subgrid_points, false, m_grid),
      m_drj(m_i0, m_stride, m_subgrid_points, false, m_grid) {

  if (verbose) {
    std::cout << "\nFeynman diagrams for correlation potential\n";
    std::cout << "lmax = " << Angular::lFromIndex(m_max_ki)
              << " (internal lines)\n";
    std::cout << "Including from core n=" << m_min_core_n << "\n";
    if (m_screen_Coulomb) {
      std::cout << "Including all-orders Coulomb screening\n";
    }
    if (m_only_screen) {
      std::cout << "Only including high-orders Coulomb screening (2nd "
                   "order not included!)\n";
    }
    if (m_hole_particle) {
      std::cout << "Including hole-particle interaction: "
                << (m_include_higher_order_hp ? "(all k in hp)\n" :
                                                "(only k=0)\n");
    }
    std::cout << "Re(w) = " << m_omre << "\n";
    std::cout << "Im(w) : " << m_wgrid.gridParameters();
    printf(", r=%.2f\n", m_wgrid.r(1) / m_wgrid.r(0));

    printf("Sigma sub-grid: r=(%.1e, %.1f)aB with %i points. [i0=%i, "
           "stride=%i]\n",
           m_grid->r(m_i0), m_grid->r(m_i0 + m_stride * (m_subgrid_points - 1)),
           int(m_subgrid_points), int(m_i0), int(m_stride));
  }

  if (m_HF->vBreit()) {
    // Breit was included into HF; we should not do this!
    std::cout
        << "\nWARNING: Trying to calculate Sigma in Feynman method when Breit "
           "is included will result in incorrect results! Calculate Sigma "
           "first (without Breit), then you may include Breit\n";
  }

  // work-around for polarisation issue
  check_min_n();
  // Construct qk, and dri/drj
  form_qk();
  // Construct pa, Vx
  form_pa();
  form_vx();
  // Construct polarisation operator
  form_qpiq();
}

//==============================================================================
void Feynman::check_min_n() {

  // This is a hacky workaround for the polarisation problem

  double cut_off = 137.0; // Is this a coincidence???

  double max_en_core = 0;
  for (const auto &n : m_HF->core()) {
    if (n.n() < m_min_core_n)
      continue;
    if (-n.en() > max_en_core) {
      max_en_core = -n.en();
    }
  }
  if (max_en_core > cut_off) {
    fmt2::warning();
    std::cout << "Core shell with n = " << m_min_core_n
              << " too deep with |En| = " << max_en_core
              << " -- updating n_min_core to " << m_min_core_n + 1 << "\n";
    ++m_min_core_n;
    // do recursively until OK (could be more eficient,
    // but I like the multiple warning messages when n_min is too small)
    check_min_n();
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

  m_qk = std::vector<ComplexGMatrix>(
      size, {m_i0, m_stride, m_subgrid_points, false, m_HF->grid_sptr()});

  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto fi = m_qk.front().index_to_fullgrid(i);
    // const auto dri = m_grid->drdu(fi) * m_grid->du() * double(m_stride);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto fj = m_qk.front().index_to_fullgrid(j);
      const auto drj = m_grid->drdu(fj) * m_grid->du() * double(m_stride);
      const auto rless = m_grid->r(std::min(fi, fj));
      const auto rgreater = m_grid->r(std::max(fi, fj));
      const auto ratio = rless / rgreater; // = r_< / r_>
      m_qk[0].ff(i, j) = (1.0 / rgreater) * drj;
      m_drj.ff(i, j) = drj;
      for (int k = 1; k <= m_max_k; ++k) {
        const auto sk = std::size_t(k);
        m_qk[sk].ff(i, j) = m_qk[sk - 1].ff(i, j) * ratio;
      }
    }
  }
  m_dri.ff() = m_drj.ff().transpose();
}

//==============================================================================
void Feynman::form_pa() {
  // Fill core |a><a|
  const auto &core = m_HF->core();

  m_pa.resize(core.size(), {m_i0, m_stride, m_subgrid_points, false, m_grid});

  for (auto ia = 0ul; ia < core.size(); ia++) {
    m_pa[ia] = green_single(core[ia], core[ia], std::complex<double>{1.0, 0.0});
  }
}

//==============================================================================
ComplexGMatrix Feynman::green_single(const DiracSpinor &ket,
                                     const DiracSpinor &bra,
                                     const std::complex<double> f) const {

  ComplexGMatrix Gmat(m_i0, m_stride, m_subgrid_points, false, m_grid);
  Gmat.add(ket, bra, f);
  return Gmat;
}

//==============================================================================
void Feynman::form_vx() {
  // Vx = -|a>Q<a|

  m_Vx_kappa =
      std::vector<GMatrix>(std::size_t(m_max_ki + 1),
                           {m_i0, m_stride, m_subgrid_points, false, m_grid});

  for (int kapi = 0; kapi <= m_max_ki; ++kapi) {

    const auto kappa = Angular::kappaFromIndex(kapi);
    const auto twojp1 = Angular::twoj_k(kappa) + 1;

    for (int k = 0; k <= m_max_k; ++k) {
      // GMatrix V_k(m_i0, m_stride, m_subgrid_points, false, m_grid);
      ComplexGMatrix V_k(m_i0, m_stride, m_subgrid_points, false, m_grid);

      const auto &core = m_HF->core();
      for (auto ia = 0ul; ia < core.size(); ia++) {

        const auto ck = Angular::Ck_kk(k, kappa, core.at(ia).kappa());
        if (ck == 0.0)
          continue;
        const auto c_ang = -1.0 * ck * ck / double(twojp1);
        // V_k.add(core[ia], core[ia], c_ang);
        // V_k += c_ang * m_pa[ia].real();
        V_k += c_ang * m_pa[ia];
      }

      const auto sk = std::size_t(k);

      m_Vx_kappa[std::size_t(kapi)] +=
          V_k.mult_elements_by(m_qk.at(sk).dri()).real();
    }
  }
}

//==============================================================================
ComplexGMatrix Feynman::green(int kappa, std::complex<double> en,
                              GreenStates states) const {
  if (states == GreenStates::core) {
    return green_core(kappa, en);
  } else if (states == GreenStates::excited) {
    return green_excited(kappa, en);
  }
  return m_Complex_green_method ? green_hf_v2(kappa, en) : green_hf(kappa, en);
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

  // Evaluate Wronskian at ~65% of the way to pinf. Should be inependent of r
  const auto pp = std::size_t(0.65 * double(xI.max_pt()));
  const auto w = (x0.f(pp) * xI.g(pp) - xI.f(pp) * x0.g(pp)) / alpha;

  // Get G0 (Green's function, without exchange):
  const auto g0 = construct_green_g0(x0, xI, w);

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
    return -1.0 * ((g0 * Vx - 1.0).invert_in_place() * g0).complex();
  }

  // G0 := G0(re{e}) - no exchange, only real part
  // G(e) = [1 + i*Im{e}*G0 - G0*Vx]^{-1} * G0
  // Note: differential dr is included in Vx (via Q)
  std::complex<double> iw{0.0, en.imag()};
  return (iw * g0.complex().drj_in_place() - (g0 * Vx).complex() + 1.0)
             .invert_in_place() *
         g0.complex();
}

//==============================================================================
ComplexGMatrix Feynman::green_hf_v2(int kappa, std::complex<double> en,
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

  // Evaluate Wronskian at ~65% of the way to pinf. Should be inependent of r
  const auto pp = std::size_t(0.65 * double(xI.max_pt()));
  const auto I = std::complex{0.0, 1.0};
  // f -> f_r + I*f_i and same for g..
  const auto w = ((x0.f(pp) + I * Ix0.f(pp)) * (xI.g(pp) + I * IxI.g(pp)) -
                  (xI.f(pp) + I * IxI.f(pp)) * (x0.g(pp) + I * Ix0.g(pp))) /
                 alpha;

  // Get G0 (Green's function, without exchange):
  const auto g0 = construct_green_g0(x0, Ix0, xI, IxI, w);

  // Include exchange (optionally, with hole-particle correction)
  auto Vx = get_Vx_kappa(kappa);
  if (Fc_hp != nullptr) {
    // Include hole-particle interaction (w/ [1-P]V[1-P]):
    Vx += calculate_Vhp(kappa, *Fc_hp);
  }

  // Include exchange using Dyson:
  return -1.0 * ((g0 * Vx.complex() - 1.0).invert_in_place() * g0);
}

//==============================================================================
ComplexGMatrix Feynman::green_core(int kappa, std::complex<double> en) const {
  // G_core = \sum_a |a><a|/(e_r + i*e_i-ea), for all a with a.kappa() = k
  ComplexGMatrix Gcore(m_i0, m_stride, m_subgrid_points, false, m_grid);

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
  ComplexGMatrix Gk(m_i0, m_stride, m_subgrid_points, false, m_grid);

  // Subtract core states, by forceing Gk to be orthogonal to core:
  // Gk -> Gk - \sum_a|a><a|G

  const auto g0 = m_Complex_green_method ? green_hf_v2(kappa, en, Fc_hp) :
                                           green_hf(kappa, en, Fc_hp);
  // nb: can also subtract of core part, but doesn't seem to make a difference
  // return orthogonalise_wrt_core(g0 - green_core(kappa, en), kappa);
  return orthogonalise_wrt_core(g0, kappa);
}

//==============================================================================
ComplexGMatrix Feynman::orthogonalise_wrt_core(const ComplexGMatrix &g_in,
                                               int kappa) const {
  // Force Gk to be orthogonal to the core states
  const auto &core = m_HF->core();
  // const auto &drj = get_drj();
  auto g_out = g_in;
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    if (core[ia].kappa() == kappa) {
      g_out -= m_pa[ia].drj() * g_in;
    }
  }
  return g_out;
}

//==============================================================================
GMatrix Feynman::calculate_Vhp(int kappa, const DiracSpinor &Fc) const {
  // Hole-particle interaction, extra potential
  // This works well for k=0, but works worse than using local pot for k=1
  // (k=1 is most important term!)

  GMatrix V0(m_i0, m_stride, m_subgrid_points, false, m_grid);

  const auto y0cc = Coulomb::yk_ab(0, Fc, Fc);
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = V0.index_to_fullgrid(i);
    V0.ff(i, i) = -y0cc[si];
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
      }
    }
  }

  V0.drj_in_place();

  GMatrix OneNegPc(m_i0, m_stride, m_subgrid_points, false, m_HF->grid_sptr());
  const auto &core = m_HF->core();
  for (std::size_t ia = 0; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    if (a.kappa() != kappa)
      continue;
    OneNegPc.add(a, a, -1.0);
    // OneNegPc -= m_Pa[ia].real();
  }
  OneNegPc.drj_in_place();
  OneNegPc += 1.0; // (1-P)
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
  GMatrix g0I(m_i0, m_stride, m_subgrid_points, false, m_grid);

  const auto winv = 1.0 / w;

  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = g0I.index_to_fullgrid(i);
    for (auto j = 0ul; j <= i; ++j) { // j <= i
      const auto sj = g0I.index_to_fullgrid(j);
      g0I.ff(i, j) = x0.f(sj) * xI.f(si) * winv;
      // g0I is symmetric
      g0I.ff(j, i) = g0I.ff(i, j);
    } // j
  }   // i

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
  ComplexGMatrix g0I(m_i0, m_stride, m_subgrid_points, false, m_grid);

  const auto winv = 1.0 / w;

  const auto I = std::complex{0.0, 1.0};

  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = g0I.index_to_fullgrid(i);
    for (auto j = 0ul; j <= i; ++j) { // j <= i
      const auto sj = g0I.index_to_fullgrid(j);
      g0I.ff(i, j) =
          (x0.f(sj) + I * Ix0.f(sj)) * (xI.f(si) + I * IxI.f(si)) * winv;
      // g0I is symmetric
      g0I.ff(j, i) = g0I.ff(i, j);
    } // j
  }   // i

  return g0I;
}

//==============================================================================
ComplexGMatrix Feynman::polarisation_k(int k,
                                       std::complex<double> omega) const {

  ComplexGMatrix pi_k(m_i0, m_stride, m_subgrid_points, false, m_grid);

  const auto Iunit = std::complex<double>{0.0, 1.0};
  const auto &core = m_HF->core();
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    if (a.n() < m_min_core_n)
      continue;

    const auto &pa = m_pa[ia]; // |a><a|
    const auto ea_minus_w = std::complex<double>{a.en()} - omega;
    const auto ea_plus_w = std::complex<double>{a.en()} + omega;

    const auto *Fa_hp = m_hole_particle ? &a : nullptr;

    for (int in = 0; in <= m_max_ki; ++in) {
      const auto kn = Angular::kappaFromIndex(in);
      const auto ck_an = Angular::Ck_kk(k, a.kappa(), kn);
      if (ck_an == 0.0)
        continue;
      const double c_ang = ck_an * ck_an / double(2 * k + 1);

      // XXX Can make this faster:
      // 1. g at e_real
      // 2. "rotate" to each imag.
      // 3. Orthogonalise entire thing wrt core
      // 4. Sum up each kappa_n, _then_ multiply by p_a

      pi_k += c_ang * ((green_excited(kn, ea_minus_w, Fa_hp) +
                        green_excited(kn, ea_plus_w, Fa_hp))
                           .mult_elements_by(pa));
    }
  }
  pi_k *= Iunit;
  return pi_k;
}

//==============================================================================
Grid Feynman::form_w_grid(double w0, double wratio) const {

  // Find max core energy: (for w_max)
  auto wmax_core = 30.0; // don't let it go below 50
  const auto &core = m_HF->core();
  for (const auto &Fc : core) {
    if (Fc.n() < m_min_core_n)
      continue;
    if (std::abs(Fc.en()) > wmax_core)
      wmax_core = std::abs(Fc.en());
  }

  // maximum Im(w): based on core energy.
  const auto wmax = 2.0 * wratio * wmax_core;

  const std::size_t wsteps = Grid::calc_num_points_from_du(
      w0, wmax, std::log(wratio), GridType::logarithmic);

  auto wgrid = Grid(w0, wmax, wsteps, GridType::logarithmic);

  return wgrid;
}

//==============================================================================
void Feynman::form_qpiq() {

  if (m_verbose) {
    std::cout << "Forming QPQ(w,k)";
    if (m_hole_particle || m_screen_Coulomb) {
      std::cout << " (w/ " << (m_screen_Coulomb ? "scr" : "")
                << (m_hole_particle && m_screen_Coulomb ? " + " : "")
                << (m_hole_particle ? "hp" : "") << ")";
    }
    std::cout << " .. " << std::flush;
  }

  const auto num_ks = std::size_t(m_max_k + 1);

  m_qpiq_wk.resize(
      m_wgrid.num_points(),
      std::vector<ComplexGMatrix>(
          num_ks, {m_i0, m_stride, m_subgrid_points, false, m_grid}));

#pragma omp parallel for collapse(2)
  for (auto iw = 0ul; iw < m_wgrid.num_points(); ++iw) {
    for (auto k = 0ul; k < num_ks; ++k) {
      const auto omega = std::complex<double>{m_omre, m_wgrid.r(iw)};

      const auto &q = m_qk.at(k); // has drj
      const auto qdri = q.dri();  // has drj, and dri
      const auto pi = polarisation_k(int(k), omega);

      if (m_screen_Coulomb && !m_only_screen) {
        const auto X = X_screen(pi, qdri);
        m_qpiq_wk[iw][k] = q * pi * X * qdri;
      } else if (m_only_screen) {
        const auto X = X_screen(pi, qdri);
        m_qpiq_wk[iw][k] = q * pi * (X - 1.0) * qdri;
      } else {
        m_qpiq_wk[iw][k] = q * pi * qdri;
      }
    }
  }

  if (m_verbose) {
    std::cout << " done\n" << std::flush;
  }
}

//==============================================================================
ComplexGMatrix Feynman::X_screen(const ComplexGMatrix &pik,
                                 const ComplexGMatrix &qdri) const {
  // X = [1 + i q*pi(w)]^-1
  const auto Iunit = std::complex<double>{0.0, 1.0};
  return (Iunit * qdri * pik + 1).inverse();
}

//==============================================================================
GMatrix Feynman::Sigma_direct(int kv, double env,
                              std::optional<int> in_k) const {

  // If in_k is set, only calculate for single k
  // Used both for testing, and for calculating f_k factors

  GMatrix Sigma(m_i0, m_stride, m_subgrid_points, false, m_grid);

  const std::complex<double> I{0.0, 1.0};

  // Store gAs in advance
  const auto num_kappas = std::size_t(m_max_ki + 1);

#pragma omp parallel for collapse(2)
  for (auto iw = 0ul; iw < m_wgrid.num_points(); iw++) { // for omega integral
    for (auto iB = 0ul; iB < num_kappas; ++iB) {

      const auto omim = m_wgrid(iw);

      // Simpson's rule: Implicit ends (integrand zero at w=0 and w>wmax)
      const auto weight = iw % 2 == 0 ? 4.0 / 3 : 2.0 / 3;

      // I, since dw is on imag. grid; 2 from symmetric +/- w
      const auto dw = I * weight * m_wgrid.drdu(iw);

      const auto kB = Angular::kappaFromIndex(int(iB));

      // const auto &gB = gBs[iB][iw]; can be more efficient
      const auto gB =
          green(kB, env + std::complex{m_omre, omim}, GreenStates::both);

      for (auto k = 0ul; int(k) <= m_max_k; k++) {

        // For doing single k
        if (in_k && in_k != int(k))
          continue;

        const auto ck_vB = Angular::Ck_kk(int(k), kv, kB);
        if (ck_vB == 0.0)
          continue;

        const auto &qpq_dw = m_qpiq_wk[iw][k];

        const auto c_ang_dw =
            dw * ck_vB * ck_vB / double(Angular::twoj_k(kv) + 1);
        const auto C_gB_QPQ_dw = (c_ang_dw * mult_elements(gB, qpq_dw)).real();

#pragma omp critical(sum_sigma_d)
        { Sigma += C_gB_QPQ_dw; }
      }
    }
  }

  // Extra 2 from symmetric + / -w
  // 2*du/(2pi) = du / pi
  Sigma *= (m_wgrid.du() / M_PI);

  return Sigma;
} // namespace MBPT

} // namespace MBPT