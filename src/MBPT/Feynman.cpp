#include "Feynman.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracODE/include.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "MBPT/RadialMatrix.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include <cassert>
#include <memory>
#include <vector>

namespace MBPT {

//==============================================================================
Feynman::Feynman(const HF::HartreeFock *vHF, std::size_t i0, std::size_t stride,
                 std::size_t size, const FeynmanOptions &options,
                 int n_min_core, bool include_G, bool verbose)
    : m_HF(vHF),
      m_grid(vHF->grid_sptr()),
      m_i0(i0),
      m_stride(stride),
      m_subgrid_points(size),
      m_max_ki_core(2 * DiracSpinor::max_l(m_HF->core())),
      m_max_ki(std::max(2 * options.max_l_internal, m_max_ki_core)),
      m_max_k(std::max(m_max_ki_core + 1, m_max_ki + 1)),
      m_min_core_n(n_min_core),
      m_include_G(include_G),
      m_omre(options.omre),
      m_wgrid(form_w_grid(options.w0, options.w_ratio)),
      m_hole_particle(options.hole_particle == HoleParticle::include ||
                      options.hole_particle == HoleParticle::include_k0),
      m_include_higher_order_hp(options.hole_particle !=
                                HoleParticle::include_k0),
      m_screen_Coulomb(options.screening == Screening::include) {

  if (verbose) {
    std::cout << "\nFeynman diagrams:\n";
    fmt::print("lmax = {} (internal lines)\n", Angular::lFromIndex(m_max_ki));
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
    std::cout << "Im(w) : " << m_wgrid.gridParameters();
    printf(", r=%.2f\n", m_wgrid.r(1) / m_wgrid.r(0));
  }

  // Construct qk, and dri/drj
  form_qk();
  // Construct pa, Vx
  form_pa();
  form_vx();
  // Construct polarisation operator
  form_qpiq();
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

  m_pa.resize(core.size(), {m_i0, m_stride, m_subgrid_points, true, m_grid});

  for (auto ia = 0ul; ia < core.size(); ia++) {
    m_pa[ia] = green_single(core[ia], core[ia], std::complex<double>{1.0, 0.0});
  }
}

//==============================================================================
ComplexGMatrix Feynman::green_single(const DiracSpinor &ket,
                                     const DiracSpinor &bra,
                                     const std::complex<double> f) const {
  ComplexGMatrix Gmat(m_i0, m_stride, m_subgrid_points, true, m_grid);
  Gmat.add(ket, bra, f);
  return Gmat;
}

//==============================================================================
// new test way to form the radial exchange potential coordinate matrix using change of basis formula -- the complete basis that is used is the hydrogenic wave functions
void Feynman::form_vx() {

  m_Vx_kappa =
      std::vector<GMatrix>(std::size_t(m_max_ki + 1),
                           {m_i0, m_stride, m_subgrid_points, true, m_grid});

  // Use basic H-like basis to express matrices in coordinate/orbital space
  // initialises the hydrogen wave function object
  Wavefunction wfH(m_grid, {"1", 0, "Ball"});
  wfH.set_HF();
  wfH.solve_core(false);

  const auto r0 = std::max(1.0e-4, 3.0 * m_grid->r0());
  const auto rmax = std::min(90.0, 0.9 * m_grid->rmax());
  const int max_n = 90;
  const auto max_l =
      std::max(Angular::lFromIndex(m_max_ki), DiracSpinor::max_l(m_HF->core()));
  const auto l_string =
      AtomData::spectroscopic_notation.substr(0, std::size_t(max_l));
  const auto basis_string = std::to_string(max_n) + l_string;

  wfH.formBasis(SplineBasis::Parameters(
      basis_string, 90, 9, r0, 0.0, rmax, basis_string,
      SplineBasis::SplineType::Derevianko, false, false));

  // constructs the exchange matrix as:
  // V(r1,r2) = \sum_n [Vx*F_n](r1) F_n^†(r2)
  for (const auto &Fn : wfH.basis()) {

    auto VxFn = m_HF->vexFa(Fn);
    if (m_HF->vBreit()) {
      VxFn += m_HF->VBr(Fn);
    }

    m_Vx_kappa[std::size_t(Fn.k_index())].add(VxFn, Fn, 1.0);
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

//     const auto kappa = Angular::kappaFromIndex(kapi);
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

  // Evaluate Wronskian at ~65% of the way to pinf. Should be independent of r
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

  // Evaluate Wronskian at ~65% of the way to pinf. Should be independent of r
  const auto pp = std::size_t(0.65 * double(xI.max_pt()));
  const auto I = std::complex{0.0, 1.0};
  // f -> f_r + I*f_i and same for g..
  const auto w = ((x0.f(pp) + I * Ix0.f(pp)) * (xI.g(pp) + I * IxI.g(pp)) -
                  (xI.f(pp) + I * IxI.f(pp)) * (x0.g(pp) + I * Ix0.g(pp))) /
                 alpha;
  // Should include conjugate??

  // Get G0 (Green's function, without exchange):
  // XXX Must be checked
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
  ComplexGMatrix Gcore(m_i0, m_stride, m_subgrid_points, true, m_grid);

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
  ComplexGMatrix Gk(m_i0, m_stride, m_subgrid_points, true, m_grid);

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

  GMatrix V0(m_i0, m_stride, m_subgrid_points, true, m_grid);

  const auto y0cc = Coulomb::yk_ab(0, Fc, Fc);
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = V0.index_to_fullgrid(i);
    V0.ff(i, i) = -y0cc[si];
    // XXX g-parts!
    V0.gg(i, i) = -y0cc[si]; // check! XXX
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
        V0.gg(i, i) += coef * ykcc[si]; // check! XXX
      }
    }
  }

  V0.drj_in_place();

  GMatrix OneNegPc(m_i0, m_stride, m_subgrid_points, true, m_HF->grid_sptr());
  const auto &core = m_HF->core();
  for (std::size_t ia = 0; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    if (a.kappa() != kappa)
      continue;
    OneNegPc.add(a, a, -1.0);
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

  // G(ri,rj) = xI(ri) x0^†(rj) when ri > rj
  // G(ri,rj) = x0(ri) xI^†(rj) when ri < rj
  // ==>, when rj < ri
  // G(rj,ri) = x0(rj) xI^†(ri)
  //          = (xI(ri) x0^†(rj))^†
  //          = G(ri,rj)^†
  // ^† refers to spinor space

  GMatrix g0I(m_i0, m_stride, m_subgrid_points, true, m_grid);

  const auto winv = 1.0 / w;

  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = g0I.index_to_fullgrid(i);
    for (auto j = 0ul; j <= i; ++j) {
      // i >= j
      const auto sj = g0I.index_to_fullgrid(j);

      g0I.ff(i, j) = xI.f(si) * x0.f(sj) * winv;
      g0I.ff(j, i) = g0I.ff(i, j);

      // G parts: -- Double Check! XXX
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
  ComplexGMatrix g0I(m_i0, m_stride, m_subgrid_points, true, m_grid);

  const auto winv = 1.0 / w;

  // XXX Doesn't work, though think it used to??

  const auto I = std::complex{0.0, 1.0};

  // for (auto i = 0ul; i < m_subgrid_points; ++i) {
  //   const auto si = g0I.index_to_fullgrid(i);
  //   for (auto j = 0ul; j <= i; ++j) {
  //     const auto sj = g0I.index_to_fullgrid(j);
  //     g0I.ff(i, j) =
  //         (x0.f(sj) + I * Ix0.f(sj)) * (xI.f(si) + I * IxI.f(si)) * winv;
  //     // g0I is symmetric
  //     g0I.ff(j, i) = g0I.ff(i, j);
  //   }
  // }

  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = g0I.index_to_fullgrid(i);
    for (auto j = 0ul; j <= i; ++j) {
      const auto sj = g0I.index_to_fullgrid(j);

      // j <= i
      // const auto x0f = x0.f(sj) + I * Ix0.f(sj);
      // const auto x0g = x0.g(sj) + I * Ix0.g(sj);
      // const auto xIf = xI.f(si) + I * IxI.f(si);
      // const auto xIg = xI.g(si) + I * IxI.g(si);

      // Large-large
      // g0I.ff(i, j) = xIf * x0f * winv;
      g0I.ff(i, j) =
          (x0.f(sj) + I * Ix0.f(sj)) * (xI.f(si) + I * IxI.f(si)) * winv;
      g0I.ff(j, i) = std::conj(g0I.ff(i, j)); //?

      // Large-small
      // g0I.fg(i, j) = xIf * x0g * winv;
      // g0I.gf(j, i) = g0I.fg(i, j);

      // // Small-large
      // g0I.gf(i, j) = xIg * x0f * winv;
      // g0I.fg(j, i) = g0I.gf(i, j);

      // // Small-small
      // g0I.gg(i, j) = xIg * x0g * winv;
      // g0I.gg(j, i) = g0I.gg(i, j);
    }
  }

  return g0I;
}

//==============================================================================
ComplexRMatrix Feynman::polarisation_k(int k, std::complex<double> omega,
                                       bool hole_particle) const {

  // polarisation operator is ~ Fa^† * [Gex(ea + w) + Gex(ea - w)] * Fa

  ComplexRMatrix pi_k(m_i0, m_stride, m_subgrid_points, m_grid);

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

    for (int in = 0; in <= m_max_ki; ++in) {
      const auto kn = Angular::kappaFromIndex(in);
      const auto ck_an = Angular::Ck_kk(k, Fa.kappa(), kn);
      if (ck_an == 0.0)
        continue;
      const double c_ang = ck_an * ck_an / double(2 * k + 1);

      ComplexGMatrix Gx_pm = green_excited(kn, ea_minus_w, Fa_hp) +
                             green_excited(kn, ea_plus_w, Fa_hp);

      // loop over coordinate indices.
      // pi is symmetric in r1, r2 so is there more efficient way to do this
      // so that I don't loop over all ri, rj?
      // pi ~ Fa^†(r1)[Gex(r1,r2,ea-w) + Gex(r1,r2,ea+w)]Fa(r2)
      for (auto i = 0ul; i < m_subgrid_points; ++i) {
        const auto si = Gx_pm.index_to_fullgrid(i);
        for (auto j = 0ul; j <= i; ++j) {
          const auto sj = Gx_pm.index_to_fullgrid(j);
          pi_k(i, j) += c_ang * (Fa.f(si) * Gx_pm.ff(i, j) * Fa.f(sj) +
                                 Fa.g(si) * Gx_pm.gf(i, j) * Fa.f(sj) +
                                 Fa.f(si) * Gx_pm.fg(i, j) * Fa.g(sj) +
                                 Fa.g(si) * Gx_pm.gg(i, j) * Fa.g(sj));
        }
      }
    }
  }

  // Fill symmetric lower half:
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    for (auto j = 0ul; j <= i; ++j) {
      pi_k(j, i) = pi_k(i, j);
    }
  }

  pi_k *= Iunit;
  return pi_k;
}

//==============================================================================
Grid Feynman::form_w_grid(double w0, double wratio) const {

  // Find max core energy: (for w_max)
  auto wmax_core = 30.0;
  const auto &core = m_HF->core();
  for (const auto &Fc : core) {
    if (Fc.n() < m_min_core_n)
      continue;
    if (std::abs(Fc.en()) > wmax_core)
      wmax_core = std::abs(Fc.en());
  }

  // Target for maximum Im(w): based on core energy.
  // XXX Check this - seems to matter (a tiny bit)
  const auto wmax_t = 2.0 * wratio * wmax_core;

  // Solve wmax < w0 * ratio^{N-1} for N
  const std::size_t wsteps =
      std::size_t(std::log(wratio * wmax_t / w0) / std::log(wratio)) + 1;

  // actual w0, to keep w_ratio exact
  const auto wmax = w0 * std::pow(wratio, int(wsteps - 1));

  return Grid(w0, wmax, wsteps, GridType::logarithmic);
}

//==============================================================================
void Feynman::form_qpiq() {
  std::cout << "Forming QPQ(w,k)";
  if (m_hole_particle || m_screen_Coulomb) {
    std::cout << " (w/ " << (m_screen_Coulomb ? "scr" : "")
              << (m_hole_particle && m_screen_Coulomb ? " + " : "")
              << (m_hole_particle ? "hp" : "") << ")";
  }
  std::cout << " .. " << std::flush;

  const auto num_ks = std::size_t(m_max_k + 1);
  const auto num_ws = m_wgrid.num_points();

  m_qpiq_wk.resize(num_ws, num_ks,
                   ComplexRMatrix{m_i0, m_stride, m_subgrid_points, m_grid});

#pragma omp parallel for collapse(2)
  for (auto iw = 0ul; iw < num_ws; ++iw) {
    for (auto k = 0ul; k < num_ks; ++k) {
      const auto omega = std::complex<double>{m_omre, m_wgrid.r(iw)};

      const auto &q = get_qk(int(k)); // has drj
      const auto qdri = q.dri();      // has drj, and dri
      const auto pi = polarisation_k(int(k), omega, m_hole_particle);

      if (m_screen_Coulomb) {
        const auto X = X_screen(pi, qdri);
        m_qpiq_wk[iw][k] = q * pi * X * qdri;
      } else {
        m_qpiq_wk[iw][k] = q * pi * qdri;
      }
    }
  }

  std::cout << " done\n" << std::flush;
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

  GMatrix Sigma(m_i0, m_stride, m_subgrid_points, m_include_G, m_grid);

  constexpr std::complex<double> I{0.0, 1.0};
  const auto num_kappas = std::size_t(m_max_ki + 1);

// Tell OpenMP how to reduce GMatrix
#pragma omp declare reduction(+ : GMatrix : omp_out += omp_in)                 \
    initializer(omp_priv = GMatrix(omp_orig))

#pragma omp parallel for collapse(2) reduction(+ : Sigma)
  for (auto iw = 0ul; iw < m_wgrid.num_points(); iw++) {
    for (auto iB = 0ul; iB < num_kappas; ++iB) {

      const auto omega = std::complex{m_omre, m_wgrid(iw)};

      // Simpson's rule: Implicit ends (integrand zero at w=0 and w>wmax)
      const auto weight = iw % 2 == 0 ? 4.0 / 3 : 2.0 / 3;

      // I, since dw is on imag. grid; 2 from symmetric +/- w
      const auto dw = I * weight * m_wgrid.drdu(iw);

      const auto kB = Angular::kappaFromIndex(int(iB));

      // Silly, but ig gB includes G, then so will gB_QPQ
      const auto gB = m_include_G ? green(kB, env + omega) :
                                    green(kB, env + omega).drop_g();

      for (auto k = 0ul; int(k) <= m_max_k; k++) {

        // For doing single k (tests and for fk factors)
        if (in_k && *in_k != int(k))
          continue;

        const auto ck_vB = Angular::Ck_kk(int(k), kv, kB);
        if (ck_vB == 0.0)
          continue;

        const auto &qpq_dw = m_qpiq_wk[iw][k];

        const auto c_ang_dw =
            dw * ck_vB * ck_vB / double(Angular::twoj_k(kv) + 1);

        Sigma += (c_ang_dw * mult_elements(gB, qpq_dw)).real();
      }
    }
  }

  // Extra 2 from symmetric + / -w
  Sigma *= (m_wgrid.du() / M_PI);

  return Sigma;
}

} // namespace MBPT