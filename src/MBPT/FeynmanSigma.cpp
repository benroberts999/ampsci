#include "MBPT/FeynmanSigma.hpp"
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/YkTable.hpp"
#include "DiracODE/DiracODE.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/SafeProfiler.hpp"
#include "LinAlg/LinAlg.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "MBPT/GreenMatrix.hpp"
#include "Maths/Grid.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cassert>
#include <numeric>
#include <optional>

// omp_get_thread_num() is not defined if not using -fopenmp
// Also: helps protect in case that "omp.h" is not available
#if defined(_OPENMP)
#include <omp.h>
constexpr bool use_omp = true;
#else
constexpr bool use_omp = true;
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

//! Many-body perturbation theory
namespace MBPT {

//******************************************************************************
FeynmanSigma::FeynmanSigma(const HF::HartreeFock *const in_hf,
                           const std::vector<DiracSpinor> &basis,
                           const Sigma_params &sigp,
                           const rgrid_params &subgridp,
                           const std::string &fname)
    : CorrelationPotential(in_hf, basis, sigp, subgridp),
      m_screen_Coulomb(sigp.screenCoulomb),
      m_holeParticle(sigp.holeParticle),
      m_omre(sigp.real_omega),
      m_w0(sigp.w0),
      m_w_ratio(sigp.w_ratio),
      m_Green_method(sigp.GreenBasis ? GrMethod::basis : GrMethod::Green),
      m_Pol_method(sigp.PolBasis ? GrMethod::basis : GrMethod::Green),
      p_hf(in_hf),
      m_min_core_n(sigp.min_n_core),
      m_max_kappaindex_core(2 * DiracSpinor::max_l(p_hf->get_core())),
      m_max_kappaindex(2 * sigp.max_l_excited) {

  std::cout << "\nCorrelation potential (Sigma): Feynman\n";

  // note: if l_max (for Green's fn) larger than l in basis, need to increase
  // m_kmax, and therefore extend the 6j and Ck tables
  const auto kmax_ex = m_max_kappaindex + 1;
  if (kmax_ex > m_maxk)
    m_maxk = kmax_ex;

  // io file name:
  const bool read_ok = read_write(fname, IO::FRW::read);

  if (!read_ok) {

    prep_Feynman();
  }
}

//******************************************************************************
void FeynmanSigma::formSigma(int kappa, double en, int n) {
  // Calc dir + exchange
  // Print D, X, (D+X) energy shift
  // Add (D+X) to m_Sigma, and (n,k) to lookup_list
  // most of this is the same between each...?

  // already exists:
  const auto index = getSigmaIndex(n, kappa);
  if (index < m_Sigma_kappa.size()) {
    const auto [n2, k2, en2] = m_nk[index];
    // already have this (exact) potential?
    if (n == n2)
      return;
  }
  // XXX Need to read/write QPQ etc!!! for this to work ?
  // Temporary solution:
  if (m_qpq_wk.size() == 0) {
    prep_Feynman();
  }

  if (p_hf->get_Breit() != nullptr) {
    // Breit was included into HF; we should not do this!
    std::cout
        << "\nWARNING: Trying to calculate Sigma in Feynman method when Breit "
           "is included will result in incorrect results! Calculate Sigma "
           "first (without Breit), then you may include Breit\n";
  }

  m_nk.emplace_back(n, kappa, en);
  auto &Sigma = m_Sigma_kappa.emplace_back(m_subgrid_points, m_include_G);

  // if v.kappa > basis, then Ck angular factor won't exist!
  if (Angular::twoj_k(kappa) > m_yeh.Ck().max_tj()) {
    std::cout << "Warning: angular not good\n";
    return;
  }

  printf("k=%2i at en=%8.5f.. ", kappa, en);
  std::cout << std::flush;

  // find lowest excited state, output <v|S|v> energy shift:
  const auto find_kappa = [kappa, n](const auto &a) {
    return a.k == kappa && (a.n == n || n == 0);
  };
  const auto vk = std::find_if(cbegin(m_excited), cend(m_excited), find_kappa);

  // Exchange part:
  if (m_print_each_k) {
    // TEMPORARY: Print each k for direct part: for testing
    std::cout << "\n";
    const auto max_k = std::min(m_maxk, m_k_cut);
    for (int k = 0; k <= max_k; ++k) {
      const auto Sigma_k = FeynmanDirect(kappa, en, k);

      // Print out the direct energy shift:
      if (vk != cend(m_excited)) {
        const auto deD = *vk * act_G_Fv(Sigma_k, *vk);
        printf(" k=%i de(k)=%9.3f \n", k, deD * PhysConst::Hartree_invcm);
        std::cout << std::flush;
      }

      Sigma += Sigma_k;
    }
  } else {
    Sigma = FeynmanDirect(kappa, en);
  }

  // Exchange part:
  const auto Gmat_X = m_ex_method == ExchangeMethod::none ?
                          0.0 * Sigma :
                          m_ex_method == ExchangeMethod::Goldstone ?
                          Exchange_Goldstone(kappa, en) :
                          m_ex_method == ExchangeMethod::w1 ?
                          FeynmanEx_1(kappa, en) :
                          FeynmanEx_w1w2(kappa, en);

  // Print energy shifts:
  if (vk != cend(m_excited)) {
    auto deD = *vk * act_G_Fv(Sigma, *vk);
    auto deX = *vk * act_G_Fv(Gmat_X, *vk);
    // nb: just approximate (uses splines)
    printf("de= %7.1f + %5.1f = ", deD * PhysConst::Hartree_invcm,
           deX * PhysConst::Hartree_invcm);
    printf("%7.1f", (deD + deX) * PhysConst::Hartree_invcm);
  }

  Sigma += Gmat_X;

  std::cout << "\n";
}

//******************************************************************************
void FeynmanSigma::prep_Feynman() {

  // Extand 6j and Ck
  m_6j.fill(m_maxk); //?
  // m_yeh.extend_Ck(m_maxk); // XXX Check max k OK in Ck tables!?

  if (m_screen_Coulomb)
    std::cout << "Including Coulomb screening\n";

  if (m_holeParticle)
    std::cout << "Including hole-particle interaction\n";

  if (!m_screen_Coulomb && !m_holeParticle)
    std::cout << "Second-order Feynman\n";

  std::cout << "lmax = " << Angular::lFromIndex(m_max_kappaindex) << "\n";
  std::cout << "Using " << ParseEnum(m_Green_method)
            << " method for Green's functions\n";
  std::cout << "Using " << ParseEnum(m_Pol_method)
            << " method for Polarisation operator\n";

  std::cout << "Including from core n=" << m_min_core_n << "\n";

  if (m_holeParticle && m_Pol_method == GrMethod::basis) {
    std::cout << "WARNING: Cannot include hole-particle using basis for "
                 "Polarisation operator!\n --> hp NOT included\n";
  }

  std::cout << "Using " << ParseEnum(m_ex_method) << " for exchange\n";

  // Print effective screening factors (if used)
  if (m_ex_method == ExchangeMethod::Goldstone && !m_fk.empty()) {
    std::cout << "w/ fk=";
    std::for_each(cbegin(m_fk), cend(m_fk),
                  [](auto fk) { std::cout << fk << ", "; });
    std::cout << "\n";
  }

  if (m_include_G)
    std::cout << "(Including FG/GF and GG)\n";

  if (m_Pol_method == GrMethod::basis || m_Green_method == GrMethod::basis ||
      m_ex_method == ExchangeMethod::Goldstone) {
    std::cout << "Basis: " << DiracSpinor::state_config(m_holes) << "/"
              << DiracSpinor::state_config(m_excited) << "\n";
  }

  form_Q_dr();
  form_Vx();
  form_Pa_core();
  setup_omega_grid();

  print_subGrid();

  const auto max_k = std::min(m_maxk, m_k_cut);
  m_qpq_wk = form_QPQ_wk(max_k, m_Pol_method, m_omre, *m_wgridD);
}

//------------------------------------------------------------------------------
void FeynmanSigma::form_Q_dr() {

  std::vector<LinAlg::Matrix<double>> qhat(
      std::size_t(m_maxk + 1), LinAlg::Matrix<double>{m_subgrid_points});

  LinAlg::Matrix<double> t_dri{m_subgrid_points};

  // calculate (real) Qhat and Jacobian
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = ri_subToFull(i);
    const auto dri = dr_subToFull(i);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto sj = ri_subToFull(j);
      const auto drj = dr_subToFull(j);
      const auto rl = p_gr->r()[std::min(si, sj)];
      const auto rm = p_gr->r()[std::max(si, sj)];
      const auto ratio = rl / rm; // = r_< / r_>
      // q0 = (1.0 / rm) = r_<^k / r_>^k+1 ,  for k=0
      qhat[0][i][j] = (1.0 / rm) * dri * drj;
      t_dri[i][j] = dri;
      for (int k = 1; k <= m_maxk; ++k) {
        qhat[std::size_t(k)][i][j] = qhat[std::size_t(k - 1)][i][j] * ratio;
      }
    } // j
  }   // i

  // fill (complex) Qhat
  m_qhat.resize(std::size_t(m_maxk) + 1, {m_subgrid_points, m_include_G});
  for (auto k = 0; k <= m_maxk; ++k) {
    m_qhat[std::size_t(k)].ff = qhat[std::size_t(k)].complex();
    // m_qhat[std::size_t(k)].ff = OldLinAlg::ComplexSqMatrix::make_complex(
    //     {1.0, 0.0}, qhat[std::size_t(k)]);
    if (m_include_G) {
      m_qhat[std::size_t(k)].gg = m_qhat[std::size_t(k)].ff;
    }
  }

  // fill (complex) Jacobian
  m_dri = std::make_unique<ComplexGMatrix>(m_subgrid_points, m_include_G);
  m_drj = std::make_unique<ComplexGMatrix>(m_subgrid_points, m_include_G);
  m_dri->ff = t_dri.complex();
  m_drj->ff = m_dri->ff.transpose();
}

//------------------------------------------------------------------------------
void FeynmanSigma::form_Vx() {
  // Fill Vx:
  const auto max_kappa_index = std::size_t(m_max_kappaindex);
  m_Vxk.resize(max_kappa_index + 1, {m_subgrid_points, m_include_G});

  // If using a local method, exclude exchange here. Note: Only used for
  // testing!
  if (p_hf->excludeExchangeQ())
    return;

  for (auto ik = 0ul; ik <= max_kappa_index; ik++) {
    m_Vxk[ik] = calculate_Vx_kappa(Angular::kappaFromIndex(int(ik)));
  }
}

//------------------------------------------------------------------------------
GMatrix FeynmanSigma::calculate_Vx_kappa(int kappa) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  assert(m_dri != nullptr && "⚠️ form_Q_dr() must be run before form_Vx()");

  GMatrix Vx(m_subgrid_points, m_include_G);
  const auto tj = Angular::twoj_k(kappa);
  const auto kmax = Angular::twojFromIndex(m_max_kappaindex_core);
  const auto &core = p_hf->get_core();

  for (int k = 0; k <= kmax; ++k) {
    GMatrix Vx_k(m_subgrid_points, m_include_G);
    for (const auto &a : core) {
      const auto ck = Angular::Ck_kk(k, kappa, a.k);
      if (ck == 0.0)
        continue;
      const auto c_ang = -1.0 * ck * ck / double(tj + 1);
      addto_G(&Vx_k, a, a, c_ang);
    }
    Vx += Vx_k.mult_elements_by(get_qk(k).get_real());
  }

  return Vx;
}

//------------------------------------------------------------------------------
GMatrix FeynmanSigma::calculate_Vhp(const DiracSpinor &Fc) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // Hole-particle interaction, extra potential
  // This works well for k=0, but works worse than using local pot for k=1
  // (k=1 is most important term!)

  GMatrix V0(m_subgrid_points, m_include_G);

  const auto y0cc = Coulomb::yk_ab(Fc, Fc, 0);
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = ri_subToFull(i);
    V0.ff[i][i] = -y0cc[si];
    if (m_include_G)
      V0.gg[i][i] = V0.ff[i][i];
  }

  V0.mult_elements_by(m_drj->get_real());

  // return V0;

  GMatrix OneNegPc(m_subgrid_points, m_include_G);
  const auto &core = p_hf->get_core();
  for (const auto &Fa : core) {
    if (Fa.k != Fc.k)
      continue;
    // const auto tjp1 = Fa.twoj() + 1;
    addto_G(&OneNegPc, Fa, Fa, -1.0);
  }

  //*****
  OneNegPc.mult_elements_by(m_drj->get_real());
  OneNegPc.plusIdent(); // (1-P)
  auto out = OneNegPc * V0 * OneNegPc;

  // out.mult_elements_by(m_dri->get_real());

  //*****
  // OneNegPc.mult_elements_by(m_drj->get_real());
  // auto OneNegPcL = OneNegPc;
  // OneNegPcL.mult_elements_by(m_dri->get_real());
  // OneNegPc.plusIdent();  // (1-P)
  // OneNegPcL.plusIdent(); // (1-P)
  // auto out = OneNegPcL * V0 * OneNegPc;

  return out;
}

//------------------------------------------------------------------------------
void FeynmanSigma::form_Pa_core() {
  // Fill core |a><a|
  const auto &core = p_hf->get_core();
  m_Pa.resize(core.size(), {m_subgrid_points, m_include_G});
  for (auto ia = 0ul; ia < core.size(); ia++) {
    m_Pa[ia] = G_single(core[ia], core[ia], std::complex<double>{1.0, 0.0});
  }
}

//------------------------------------------------------------------------------
void FeynmanSigma::setup_omega_grid() {
  // Set up imaginary frequency grid:
  std::cout << "Re(w) = " << m_omre << "\n";

  // Find max core energy: (for w_max)
  auto wmax_core = 30.0; // don't let it go below 50
  const auto &core = p_hf->get_core();
  for (const auto &Fc : core) {
    if (Fc.n < m_min_core_n)
      continue;
    if (std::abs(Fc.en()) > wmax_core)
      wmax_core = std::abs(Fc.en());
  }

  // Make inputs:
  const auto w0 = m_w0;
  const double wratio = m_w_ratio;

  // allow -ve Im(w) grid [for testing]
  if (w0 < 0)
    wmax_core *= -1;

  // maximum Im(w): based on core energy.
  const auto wmax = 2.0 * wratio * wmax_core;

  const std::size_t wsteps = Grid::calc_num_points_from_du(
      w0, wmax, std::log(wratio), GridType::logarithmic);
  m_wgridD = std::make_unique<Grid>(w0, wmax, wsteps, GridType::logarithmic);
  std::cout << "Im(w) " << m_wgridD->gridParameters();
  printf(". r=%.2f\n", m_wgridD->r()[1] / m_wgridD->r()[0]);
  if (m_wX_stride != 1 && (m_ex_method == ExchangeMethod::w1 ||
                           m_ex_method == ExchangeMethod::w1w2)) {
    std::cout << "Exchange Im(w) uses stride: " << m_wX_stride << "\n";
  }

  // for (const auto &w : m_wgridD->r() ) {
  //   std::cout << w << ", ";
  // }
  // std::cout << "\n";
}

//******************************************************************************
const GMatrix &FeynmanSigma::get_Vx_kappa(int kappa) const {
  const auto kappa_index = std::size_t(Angular::indexFromKappa(kappa));
  assert(kappa_index < m_Vxk.size());
  return m_Vxk[kappa_index];
}

const ComplexGMatrix &FeynmanSigma::get_qk(int k) const {
  assert(std::size_t(k) < m_qhat.size());
  return m_qhat[std::size_t(k)];
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
ComplexGMatrix FeynmanSigma::Green(int kappa, std::complex<double> en,
                                   States states, GrMethod method) const {
  if (states == States::core) {
    return Green_core(kappa, en);
  } else if (states == States::excited) {
    return Green_ex(kappa, en, method);
  }
  return (method == GrMethod::Green) ? Green_hf(kappa, en) :
                                       Green_hf_basis(kappa, en);
}

//******************************************************************************
ComplexGMatrix FeynmanSigma::G_single(const DiracSpinor &ket,
                                      const DiracSpinor &bra,
                                      const std::complex<double> f) const {
  ComplexGMatrix Gmat(m_subgrid_points, m_include_G);
  // const auto [x, iy] = f;
  auto x = f.real();
  auto iy = f.imag();

  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = ri_subToFull(i);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto sj = ri_subToFull(j);
      // Gmat.ff[i][j] = {x * ket.f(si) * bra.f(sj), iy * ket.f(si) *
      // bra.f(sj)};
      Gmat.ff[i][j] = f * ket.f(si) * bra.f(sj);
    } // j
  }   // i

  if (m_include_G) {
    for (auto i = 0ul; i < m_subgrid_points; ++i) {
      const auto si = ri_subToFull(i);
      for (auto j = 0ul; j < m_subgrid_points; ++j) {
        const auto sj = ri_subToFull(j);
        Gmat.fg[i][j] = {x * ket.f(si) * bra.g(sj), iy * ket.f(si) * bra.g(sj)};
        Gmat.gf[i][j] = {x * ket.g(si) * bra.f(sj), iy * ket.g(si) * bra.f(sj)};
        Gmat.gg[i][j] = {x * ket.g(si) * bra.g(sj), iy * ket.g(si) * bra.g(sj)};
      } // j
    }   // i
  }

  return Gmat;
}

//******************************************************************************
ComplexGMatrix FeynmanSigma::Green_core(int kappa,
                                        std::complex<double> en) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__, "complex");
  // G_core = \sum_a |a><a|/(e_r + i*e_i-ea), for all a with a.k = k
  ComplexGMatrix Gcore(m_subgrid_points, m_include_G);

  // loop over HF core, not Sigma core (used in subtraction to get G^excited)
  const auto &core = p_hf->get_core();
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    if (a.k != kappa)
      continue;
    const auto inv_de = 1.0 / (en - std::complex<double>{a.en()});
    Gcore += inv_de * m_Pa[ia]; // Pa = |a><a|
  }
  return Gcore;
}

//------------------------------------------------------------------------------
ComplexGMatrix FeynmanSigma::Green_ex(int kappa, std::complex<double> en,
                                      GrMethod method, const DiracSpinor *Fc_hp,
                                      int k_hp) const {
  ComplexGMatrix Gk(m_subgrid_points, m_include_G);

  if (method == GrMethod::basis) {
    Gk = Green_hf_basis(kappa, en, true);
  } else {
    // Subtract core states, by forceing Gk to be orthogonal to core:
    // Gk -> Gk - \sum_a|a><a|G
    Gk = Green_hf(kappa, en, Fc_hp, k_hp); // - Green_core(kappa, en);
    makeGOrthogCore(&Gk, kappa);
  }

  return Gk;
}

//------------------------------------------------------------------------------
void FeynmanSigma::makeGOrthogCore(ComplexGMatrix *Gk, int kappa) const {
  // Force Gk to be orthogonal to the core states
  const auto &core = p_hf->get_core();
  const auto &drj = get_drj();
  const auto Gk_old = *Gk;
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    if (core[ia].k == kappa)
      *Gk -= mult_elements(m_Pa[ia], drj) * Gk_old;
  }
}

//------------------------------------------------------------------------------
ComplexGMatrix FeynmanSigma::Green_hf_basis(int kappa, std::complex<double> en,
                                            bool ex_only) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  ComplexGMatrix Gc(m_subgrid_points, m_include_G);
  // const auto &core = p_hf->get_core(); // ?????
  const auto &core = m_holes; // p_hf->get_core(); // ?????
  // XXX Should include all states, even below n_min_core?
  // XXX BUT, total greens fn should be orthog!
  // No matter, since we don't actually use this method for green.....
  // (except for excited)
  const auto &ex = m_excited;
  for (const auto orbs : {&core, &ex}) {
    if (ex_only && orbs == &core)
      continue;
    for (const auto &a : *orbs) {
      if (a.k != kappa)
        continue;

      const auto inv_de = 1.0 / (en - std::complex<double>{a.en()});
      Gc += G_single(a, a, inv_de);
    }
  }
  return Gc;
}

//------------------------------------------------------------------------------
ComplexGMatrix FeynmanSigma::Green_hf(int kappa, std::complex<double> en,
                                      const DiracSpinor *Fc_hp,
                                      int k_hp) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  /*
NOTE: k_hp is "dodgy" parameter.
For including hole-particle interaction:
the 'local' method works best for most k's, but the matrix method (with [1-P])
works better for k=0 (and k>=5 ?)
  */

  // Solve DE (no exchange), regular at 0, infinity ("pinf")
  DiracSpinor x0(0, kappa, p_gr);
  DiracSpinor xI(0, kappa, p_gr);
  const auto alpha = p_hf->get_alpha();
  const auto &Hmag = p_hf->get_Hrad_mag(x0.l());

  auto vl = p_hf->get_vlocal(x0.l());
  if (Fc_hp != nullptr && k_hp != 0) {
    // Include hole-particle interaction (simple):
    // This way works better for k=1, but ruins k=0
    auto y0cc = Coulomb::yk_ab(*Fc_hp, *Fc_hp, 0);
    qip::compose(std::minus{}, &vl, y0cc);
  }

  DiracODE::regularAtOrigin(x0, en.real(), vl, Hmag, alpha);
  DiracODE::regularAtInfinity(xI, en.real(), vl, Hmag, alpha);

  // Evaluate Wronskian at ~65% of the way to pinf. Should be inependent of r
  const auto pp = std::size_t(0.65 * double(xI.max_pt()));
  // Not sure why -ve sign here... ??? But needed to agree w/ basis version;
  const auto w = -1.0 * (xI.f(pp) * x0.g(pp) - x0.f(pp) * xI.g(pp)) / alpha;

  // Get G0 (Green's function, without exchange):
  const auto g0 = MakeGreensG0(x0, xI, w);
  auto Vx = get_Vx_kappa(kappa);

  if (Fc_hp != nullptr && k_hp == 0) {
    // Include hole-particle interaction (w/ [1-P]V[1-P]):
    // This way works better for k=0, but ruins k=1
    Vx += calculate_Vhp(*Fc_hp);
  }

  const std::complex<double> one{1.0, 0.0}; // to convert real to complex

  // Include exchange, and imaginary energy part:

  if (en.imag() == 0.0) {
    // G = [1 - G0*Vx]^{-1} * G0 = -[G0*Vx-1]^{-1} * G0
    // nb: much faster to invert _before_ make complex!
    // (but, only if imag. part is zero)
    return (-1.0 * one) * ((g0 * Vx).plusIdent(-1.0).invert() * g0);
  }

  // G0 := G0(re{e}) - no exchange, only real part
  // G(e) = [1 + i*Im{e}*G0 - G0*Vx]^{-1} * G0
  // Note: differential dr is included in Vx (via Q)
  std::complex<double> iw{0.0, en.imag()};
  return ((iw * g0).mult_elements_by(*m_drj) - one * (g0 * Vx))
             .plusIdent(1.0)
             .invert() *
         (one * g0);
}

//------------------------------------------------------------------------------
ComplexGMatrix FeynmanSigma::GreenAtComplex(const ComplexGMatrix &Gr,
                                            double om_imag) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // Given G(wr) and wi, returns G(wr+i*wi)
  // G(w) =  G(re(w)+im(w)) ;  Gr = G(re(w)), G = G(w),   im(w) = wi
  // G = Gr * [1 + i*wi*Gr]^-1 = Gr * iX
  const std::complex<double> iw{0.0, om_imag};
  return ((iw * Gr).mult_elements_by(*m_drj).plusIdent(1.0).invert()) * Gr;
}

//------------------------------------------------------------------------------
GMatrix FeynmanSigma::MakeGreensG0(const DiracSpinor &x0, const DiracSpinor &xI,
                                   const double w) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // Takes sub-grid into account; ket,bra are on full grid, G on sub-grid
  // G(r1,r2) = x0(rmin)*xI(imax)/w
  GMatrix g0I(m_subgrid_points, m_include_G);

  const auto winv = 1.0 / w;

  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = ri_subToFull(i);
    for (auto j = 0ul; j <= i; ++j) { // j <= i
      const auto sj = ri_subToFull(j);
      g0I.ff[i][j] = x0.f(sj) * xI.f(si) * winv;
      // g0I is symmetric
      g0I.ff[j][i] = g0I.ff[i][j];
    } // j
  }   // i

  if (m_include_G) {
    for (auto i = 0ul; i < m_subgrid_points; ++i) {
      const auto si = ri_subToFull(i);
      for (auto j = 0ul; j < m_subgrid_points; ++j) {
        const auto sj = ri_subToFull(j);
        const auto irmin = std::min(sj, si);
        const auto irmax = std::max(sj, si);
        g0I.fg[i][j] = x0.f(irmin) * xI.g(irmax) * winv;
        // fg = gf?
        g0I.gf[i][j] = x0.g(irmin) * xI.f(irmax) * winv;
        g0I.gg[i][j] = x0.g(irmin) * xI.g(irmax) * winv;
      } // j
    }   // i
  }

  return g0I;
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
ComplexGMatrix FeynmanSigma::Polarisation_k(int k, std::complex<double> omega,
                                            GrMethod method) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  ComplexGMatrix pi_k(m_subgrid_points, m_include_G);

  const auto Iunit = std::complex<double>{0.0, 1.0};
  const auto &core = p_hf->get_core();
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    if (a.n < m_min_core_n)
      continue;

    const auto &pa = m_Pa[ia]; // |a><a|
    const auto ea_minus_w = std::complex<double>{a.en()} - omega;
    const auto ea_plus_w = std::complex<double>{a.en()} + omega;

    const auto *Fa_hp = m_holeParticle ? &a : nullptr;

    for (int in = 0; in <= m_max_kappaindex; ++in) {
      const auto kn = Angular::kappaFromIndex(in);
      const auto ck_an = Angular::Ck_kk(k, a.k, kn);
      if (ck_an == 0.0)
        continue;
      const double c_ang = ck_an * ck_an / double(2 * k + 1);

      pi_k += c_ang * ((Green_ex(kn, ea_minus_w, method, Fa_hp, k) +
                        Green_ex(kn, ea_plus_w, method, Fa_hp, k))
                           .mult_elements_by(pa));
    }
  }
  pi_k *= Iunit;
  return pi_k;
}

//******************************************************************************
ComplexGMatrix FeynmanSigma::X_PiQ(const ComplexGMatrix &pik,
                                   const ComplexGMatrix &qk) const {
  // Calculates [1-pi*q]^{-1} for single k (i.e., no sum over k inside)

  const auto Iunit = std::complex<double>{0.0, 1.0};
  ComplexGMatrix X_piq = +1.0 * Iunit * pik * qk;
  // Extra factor of (-i) -- where from ??
  X_piq.plusIdent(1.0).invert();
  return X_piq;
}

//******************************************************************************
std::vector<std::vector<ComplexGMatrix>>
FeynmanSigma::make_pi_wk(int max_k, GrMethod pol_method, double omre,
                         const Grid &wgrid) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  const auto num_ks = std::size_t(max_k + 1);
  std::vector<std::vector<ComplexGMatrix>> pi_wk(wgrid.num_points());

#pragma omp parallel for
  for (auto iw = 0ul; iw < wgrid.num_points(); ++iw) {
    const auto omega = std::complex<double>{omre, wgrid.r()[iw]};
    pi_wk[iw].reserve(num_ks);
    for (auto k = 0ul; k < num_ks; ++k) {
      pi_wk[iw].push_back(Polarisation_k(int(k), omega, pol_method));
    }
  }
  return pi_wk;
}

//******************************************************************************
std::vector<std::vector<ComplexGMatrix>>
FeynmanSigma::form_QPQ_wk(int max_k, GrMethod pol_method, double omre,
                          const Grid &wgrid) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // Returns QPQ, function of w and k

  std::vector<std::vector<ComplexGMatrix>> qpq(wgrid.num_points());
  std::cout << "Forming QPQ(w,k) matrix.." << std::flush;

  const auto pi_wk = make_pi_wk(max_k, pol_method, omre, wgrid);
  std::cout << "." << std::flush;

  std::cout << "." << std::flush;

  const auto num_ks = std::size_t(max_k + 1);
#pragma omp parallel for
  for (auto iw = 0ul; iw < wgrid.num_points(); ++iw) {
    qpq[iw].reserve(num_ks);
    for (auto k = 0ul; k < num_ks; ++k) {
      // q*p*q => q*X*p*q, x = [1-Pi*Q]^(-1)
      const auto &qk = get_qk(int(k));
      const auto &pi = pi_wk[iw][k];
      if (m_screen_Coulomb) {
        // This way: works!
        const auto X = X_PiQ(pi, qk);
        qpq[iw].emplace_back(qk * X * pi * qk);
      } else {
        qpq[iw].emplace_back(qk * pi * qk);
      }
    }
  }
  std::cout << "..done\n";
  return qpq;
}

//******************************************************************************
std::vector<std::vector<ComplexGMatrix>>
FeynmanSigma::form_Greens_kapw(int max_kappa_index, GrMethod method,
                               double en_re, const Grid &wgrid) const {
  // G(en_re+iw) for each kappa, w
  // nb: en_re = en_v + omre

  const auto num_kappas = std::size_t(max_kappa_index + 1);
  std::vector<std::vector<ComplexGMatrix>> gs(num_kappas);
#pragma omp parallel for
  for (auto ik = 0ul; ik < num_kappas; ++ik) {
    const auto kappa = Angular::kappaFromIndex(int(ik));
    gs[ik].reserve(wgrid.num_points());
    for (auto iw = 0ul; iw < wgrid.num_points(); iw++) {
      std::complex<double> evpw{en_re, wgrid.r()[iw]};
      gs[ik].push_back(Green(kappa, evpw, States::both, method));
    }
  }
  return gs;
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
GMatrix FeynmanSigma::FeynmanDirect(int kv, double env, int in_k) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  /* TEMPORARY: in_k ionly for testing; useful for comparing to Dzuba
   * though..*/

  GMatrix Sigma(m_subgrid_points, m_include_G);

  const std::complex<double> I{0.0, 1.0};

  // // Set up imaginary frequency grid:
  const double omre = m_omre;
  const auto &wgrid = *m_wgridD;
  const auto max_k = std::min(m_maxk, m_k_cut);

  // Store gAs in advance
  const auto num_kappas = std::size_t(m_max_kappaindex + 1);
  const auto gBs =
      form_Greens_kapw(m_max_kappaindex, m_Green_method, env + omre, wgrid);

  // If Im(w) grid is -ve, we integrate "wrong" way around contour; extra -ve
  const auto sw = wgrid.r(0) > 0.0 ? 1.0 : -1.0;

#pragma omp parallel for
  for (auto iw = 0ul; iw < wgrid.num_points(); iw++) { // for omega integral

    // Simpson's rule: Implicit ends (integrand zero at w=0 and w>wmax)
    const auto weight = iw % 2 == 0 ? 4.0 / 3 : 2.0 / 3;

    // I, since dw is on imag. grid; 2 from symmetric +/- w
    const auto dw = I * weight * wgrid.drdu()[iw];

    for (auto k = 0ul; int(k) <= max_k; k++) {

      // For testing only:
      if (in_k >= 0 && in_k != int(k))
        continue;

      const auto qpq_dw = dw * m_qpq_wk[iw][k];

      for (auto iB = 0ul; iB < num_kappas; ++iB) {
        const auto kB = Angular::kappaFromIndex(int(iB));
        const auto ck_vB = Angular::Ck_kk(int(k), kv, kB);
        if (ck_vB == 0.0)
          continue;

        const auto c_ang = ck_vB * ck_vB / double(Angular::twoj_k(kv) + 1);
        const auto C_gB_QPQ_dw =
            c_ang * (mult_elements(gBs[iB][iw], qpq_dw)).get_real();

#pragma omp critical(sum_sigma_d)
        { Sigma += C_gB_QPQ_dw; }

      } // k
    }   // omega
  }     // beta

  // Extra 2 from symmetric + / -w
  Sigma *= (2.0 * sw * wgrid.du() / (2 * M_PI));

  // devide through by dri, drj [these included in q's, but want
  // differential operator for sigma] or.. include one of these in
  // definition of opertion S|v> ?
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto dri = dr_subToFull(i);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto drj = dr_subToFull(j);
      Sigma.ff[i][j] /= (dri * drj);
      if (m_include_G) {
        Sigma.fg[i][j] /= (dri * drj);
        Sigma.gf[i][j] /= (dri * drj);
        Sigma.gg[i][j] /= (dri * drj);
      }
    }
  }

  return Sigma;
} // namespace MBPT

//******************************************************************************
GMatrix FeynmanSigma::FeynmanEx_w1w2(int kv, double en_r) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  // Note: Takes really long time, and doesn't work! XXX

  GMatrix Sx(m_subgrid_points, m_include_G);

  // Set up imaginary frequency grid:
  const double omre = m_omre;
  const auto &wgrid = *m_wgridD;
  const std::complex<double> I{0.0, 1.0};
  const std::complex<double> en{en_r, 0.0};
  const auto tjvp1 = Angular::twoj_k(kv) + 1;

  const auto max_k = std::min(m_maxk, m_k_cut);

  // If Im(w) grid is -ve, we integrate "wrong" way around contour; extra -ve
  const auto sw = wgrid.r(0) > 0.0 ? 1.0 : -1.0;

  std::cout << std::endl;

  const auto num_kappas = std::size_t(m_max_kappaindex + 1);

  // const std::size_t num_para_threads =
  //     use_omp ? num_kappas * num_kappas / 4 : 1;
  const std::size_t num_para_threads =
      use_omp ? std::min(num_kappas * num_kappas,
                         std::size_t(4 * omp_get_max_threads())) :
                1;

  // Store parts of Sx seperately, for more efficient parallelisation
  std::vector<GMatrix> Sxs(num_para_threads, {m_subgrid_points, m_include_G});

  const auto wmax = 40.0;

#pragma omp parallel for num_threads(num_para_threads) collapse(2)
  for (auto iA = 0ul; iA < num_kappas; ++iA) {   // alpha
    for (auto iB = 0ul; iB < num_kappas; ++iB) { // beta

      const auto tid = std::size_t(omp_get_thread_num());
      auto &Sc_i = Sxs[tid];

      for (auto iG = 0ul; iG < num_kappas; ++iG) { // gamma
        const auto [kA, kB, kG] = qip::apply_to(
            Angular::kappaFromIndex, std::array{int(iA), int(iB), int(iG)});

        // The part where the imaginary part of w1 and w2 have same sign is
        // symmetric, and the part where they have opposite sign is also
        // symmetric Therefore, calculate 2*[X(w1,w2)+X(w1,-w2)].
        // w2 here refers just to imaginary part! real part is still + 2 is
        // included below
        for (auto iw1 = 0ul; iw1 < wgrid.num_points(); iw1 += m_wX_stride) {
          const auto dw1 = wgrid.drdu()[iw1];
          const std::complex<double> evpw1{en_r + omre, wgrid.r()[iw1]};

          if (std::abs(wgrid.r()[iw1]) > wmax)
            continue;

          const auto &gA = Green(kA, evpw1, States::both, m_Green_method);

          for (auto iw2 = 0ul; iw2 < wgrid.num_points(); iw2 += m_wX_stride) {
            const auto dw2 = wgrid.drdu()[iw2]; //-ve for -Im(w)
            const std::complex<double> ev_p_w2{en_r + omre, +wgrid.r()[iw2]};
            const std::complex<double> ev_m_w2{en_r + omre, -wgrid.r()[iw2]};

            if (std::abs(wgrid.r()[iw2]) > wmax)
              continue;

            // This seems to be a very expensive random number generator...

            // note: sumkl_gqgqg is the slow part: generateing Green's
            // function here is very ineficient, but doesn't actually take any
            // longer
            // if (std::abs((evpw1 + ev_p_w2).cimag()) <= 1.5 * wmax)
            {
              // (w1 + w2) case:
              const auto &gB_p =
                  Green(kB, evpw1 + ev_p_w2, States::both, m_Green_method);
              const auto &gG_p =
                  Green(kG, ev_p_w2, States::both, m_Green_method);
              const auto gqgqg_p =
                  sumkl_gqgqg(gA, gB_p, gG_p, kv, kA, kB, kG, max_k);
              Sc_i += (dw1 * dw2) * (gqgqg_p);
            }

            // nb: I thought these should be Sc_i -= ...
            // because we have (i*i)=-1 from dw1*dw2, and sumkl_gqgqg is real
            // But, I included an extra i inside sumkl_gqgqg instead (i.e.,
            // assume only get 1 i)?
            // Still, factor of 10 too big!?
            // if (std::abs((evpw1 + ev_m_w2).cimag()) <= 1.5 * wmax)
            {
              // (w1 - w2) case:
              const auto &gB_m =
                  Green(kB, evpw1 + ev_m_w2, States::both, m_Green_method);
              const auto &gG_m =
                  Green(kG, ev_m_w2, States::both, m_Green_method);

              const auto gqgqg_m =
                  sumkl_gqgqg(gA, gB_m, gG_m, kv, kA, kB, kG, max_k);
              // -ve for 'm', since we go wrong direction around w2 contour??
              Sc_i += (dw1 * (-dw2)) * (gqgqg_m);
            }

            // im(gqgqg_m) is small, but real part is v. large
            // im(gqgqg_p) is roughly ok, re(gqgqg_p) is large
            // Note: sumkl_gqgqg already returns only real part.

          } // w2
        }   // w1
      }     // gamma
    }       // beta
  }         // alpha

  for (const auto &Sc_i : Sxs) {
    Sx += Sc_i;
  }

  // sw?
  const double dw_const = sw * double(m_wX_stride) * wgrid.du() / (2.0 * M_PI);
  Sx *= 2.0 * (dw_const * dw_const / tjvp1);

  // devide through by dri, drj [these included in q's, but want
  // differential operator for sigma] or.. include one of these in
  // definition of opertion S|v> ?
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto dri = dr_subToFull(i);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto drj = dr_subToFull(j);
      Sx.ff[i][j] /= (dri * drj);
      // no g?
    }
  }

  return Sx;
}

//******************************************************************************
GMatrix FeynmanSigma::FeynmanEx_1(int kv, double env) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  GMatrix Sx1(m_subgrid_points, m_include_G);

  // Set up imaginary frequency grid:
  const double omre = m_omre; // does seem to depend on this..
  const auto &wgrid = *m_wgridD;
  const auto &core = p_hf->get_core();
  const auto tjvp1 = Angular::twoj_k(kv) + 1;

  // If Im(w) grid is -ve, we integrate "wrong" way around contour; extra -ve
  const auto sw = wgrid.r(0) > 0.0 ? 1.0 : -1.0;

  // Store gAs in advance, since these depend only on w and kA, not a, kB or k
  const auto num_kappas = std::size_t(m_max_kappaindex + 1);
  const auto gAs =
      form_Greens_kapw(m_max_kappaindex, m_Green_method, env + omre, wgrid);

  // // const std::size_t num_para_threads = 12;
  // const std::size_t num_para_threads =
  //     use_omp ? num_kappas * wgrid.num_points() / m_wX_stride / 4 : 1;
  const std::size_t num_para_threads =
      use_omp ? std::min(num_kappas * wgrid.num_points() / m_wX_stride,
                         std::size_t(4 * omp_get_max_threads())) :
                1;

  std::vector<GMatrix> Sx_k(num_para_threads, {m_subgrid_points, m_include_G});

  const auto wmax = 100.0; // XXX Temp?

// #pragma omp parallel for collapse(2)
#pragma omp parallel for num_threads(num_para_threads) collapse(2)
  for (auto iB = 0ul; iB < num_kappas; ++iB) {
    for (auto iw = 0ul; iw < wgrid.num_points(); iw += m_wX_stride) {
      const auto kB = Angular::kappaFromIndex(int(iB));

      if (std::abs(wgrid.r()[iw]) > wmax)
        continue;

      const auto tid = std::size_t(omp_get_thread_num());

      auto omim = wgrid.r()[iw]; // XXX Symmetric?? Or Not??
      const auto omega = std::complex<double>{omre, omim};
      const auto dw1 = wgrid.drdu()[iw]; // rest in 'factor'

      const auto *const qpqw_k = m_screen_Coulomb ? &m_qpq_wk[iw] : nullptr;

      for (auto ia = 0ul; ia < core.size(); ++ia) {
        const auto &Fa = core[ia];
        if (Fa.n < m_min_core_n)
          continue;
        const auto &pa = m_Pa[ia];
        const auto ea = std::complex<double>{Fa.en(), 0.0};

        // hp here too?
        const auto gxBm = Green_ex(kB, ea - omega, m_Green_method);
        const auto gxBp = Green_ex(kB, ea + omega, m_Green_method);

        for (auto iA = 0ul; iA < num_kappas; ++iA) {
          const auto kA = Angular::kappaFromIndex(int(iA));

          const auto &gA = gAs[iA][iw];

          const auto gqpg =
              sumkl_GQPGQ(gA, gxBm, gxBp, pa, kv, kA, kB, Fa.k, qpqw_k);

          Sx_k[tid] += dw1 * gqpg;

        } // alpha
      }   // a
    }
  }

  Sx1 = std::accumulate(Sx_k.cbegin(), Sx_k.cend(), Sx1);

  const auto factor =
      -2.0 * sw * double(m_wX_stride) * wgrid.du() / (2.0 * M_PI) / tjvp1;
  Sx1 *= factor;

  // devide through by dri, drj [these included in q's, but want
  // differential operator for sigma] or.. include one of these in
  // definition of opertion S|v> ?
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto dri = dr_subToFull(i);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto drj = dr_subToFull(j);
      Sx1.ff[i][j] /= (dri * drj);
      // No g
    }
  }

  return Sx1;
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
GMatrix FeynmanSigma::sumkl_gqgqg(const ComplexGMatrix &gA,
                                  const ComplexGMatrix &gB,
                                  const ComplexGMatrix &gG, int kv, int kA,
                                  int kB, int kG, int kmax) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // EXCHANGE part, used in w1w2 version

  auto gqgqg = GMatrix(m_subgrid_points, m_include_G);

  for (int k = 0; k <= kmax; k++) { // k (k1)
    const auto &qk = get_qk(k);

    for (int l = 0; l <= kmax; l++) { // l (k2)
      const auto &ql = get_qk(l);

      // tensor_5_product adds the real part of below to result
      // Sum_ij [ factor * a1j * bij * cj2 * (d_1i * e_i2) ]
      const auto Lkl = Lkl_abcd(k, l, kv, kB, kA, kG);
      if (Lkl != 0.0) {
        const auto s = Angular::neg1pow(k + l);
        // const auto sLkl = std::complex<double>{s * Lkl, 0.0};
        // // XXX extra factor of i ??: - pretty sure this is wrong
        const auto sLkl = std::complex<double>{0.0, s * Lkl};
        tensor_5_product(&gqgqg, sLkl, qk, gB, gG, gA, ql);
      }
    } // l
  }   // k

  return gqgqg;
}

//******************************************************************************
GMatrix FeynmanSigma::sumkl_GQPGQ(
    const ComplexGMatrix &gA, const ComplexGMatrix &gxBm,
    const ComplexGMatrix &gxBp, const ComplexGMatrix &pa, int kv, int kA,
    int kB, int ka, const std::vector<ComplexGMatrix> *const qpqw_k) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // EXCHANGE part, used in w1 version
  // (w2 was integrated over analytically)

  GMatrix sum_GQPG{m_subgrid_points, m_include_G};
  const std::complex<double> I{0.0, 1.0}; // to convert real to complex

  const auto kmax = std::min(m_maxk, m_k_cut);

  for (auto k = 0; k <= kmax; ++k) {
    // qk -> qk - 2i * QPxQ
    const auto &tqk = get_qk(k);
    // Screen q^k(w1) {not checked}
    const auto qk =
        qpqw_k != nullptr ? tqk - 2.0 * I * (*qpqw_k)[std::size_t(k)] : tqk;

    for (auto l = 0; l <= kmax; ++l) {
      const auto &ql = get_qk(l);

      const auto s0 = Angular::neg1pow(k + l);

      const auto L1 = Lkl_abcd(k, l, kv, kB, kA, ka);
      const auto L2 = Lkl_abcd(k, l, kv, ka, kA, kB);

      const auto ic1 = std::complex<double>(s0 * L1, 0.0);
      const auto ic2 = std::complex<double>(s0 * L2, 0.0);

      // tensor_5_product adds the real part of below to result
      // Sum_ij [ factor * a1j * bij * cj2 * (d_1i * e_i2) ]

      if (L1 != 0.0)
        tensor_5_product(&sum_GQPG, ic1, qk, gxBp, pa, gA, ql);

      if (L2 != 0.0)
        tensor_5_product(&sum_GQPG, ic2, qk, pa, gxBm, gA, ql);

    } // l
  }   // k

  return sum_GQPG;
}

//******************************************************************************
void FeynmanSigma::tensor_5_product(
    GMatrix *result, const std::complex<double> &factor,
    const ComplexGMatrix &a, const ComplexGMatrix &b, const ComplexGMatrix &c,
    const ComplexGMatrix &d, const ComplexGMatrix &e) const {
  // Adds real part of below to result
  // Sum_ij [ factor * a1j * bij * cj2 * (d_1i * e_i2) ]
  const auto size = result->size;

  // The a*c part depends only on j (not i).
  // Doing this mult early saves factor of 'size' complex multiplications
  // Leads to a >2x speed-up!
  std::vector<std::complex<double>> ac(size); // see below

  for (auto r1 = 0ul; r1 < size; ++r1) {
    for (auto r2 = 0ul; r2 < size; ++r2) {

      for (auto j = 0ul; j < size; ++j) {
        // early a*c mult
        // Note: a and c are actually real (only half the time..)
        ac[j] = std::complex<double>(a.ff[r1][j]) * c.ff[j][r2];
      }
      std::complex<double> sum_ij{0.0, 0.0};
      for (auto i = 0ul; i < size; ++i) {
        std::complex<double> sum_j{0.0, 0.0};
        for (auto j = 0ul; j < size; ++j) {
          sum_j += ac[j] * b.ff[i][j];
        }
        sum_ij += sum_j * d.ff[r1][i] * e.ff[i][r2];
      }
      result->ff[r1][r2] += (factor * sum_ij).real();

    } // r2
  }   // r1
}

//******************************************************************************
double FeynmanSigma::Lkl_abcd(int k, int l, int ka, int kb, int kc,
                              int kd) const {

  const auto &Ck = m_yeh.Ck();

  const auto Ckac = Ck(k, ka, kc);
  const auto Ckbd = Ck(k, kb, kd);
  const auto Clad = Ck(l, ka, kd);
  const auto Clbc = Ck(l, kb, kc);
  if (Ckac == 0.0 || Ckbd == 0.0 || Clad == 0.0 || Clbc == 0.0)
    return 0.0;

  const auto tja = Angular::twoj_k(ka);
  const auto tjb = Angular::twoj_k(kb);
  const auto tjc = Angular::twoj_k(kc);
  const auto tjd = Angular::twoj_k(kd);
  // return m_6j(tja, tjc, tjb, tjd, k, l) * (Ckac * Ckbd * Clad * Clbc);
  return m_6j.get_2(tja, tjc, 2 * k, tjb, tjd, 2 * l) *
         (Ckac * Ckbd * Clad * Clbc);
}

//******************************************************************************
GMatrix FeynmanSigma::Exchange_Goldstone(const int kappa,
                                         const double en) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  // XXX Update Goldstone to work the same way as this one.
  // And, make me able to call it!!

  // Two second-order exchange diagrams:
  // Diagram (b) (exchange):
  // |Q^k_amn><P^k_amn| / de_amn / [k][j]
  // Diagram (d) (exchange):
  // |Q^k_nba><P^k_nba| / de_nba / [k][j]
  // where:
  // All indeces are summed over,
  // a & b are core states, n & m are virtual excited states,
  // k is multipolarity [Coloulmb expansion]
  // de_xyz = e_v + e_x - e_y - e_z

  GMatrix SxG(m_subgrid_points, m_include_G);

  const auto &Ck = m_yeh.Ck();

  std::vector<GMatrix> Gxs(m_holes.size(), {m_subgrid_points, m_include_G});
#pragma omp parallel for
  for (auto ia = 0ul; ia < m_holes.size(); ia++) {
    const auto &a = m_holes[ia];
    auto &Ga_x = Gxs[ia];
    auto Qkv = DiracSpinor(0, kappa, p_gr); // re-use to reduce alloc'ns
    auto Pkv = DiracSpinor(0, kappa, p_gr); // re-use to reduce alloc'ns
    for (const auto &n : m_excited) {
      const auto [kmin_nb, kmax_nb] = Coulomb::k_minmax(n, a);
      for (int k = kmin_nb; k <= kmax_nb; ++k) {
        if (Ck(k, a.k, n.k) == 0)
          continue;
        const auto f_kkjj = (2 * k + 1) * (Angular::twoj_k(kappa) + 1);
        // const auto &yknb = *m_yeh.get(k, n, a); // check null!

        // Effective screening parameter:
        const auto fk = get_fk(k);

        // Diagrams (b) [exchange]
        for (const auto &m : m_excited) {
          if (Ck(k, kappa, m.k) == 0)
            continue;
          // Coulomb::Qkv_bcd(&Qkv, a, m, n, k, yknb, Ck);
          Qkv = m_yeh.Qkv_bcd(Qkv.k, a, m, n, k);
          // Pkv_bcd_2 allows different screening factor for each 'k2' in
          // exch.
          // Coulomb::Pkv_bcd_2(&Pkv, a, m, n, k, m_yeh(m, a), Ck, m_6j, m_fk);
          // m_yeh.Pkv_bcd_2(&Pkv, a, m, n, k, m_fk);
          Pkv = m_yeh.Pkv_bcd(Pkv.k, a, m, n, k, m_fk);
          const auto dele = en + a.en() - m.en() - n.en();
          const auto factor = fk / (f_kkjj * dele);
          addto_G(&Ga_x, Qkv, Pkv, factor);
        } // m

        // Diagrams (d) [exchange]
        for (const auto &b : m_holes) {
          if (Ck(k, kappa, b.k) == 0)
            continue;
          // Coulomb::Qkv_bcd(&Qkv, n, b, a, k, yknb, Ck);
          Qkv = m_yeh.Qkv_bcd(Qkv.k, n, b, a, k);
          // Coulomb::Pkv_bcd_2(&Pkv, n, b, a, k, m_yeh(n, b), Ck, m_6j, m_fk);
          // m_yeh.Pkv_bcd_2(&Pkv, n, b, a, k, m_fk);
          Pkv = m_yeh.Pkv_bcd(Pkv.k, n, b, a, k, m_fk);
          const auto dele = en + n.en() - b.en() - a.en();
          const auto factor = fk / (f_kkjj * dele); // XXX
          addto_G(&Ga_x, Qkv, Pkv, factor);
        } // b

      } // k
    }   // n
  }     // a

  for (const auto &Gx : Gxs)
    SxG += Gx;

  return SxG;
}

} // namespace MBPT
