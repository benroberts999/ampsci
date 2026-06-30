#include "CorrelationPotential.hpp"
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

namespace MBPT {

//==============================================================================
CorrelationPotential::CorrelationPotential(
  const std::string &fname, const HF::HartreeFock *vHF,
  const std::vector<DiracSpinor> &basis, double r0, double rmax,
  std::size_t stride, int n_min_core, SigmaMethod method, bool include_g,
  bool include_Breit_b2, int n_max_breit, const FeynmanOptions &Foptions,
  bool calculate_fk, const std::vector<double> &fk,
  const std::vector<double> &etak, std::optional<LadderOptions> ladder_opts,
  const std::string &ladder_sigma_file)
  : m_HF(vHF),
    m_basis(basis),
    m_r0(r0),
    m_rmax(rmax),
    m_stride(stride),
    m_i0(m_HF->grid().getIndex(r0)),
    m_size((m_HF->grid().getIndex(rmax) - m_i0) / m_stride + 1),
    m_method(method),
    m_n_min_core(n_min_core),
    m_includeG(include_g),
    m_includeBreit_b2(include_Breit_b2),
    m_n_max_breit(n_max_breit),
    m_Foptions(Foptions),
    m_calculate_fk(calculate_fk),
    m_fk(fk),
    m_etak(etak),
    m_fname(fname),
    m_ladder_opts(std::move(ladder_opts)),
    m_use_ladder(m_ladder_opts.has_value()),
    m_ladder_sigma_file(ladder_sigma_file) {

  std::cout << "\nConstruct Correlation Potential\n";

  // attempt to read in Sigma file:
  // (Just contains Sigma matrix, nothing else)
  const bool read_ok = read_write(fname, IO::FRW::read);

  // attempt to read the (separate) ladder Sigma_L file, if ladder requested and
  // the base Sigma was read (so entries exist to attach Sigma_L to):
  if (m_use_ladder && read_ok)
    read_write_ladder(m_ladder_sigma_file, IO::FRW::read);

  if (!read_ok) {

    if (m_method == SigmaMethod::Feynman) {
      std::cout << "Using Feynman method for direct diagrams, Goldstone "
                   "for exchange\n";
      if (m_calculate_fk && m_Foptions.screening == Screening::include) {
        std::cout << "Calculating f_k from scratch for exchange screening\n";
      } else {
        if (!m_fk.empty()) {
          std::cout << "Exchange screening with: fk = {";
          for (auto &tfk : m_fk) {
            printf("%.3f, ", tfk);
          }
          std::cout << "}\n";
        }
      }
    }

    if (m_method == SigmaMethod::Goldstone) {
      std::cout << "Using Goldstone method for direct/exchnage\n";
      if (!m_fk.empty()) {
        std::cout << "Approximate screening with: fk = {";
        for (auto &tfk : m_fk) {
          printf("%.3f, ", tfk);
        }
        std::cout << "}\n";
      }
      if (!m_etak.empty()) {
        std::cout << "Approx hole-particle: etak = {";
        for (auto &tetak : m_etak) {
          printf("%.3f, ", tetak);
        }
        std::cout << "}\n";
      }
    }

    if (m_includeG) {
      std::cout << "Including G parts of matrix\n";
    }
    if (m_HF->vBreit()) {
      std::cout << "Including one-body Breit (via basis/Green's function)\n";
    }
    if (m_HF->vBreit() && m_includeBreit_b2) {
      std::cout << "Including two-body Breit [B2] correction: up to n="
                << m_n_max_breit << ", using Goldstone\n";
    }

    printf("Sigma sub-grid: r=(%.1e, %.1f)aB with %i points. [i0=%i, "
           "stride=%i]\n",
           vHF->grid().r(m_i0), vHF->grid().r(m_i0 + m_stride * (m_size - 1)),
           int(m_size), int(m_i0), int(m_stride));

    // If didn't read, setup Goldstone/Feynman (create Yk, pol operator etc.)
    m_Gold =
      Goldstone(basis, m_HF->core(), m_i0, m_stride, m_size, n_min_core,
                m_includeG, m_includeBreit_b2 ? m_HF->vBreit() : nullptr);
    if (m_method == SigmaMethod::Feynman) {
      setup_Feynman();
    }
  }
}

//==============================================================================
void CorrelationPotential::formSigma(int kappa, double ev, int n,
                                     const DiracSpinor *Fv) {

  // 1. check if exists. If so, do nothing.

  const auto it =
    std::find_if(m_Sigmas.cbegin(), m_Sigmas.cend(), [kappa, n](const auto &s) {
      return s.kappa == kappa && (s.n == n || n <= 0);
    });
  if (it != m_Sigmas.cend()) {
    // have sigma already!
    // print deets!
    auto de = Fv ? *Fv * (it->Sigma * *Fv) : 0.0;
    fmt::print("Have Sigma: kappa={}, en={:.5f}, de={:.5e}\n", it->kappa,
               it->en, de);
    return;
  }
  if (Fv) {
    assert(Fv->kappa() == kappa);
    fmt::print("Form Σ for {} at e = {:.4f} au = {:.2f} /cm\n",
               Fv->shortSymbol(), ev, ev * PhysConst::Hartree_invcm);
  } else {
    fmt::print("Form Σ for kappa={} at e = {:.4f}\n", kappa, ev);
  }

  auto S = m_method == SigmaMethod::Feynman ? formSigma_F(kappa, ev, Fv) :
                                              formSigma_G(kappa, ev, Fv);

  m_Sigmas.push_back({kappa, ev, std::move(S), n, 1.0});
}

//==============================================================================
void CorrelationPotential::prepare_ladder(const Wavefunction &wf) {
  if (!m_use_ladder || m_Ladder.has_value())
    return;
  // If every stored Sigma already has its Sigma_L (read from file), there is
  // nothing to compute -- don't build the (expensive) Lk/Qk machinery.
  const bool all_cached =
    !m_Sigmas.empty() &&
    std::all_of(m_Sigmas.cbegin(), m_Sigmas.cend(),
                [](const auto &s) { return s.Sigma_L.has_value(); });
  if (all_cached)
    return;
  m_Ladder.emplace(wf, m_r0, m_rmax, m_stride, *m_ladder_opts);
}

//==============================================================================
void CorrelationPotential::formSigma_L(int kappa, double ev, int n,
                                       const DiracSpinor *Fv) {
  if (!m_use_ladder)
    return;

  // Sigma_L is stored alongside the base Sigma for this (kappa, n).
  auto it =
    std::find_if(m_Sigmas.begin(), m_Sigmas.end(), [kappa, n](const auto &s) {
      return s.kappa == kappa && (s.n == n || n <= 0);
    });
  if (it == m_Sigmas.end())
    return;

  if (it->Sigma_L.has_value()) {
    // Read from file -- nothing to compute.
    const auto deL = Fv ? *Fv * (*it->Sigma_L * *Fv) : 0.0;
    fmt::print("Have Sigma_L: kappa={}, de={:.5e}\n", kappa, deL);
    return;
  }

  // prepare_ladder() must be called first
  if (!m_Ladder.has_value())
    return;

  // Build Sigma_L with the same g-inclusion as the rest of Sigma (m_includeG).
  auto SL = m_Ladder->Sigma_ladder(kappa, ev, Fv, m_includeG);
  if (Fv) {
    const auto deL = (*Fv) * (SL * *Fv);
    fmt::print("  de[ladder] = {:.2f}\n", deL * PhysConst::Hartree_invcm);
  }
  it->Sigma_L = std::move(SL);
}

//==============================================================================
GMatrix CorrelationPotential::formSigma_F(int kappa, double ev,
                                          const DiracSpinor *Fv) {

  if (!m_Fy) {
    setup_Feynman();
  }

  std::vector<double> vfk;

  // Calculate screening factors:
  if (m_calculate_fk && m_Fy->screening()) {
    assert(Fv != nullptr && "Cannot calculate fk without Fv");
    vfk = calculate_fk(ev, *Fv);
    std::cout << "  fk   = {";
    for (auto &e : vfk) {
      printf("%.3f, ", e);
    }
    std::cout << "}\n";
    // If not stored, store first screening factors
    if (m_fk.empty()) {
      m_fk = vfk;
    }
  } else {
    vfk = m_fk;
  }

  if (Fv) {
    fmt::print("  de({}) = ", Fv->shortSymbol());
    std::cout << std::flush;
  }

  auto Sd = m_Fy->Sigma_direct(kappa, ev);

  double deD{0.0};
  if (Fv) {
    deD = (*Fv) * (Sd * *Fv);
    fmt::print("{:.2f} + ", deD * PhysConst::Hartree_invcm);
    std::cout << std::flush;
  }

  const auto Sx = m_Gold->Sigma_exchange(kappa, ev, vfk);

  if (Fv) {
    const auto deX = (*Fv) * (Sx * *Fv);
    fmt::print("{:.2f} = {:.2f}\n", deX * PhysConst::Hartree_invcm,
               (deD + deX) * PhysConst::Hartree_invcm);
    std::cout << std::flush;
  }

  if (m_includeBreit_b2) {
    // nb: do some extra work to calculate it seperately (Qk and Pk)..
    // But, since selection rules are different, it's better this way
    fmt::print("  de[B2]  = ");
    std::cout << std::flush;
    const auto dS =
      m_Gold->dSigma_Breit2(kappa, ev, m_fk, m_etak, 99, m_n_max_breit);
    const auto deB2 = (*Fv) * (dS * *Fv);
    fmt::print("{:.2f}\n", deB2 * PhysConst::Hartree_invcm);
    std::cout << std::flush;
    Sd += dS;
  }

  return Sd + Sx;
}

//==============================================================================
std::vector<double>
CorrelationPotential::calculate_fk(double ev, const DiracSpinor &v) const {

  assert(m_Fy0 && m_FyX);

  // "Clamp" screening factors: don't allow |fk| to be too large
  // Occurs very rarely when diagram is very small;
  // Division is numerically unstable, but overall diagram is insignificant
  const auto max_fk = 10.0;
  assert(max_fk > 1.0 && "max_fk must be >1, otherwise meaningless.");
  // flag for printing clamping warning
  int N_clamped = 0;

  std::vector<double> vfk;
  for (int k = 0; k <= 8; ++k) {
    const auto Sd0 = m_Fy0->Sigma_direct(v.kappa(), ev, k);
    const auto SdX = m_FyX->Sigma_direct(v.kappa(), ev, k);
    const auto de0 = v * (Sd0 * v);
    const auto deX = v * (SdX * v);
    auto fk = de0 != 0.0 ? deX / de0 : 1.0;

    // clamp fk:
    if (std::abs(fk) > max_fk) {
      N_clamped++;
      fk = std::copysign(max_fk, fk);
    }

    vfk.push_back(fk);
  }

  if (N_clamped > 0) {
    fmt::print("  (* Warning: clamped {} screening factor to |fk|<={})\n",
               N_clamped, max_fk);
  }
  return vfk;
}

//==============================================================================
std::vector<double>
CorrelationPotential::calculate_etak(double ev, const DiracSpinor &v) const {
  assert(m_Fy0 && m_FyH);
  std::vector<double> vetak;
  for (int k = 0; k <= 6; ++k) {
    // Include screening when calc eta:
    const auto Sd0 = m_FyX->Sigma_direct(v.kappa(), ev, k);
    const auto SdX = m_Fy->Sigma_direct(v.kappa(), ev, k);
    const auto de0 = v * (Sd0 * v);
    const auto deH = v * (SdX * v);
    const auto etak = de0 != 0.0 ? deH / de0 : 1.0;
    vetak.push_back(etak);
    if (std::abs(etak - 1.0) < 0.01 || de0 == 0.0)
      break;
  }
  return vetak;
}

//==============================================================================
GMatrix CorrelationPotential::formSigma_G(int kappa, double ev,
                                          const DiracSpinor *Fv) {

  // faster to calculate direct and exchange together;
  // ...but then you lose info on the relative contributions
  bool exchange_seperately = true;

  if (!m_Gold) {
    m_Gold =
      Goldstone(m_basis, m_HF->core(), m_i0, m_stride, m_size, m_n_min_core,
                m_includeG, m_includeBreit_b2 ? m_HF->vBreit() : nullptr);
  }

  auto Sd = exchange_seperately ?
              m_Gold->Sigma_direct(kappa, ev, m_fk, m_etak) :
              m_Gold->Sigma_both(kappa, ev, m_fk, m_etak);

  double deD{0.0};
  if (Fv) {
    deD = (*Fv) * (Sd * *Fv);
    fmt::print("  de({}) = {:.2f} ", Fv->shortSymbol(),
               deD * PhysConst::Hartree_invcm);
  }

  if (exchange_seperately) {
    const auto Sx = m_Gold->Sigma_exchange(kappa, ev, m_fk);
    if (Fv) {
      const auto deX = (*Fv) * (Sx * *Fv);
      fmt::print("+ {:.2f} = {:.2f}\n", deX * PhysConst::Hartree_invcm,
                 (deD + deX) * PhysConst::Hartree_invcm);
    }
    Sd += Sx;
  } else {
    std::cout << "\n" << std::flush;
  }

  if (m_includeBreit_b2) {
    // nb: do some extra work to calculate it seperately (Qk and Pk)..
    // But, since selection rules are different, it's better this way
    fmt::print("  de[B2]  = ");
    std::cout << std::flush;
    const auto dS =
      m_Gold->dSigma_Breit2(kappa, ev, m_fk, m_etak, 99, m_n_max_breit);
    const auto deB2 = (*Fv) * (dS * *Fv);
    fmt::print("{:.2f}\n", deB2 * PhysConst::Hartree_invcm);
    std::cout << std::flush;
    Sd += dS;
  }

  return Sd;
}

//==============================================================================
void CorrelationPotential::setup_Feynman() {

  if (!m_Gold) {
    // Also need Goldstone for Feynman (exchange)
    m_Gold = Goldstone(m_basis, m_HF->core(), m_i0, m_stride, m_size,
                       m_n_min_core, m_includeG);
  }

  if (!m_Fy) {
    m_Fy = Feynman(m_HF, m_i0, m_stride, m_size, m_Foptions, m_n_min_core,
                   m_includeG, true, m_fname);
  }

  if (m_calculate_fk && !m_Fy0) {

    // Allow slightly larger stride for calculating effective screening
    const auto stride_scale = 1.0;
    const auto t_stride = std::size_t(stride_scale * double(m_stride));
    const auto t_size = (m_HF->grid().getIndex(m_rmax) - m_i0) / t_stride + 1;

    auto t_Foptions0{m_Foptions};
    auto t_FoptionsX{m_Foptions};

    if (m_Fy->screening()) {
      // Fy with no screening (+no hp)
      t_Foptions0.screening = Screening::exclude;
      t_Foptions0.hole_particle = HoleParticle::exclude;
      m_Fy0 = Feynman(m_HF, m_i0, t_stride, t_size, t_Foptions0, m_n_min_core,
                      m_includeG, false, m_fname);

      // Fy with screening (but no hp)
      t_FoptionsX.screening = Screening::include;
      t_FoptionsX.hole_particle = HoleParticle::exclude;
      m_FyX = m_Fy->hole_particle() ?
                Feynman(m_HF, m_i0, t_stride, t_size, t_FoptionsX, m_n_min_core,
                        m_includeG, false, m_fname) :
                m_Fy;
    }
  }
}

//==============================================================================
const SigmaData *CorrelationPotential::get(int kappa, int n) const {
  if (n <= 0) {
    // returns FIRST sigma that has correct kappa, order matters!
    const auto it =
      std::find_if(m_Sigmas.cbegin(), m_Sigmas.cend(),
                   [kappa](const auto &s) { return s.kappa == kappa; });
    return it != m_Sigmas.cend() ? &(*it) : nullptr;
  } else {
    // Find first Sigma that matches kappa _and_ n
    const auto it = std::find_if(
      m_Sigmas.cbegin(), m_Sigmas.cend(),
      [kappa, n](const auto &s) { return s.kappa == kappa && s.n == n; });
    // If not found, look (recursively) for next lowest n
    return it != m_Sigmas.cend() ? &(*it) : get(kappa, n - 1);
  }
}

//==============================================================================
const GMatrix *CorrelationPotential::getSigma(int kappa, int n) const {
  const auto Sig = get(kappa, n);
  return Sig ? &(Sig->Sigma) : nullptr;
}

//==============================================================================
double CorrelationPotential::getLambda(int kappa, int n) const {
  const auto Sig = get(kappa, n);
  return Sig ? (Sig->lambda) : 1.0;
}

//==============================================================================
//! returns Spinor: Sigma|Fv>
//! @details If Sigma for kappa_v doesn't exist, returns |0>.
DiracSpinor CorrelationPotential::SigmaFv(const DiracSpinor &Fv) const {
  const auto Sv = get(Fv.kappa(), Fv.n());
  if (!Sv)
    return 0.0 * Fv;
  auto SF = Sv->Sigma * Fv;
  if (m_use_ladder && Sv->Sigma_L.has_value())
    SF += *Sv->Sigma_L * Fv;
  return Sv->lambda * SF;
}

//==============================================================================
//! Stores scaling factors, lambda, for each kappa (Sigma -> lamda*Sigma)
void CorrelationPotential::scale_Sigma(const std::vector<double> &lambdas) {
  for (std::size_t i = 0; i < m_Sigmas.size() && i < lambdas.size(); ++i) {
    m_Sigmas.at(i).lambda = lambdas.at(i);
  }
}

//==============================================================================
// if n=0, scales _all_
void CorrelationPotential::scale_Sigma(double lambda, int kappa, int n) {
  for (auto &Sig : m_Sigmas) {
    if (Sig.kappa == kappa && (Sig.n == n || n <= 0)) {
      Sig.lambda = lambda;
    }
  }
}

//==============================================================================
//! Prints the scaling factors to screen
void CorrelationPotential::print_scaling() const {

  bool print = false;
  for (const auto &Sig : m_Sigmas) {
    if (std::abs(Sig.lambda - 1.0) > 0.00001)
      print = true;
  }

  if (print) {
    std::cout << "Scaling factors: lambda = ";
    for (const auto &Sig : m_Sigmas) {
      std::cout << Sig.lambda << ", ";
    }
    std::cout << "\n";
  }
}

//==============================================================================
//! Prints the sub-grid parameters to screen
void CorrelationPotential::print_subGrid() const {
  printf("Sigma sub-grid: r=(%.1e, %.1f)aB with %i points. [i0=%i, "
         "stride=%i]\n",
         m_HF->grid().r(m_i0), m_HF->grid().r(m_i0 + m_stride * (m_size - 1)),
         int(m_size), int(m_i0), int(m_stride));
}

//==============================================================================
bool CorrelationPotential::read_write(const std::string &fname,
                                      IO::FRW::RoW rw) {

  if (fname == "false")
    return false;
  if (rw == IO::FRW::read && !IO::FRW::file_exists(fname))
    return false;

  const auto rw_str =
    rw == IO::FRW::write ? "\nWriting to " : "\nReading from ";
  std::cout << rw_str << "Sigma file: " << fname << " ... " << std::flush;

  std::fstream iofs;
  IO::FRW::open_binary(iofs, fname, rw);

  // // write/read some grid parameters - just to check
  {
    double r0 = rw == IO::FRW::write ? m_HF->grid().r0() : 0;
    double rmax = rw == IO::FRW::write ? m_HF->grid().rmax() : 0;
    double b = rw == IO::FRW::write ? m_HF->grid().loglin_b() : 0;
    std::size_t pts = rw == IO::FRW::write ? m_HF->grid().num_points() : 0;
    rw_binary(iofs, rw, r0, rmax, b, pts);
    if (rw == IO::FRW::read) {
      const bool grid_ok = std::abs((r0 - m_HF->grid().r0()) / r0) < 1.0e-6 &&
                           std::abs(rmax - m_HF->grid().rmax()) < 0.001 &&
                           std::abs(b - m_HF->grid().loglin_b()) < 0.001 &&
                           pts == m_HF->grid().num_points();
      if (!grid_ok) {
        std::cout << "\nCannot read from:" << fname << ". Grid mismatch\n"
                  << "Read: " << r0 << ", " << rmax << " w/ N=" << pts
                  << ", b=" << b << ",\n but expected: " << m_HF->grid().r0()
                  << ", " << m_HF->grid().rmax()
                  << " w/ N=" << m_HF->grid().num_points()
                  << ", b=" << m_HF->grid().loglin_b() << "\n";
        std::cout << "Will calculate from scratch, + over-write file.\n";
        return false;
      }
    }
  }

  // Sub-grid:
  rw_binary(iofs, rw, m_r0, m_rmax, m_stride, m_i0, m_size, m_includeG);

  // Number of kappas (number of Sigma/G matrices)
  std::size_t num_Sigmas = rw == IO::FRW::write ? m_Sigmas.size() : 0;
  rw_binary(iofs, rw, num_Sigmas);

  for (std::size_t iS = 0; iS < num_Sigmas; ++iS) {

    if (rw == IO::FRW::read) {
      m_Sigmas.push_back(
        {0, 0.0, GMatrix{m_i0, m_stride, m_size, m_includeG, m_HF->grid_sptr()},
         0, 1.0}); // don't read/write lamba
    }
    auto &Sig = m_Sigmas.at(iS);
    rw_binary(iofs, rw, Sig.kappa, Sig.en, Sig.n);
    auto &Gk = Sig.Sigma;

    assert(Gk.includes_g() == m_includeG);
    assert(Gk.size() == m_size);
    assert(Gk.stride() == m_stride);
    assert(Gk.i0() == m_i0);
    for (auto i = 0ul; i < m_size; ++i) {
      for (auto j = 0ul; j < m_size; ++j) {
        rw_binary(iofs, rw, Gk.ff(i, j));
        if (m_includeG) {
          rw_binary(iofs, rw, Gk.fg(i, j));
          rw_binary(iofs, rw, Gk.gf(i, j));
          rw_binary(iofs, rw, Gk.gg(i, j));
        }
      }
    }
  }

  std::cout << "done.\n";
  if (rw == IO::FRW::read) {
    std::cout << "Read Sigma from file: " << fname << "\n";
    print_info();
  }
  return true;
}

//==============================================================================
bool CorrelationPotential::read_write_ladder(const std::string &fname,
                                             IO::FRW::RoW rw) {
  const bool reading = rw == IO::FRW::read;
  if (fname.empty() || fname == "false")
    return false;
  if (reading && !IO::FRW::file_exists(fname))
    return false;

  const auto rw_str = reading ? "\nReading from " : "\nWriting to ";
  std::cout << rw_str << "ladder Sigma_L file: " << fname << " ... "
            << std::flush;

  std::fstream iofs;
  IO::FRW::open_binary(iofs, fname, rw);

  // Sub-grid (must match the base Sigma):
  auto t_r0 = m_r0, t_rmax = m_rmax;
  auto t_stride = m_stride, t_i0 = m_i0, t_size = m_size;
  bool t_g = m_includeG;
  rw_binary(iofs, rw, t_r0, t_rmax, t_stride, t_i0, t_size, t_g);
  if (reading && (t_stride != m_stride || t_i0 != m_i0 || t_size != m_size ||
                  t_g != m_includeG)) {
    std::cout << "\nLadder Sigma_L file sub-grid mismatch; ignoring.\n";
    return false;
  }

  // Number of stored Sigma_L matrices:
  std::size_t num = 0;
  if (!reading) {
    for (const auto &Sig : m_Sigmas)
      if (Sig.Sigma_L.has_value())
        ++num;
  }
  rw_binary(iofs, rw, num);

  const auto rw_matrix = [&](GMatrix &Gk) {
    for (auto i = 0ul; i < m_size; ++i) {
      for (auto j = 0ul; j < m_size; ++j) {
        rw_binary(iofs, rw, Gk.ff(i, j));
        if (m_includeG) {
          rw_binary(iofs, rw, Gk.fg(i, j));
          rw_binary(iofs, rw, Gk.gf(i, j));
          rw_binary(iofs, rw, Gk.gg(i, j));
        }
      }
    }
  };

  if (!reading) {
    for (auto &Sig : m_Sigmas) {
      if (!Sig.Sigma_L.has_value())
        continue;
      rw_binary(iofs, rw, Sig.kappa, Sig.n);
      rw_matrix(*Sig.Sigma_L);
    }
  } else {
    for (std::size_t iS = 0; iS < num; ++iS) {
      int kappa = 0, n = 0;
      rw_binary(iofs, rw, kappa, n);
      GMatrix Gk{m_i0, m_stride, m_size, m_includeG, m_HF->grid_sptr()};
      rw_matrix(Gk);
      // Attach to the matching base Sigma (silently drop if none):
      auto it = std::find_if(
        m_Sigmas.begin(), m_Sigmas.end(),
        [kappa, n](const auto &s) { return s.kappa == kappa && s.n == n; });
      if (it != m_Sigmas.end())
        it->Sigma_L = std::move(Gk);
    }
  }

  std::cout << "done.\n";
  return true;
}

//==============================================================================
void CorrelationPotential::print_info() const {
  for (const auto &Sig : m_Sigmas) {
    fmt::print("kappa = {}, ev = {:.5f}", Sig.kappa, Sig.en);
    if (Sig.n > 0) {
      fmt::print(" (n = {})", Sig.n);
    }
    if (std::abs(Sig.lambda - 1.0) > 1.0e-16) {
      fmt::print(" ; scaled with λ = {}", Sig.lambda);
    }
    std::cout << "\n";
  }
}

} // namespace MBPT