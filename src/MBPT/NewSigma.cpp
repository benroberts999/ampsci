#include "NewSigma.hpp"
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "MBPT/RDMatrix.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/format.hpp"
#include <cassert>
#include <iostream>
#include <vector>

namespace MBPT {

//==============================================================================
NewSigma::NewSigma(const std::string &fname, const HF::HartreeFock *vHF,
                   const std::vector<DiracSpinor> &basis, double r0,
                   double rmax, std::size_t stride, int n_min_core,
                   SigmaMethod method, bool include_g,
                   const FeynmanOptions &Foptions, bool calculate_fk,
                   const std::vector<double> &fk,
                   const std::vector<double> &etak)
    : m_HF(vHF),
      m_basis(basis), // dumb
      m_r0(r0),
      m_rmax(rmax),
      m_stride(stride),
      m_i0(m_HF->grid().getIndex(r0)),
      m_size((m_HF->grid().getIndex(rmax) - m_i0) / m_stride + 1),
      m_method(method),
      m_n_min_core(n_min_core),
      m_includeG(include_g),
      m_Foptions(Foptions),
      m_calculate_fk(calculate_fk),
      m_fk(fk),
      m_etak(etak) {

  std::cout << "\nConstruct Correlation Potential\n";

  // attempt to read in Sigma file:
  // (Just contains Sigma matrix, nothing else)
  const bool read_ok = read_write(fname, IO::FRW::read);

  if (!read_ok) {

    if (m_method == SigmaMethod::Feynman) {
      std::cout << "Using Feynman method for direct diagrams, Goldstone "
                   "for exchange\n";
      if (m_includeG) {
        std::cout << "Using Goldstone for G parts of matrix\n";
      }
      if (m_calculate_fk && m_Foptions.screening == Screening::include) {
        std::cout << "Calculating f_k from scratch for exchange screening\n";
        if (m_includeG) {
          std::cout << "Calculating eta_k from scratch for G hole-particle\n";
        }
      } else {
        if (!m_fk.empty()) {
          std::cout << "Exchange screening with: fk = {";
          for (auto &tfk : m_fk) {
            printf("%.3f, ", tfk);
          }
          std::cout << "}\n";
        }
        if (m_includeG && !m_etak.empty()) {
          std::cout << "Approx hole-particle in G: etak = {";
          for (auto &tetak : m_etak) {
            printf("%.3f, ", tetak);
          }
          std::cout << "}\n";
        }
      }
    }

    if (m_method == SigmaMethod::Goldstone) {
      std::cout << "Using Goldstone method for direct/exchnage\n";
      if (m_includeG) {
        std::cout << "Including G parts of matrix\n";
      }
      if (!m_fk.empty()) {
        std::cout << "Approximate screening with: fk = {";
        for (auto &tfk : m_fk) {
          printf("%.3f, ", tfk);
        }
        std::cout << "}\n";
      }
      if (!m_etak.empty()) {
        std::cout << "Approx hole-particle in G: etak = {";
        for (auto &tetak : m_etak) {
          printf("%.3f, ", tetak);
        }
        std::cout << "}\n";
      }
    }

    // If didn't read, setup Goldstone/Feynman (create Yk, pol operator etc.)
    m_Gold = Goldstone(basis, vHF->core(), m_i0, m_stride, m_size, n_min_core,
                       m_includeG);
    if (m_method == SigmaMethod::Feynman) {
      setup_Feynman();
    }
  }
}

//==============================================================================
void NewSigma::formSigma(int kappa, double ev, int n, const DiracSpinor *Fv) {

  // 1. check if exists. If so, do nothing.

  const auto it = std::find_if(
      m_Sigmas.cbegin(), m_Sigmas.cend(), [kappa, n](const auto &s) {
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
    fmt::print("Form Σ for {} at e = {:.4f}\n", Fv->shortSymbol(), ev);
  } else {
    fmt::print("Form Σ for kappa={} at e = {:.4f}\n", kappa, ev);
  }

  auto S = m_method == SigmaMethod::Feynman ? formSigma_F(kappa, ev, Fv) :
                                              formSigma_G(kappa, ev, Fv);

  m_Sigmas.push_back({kappa, ev, std::move(S), n, 1.0});
}

//==============================================================================
GMatrix NewSigma::formSigma_F(int kappa, double ev, const DiracSpinor *Fv) {

  if (!m_Fy) {
    setup_Feynman();
  }

  std::vector<double> vfk;
  std::vector<double> vetak;

  // Calculate screening factors:
  if (m_calculate_fk && m_Fy->screening()) {
    assert(Fv != nullptr && "Cannot calculate fk without having Fv");
    vfk = calculate_fk(ev, *Fv);
    std::cout << "  fk   = {";
    for (auto &e : vfk) {
      printf("%.3f, ", e);
    }
    std::cout << "}\n";
    if (m_includeG) {
      vetak = calculate_etak(ev, *Fv);
      std::cout << "  etak = {";
      for (auto &e : vetak) {
        printf("%.3f, ", e);
      }
      std::cout << "}\n";
    }
  } else {
    vfk = m_fk;
    vetak = m_etak;
  }

  auto Sd = m_Fy->Sigma_direct(kappa, ev);
  if (m_includeG) {
    // If 'Includeing G', calculate using Goldstone technique
    auto Sd2 = m_Gold->Sigma_direct(kappa, ev, vfk, vetak);
    Sd2.ff() = Sd.ff();
    Sd = Sd2;
  }

  double deD{0.0};
  if (Fv) {
    deD = (*Fv) * (Sd * *Fv);
    fmt::print("  de({}) = {:.2f} + ", Fv->shortSymbol(),
               deD * PhysConst::Hartree_invcm);
  }

  const auto Sx = m_Gold->Sigma_exchange(kappa, ev, vfk);

  if (Fv) {
    const auto deX = (*Fv) * (Sx * *Fv);
    fmt::print("{:.2f} = {:.2f}\n", deX * PhysConst::Hartree_invcm,
               (deD + deX) * PhysConst::Hartree_invcm);
  }

  return Sd + Sx;

  //
}

//==============================================================================
std::vector<double> NewSigma::calculate_fk(double ev,
                                           const DiracSpinor &v) const {

  assert(m_Fy0 && m_FyX);
  std::vector<double> vfk;
  for (int k = 0; k <= 8; ++k) {
    const auto Sd0 = m_Fy0->Sigma_direct(v.kappa(), ev, k);
    const auto SdX = m_FyX->Sigma_direct(v.kappa(), ev, k);
    // include hp when calculate fk: (ensure m_FyH exists!)
    // const auto Sd0 = m_FyH->Sigma_direct(v.kappa(), ev, k);
    // const auto SdX = m_Fy->Sigma_direct(v.kappa(), ev, k);
    const auto de0 = v * (Sd0 * v);
    const auto deX = v * (SdX * v);
    const auto fk = de0 != 0.0 ? deX / de0 : 1.0;
    vfk.push_back(fk);
    if (std::abs(fk - 1.0) < 0.01 || de0 == 0.0)
      break;
  }
  return vfk;
}

//==============================================================================
std::vector<double> NewSigma::calculate_etak(double ev,
                                             const DiracSpinor &v) const {
  assert(m_Fy0 && m_FyH);
  std::vector<double> vetak;
  for (int k = 0; k <= 6; ++k) {
    // slightly ineficient.. does Sd0 twice...
    const auto Sd0 = m_Fy0->Sigma_direct(v.kappa(), ev, k);
    const auto SdX = m_FyH->Sigma_direct(v.kappa(), ev, k);
    // Include screening when calc eta:
    // const auto Sd0 = m_FyX->Sigma_direct(v.kappa(), ev, k);
    // const auto SdX = m_Fy->Sigma_direct(v.kappa(), ev, k);
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
GMatrix NewSigma::formSigma_G(int kappa, double ev, const DiracSpinor *Fv) {

  if (!m_Gold) {
    m_Gold = Goldstone(m_basis, m_HF->core(), m_i0, m_stride, m_size,
                       m_n_min_core, m_includeG);
  }

  const auto Sd = m_Gold->Sigma_direct(kappa, ev, m_fk, m_etak);

  double deD{0.0};
  if (Fv) {
    deD = (*Fv) * (Sd * *Fv);
    fmt::print("  de({}) = {:.2f} + ", Fv->shortSymbol(),
               deD * PhysConst::Hartree_invcm);
  }

  const auto Sx = m_Gold->Sigma_exchange(kappa, ev, m_fk);

  if (Fv) {
    const auto deX = (*Fv) * (Sx * *Fv);
    fmt::print("{:.2f} = {:.2f}\n", deX * PhysConst::Hartree_invcm,
               (deD + deX) * PhysConst::Hartree_invcm);
  }
  return Sd + Sx;
}

//==============================================================================
void NewSigma::setup_Feynman() {

  if (!m_Gold) {
    // Also need Goldstone for Feynman
    m_Gold = Goldstone(m_basis, m_HF->core(), m_i0, m_stride, m_size,
                       m_n_min_core, m_includeG);
  }

  if (!m_Fy) {
    m_Fy = Feynman(m_HF, m_i0, m_stride, m_size, m_Foptions, m_n_min_core);
  }

  if (m_calculate_fk && !m_Fy0) {

    auto t_n_min = m_Fy->n_min();

    const auto t_stride = std::size_t(1.0 * double(m_stride)); // ?
    const auto t_size = (m_HF->grid().getIndex(m_rmax) - m_i0) / t_stride + 1;

    auto t_Foptions0{m_Foptions};
    auto t_FoptionsX{m_Foptions};
    auto t_FoptionsH{m_Foptions};

    t_Foptions0.screening = Screening::exclude;
    t_Foptions0.hole_particle = HoleParticle::exclude;

    t_FoptionsX.screening = Screening::include;
    t_FoptionsX.hole_particle = HoleParticle::exclude;

    t_FoptionsH.screening = Screening::exclude;
    t_FoptionsH.hole_particle = HoleParticle::include;

    if (m_Fy->screening()) {
      m_Fy0 =
          Feynman(m_HF, m_i0, t_stride, t_size, t_Foptions0, t_n_min, false);
      m_FyX =
          Feynman(m_HF, m_i0, t_stride, t_size, t_FoptionsX, t_n_min, false);
    }
    if (m_includeG && m_Fy->hole_particle()) { // only need eta for 'g'
      m_FyH =
          Feynman(m_HF, m_i0, t_stride, t_size, t_FoptionsH, t_n_min, false);
    }
  }
}

//==============================================================================
const SigmaData *NewSigma::get(int kappa, int n) const {
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
    // If not found, look (recursively) for closest n
    // if (it == m_Sigmas.cend()) {
    //   std::cout << "\n***" << it->kappa << " " << it->n << " " << it->en
    //             << "\n";
    // }
    return it != m_Sigmas.cend() ? &(*it) : get(kappa, n - 1);
  }
}

//==============================================================================
const GMatrix *NewSigma::getSigma(int kappa, int n) const {
  const auto Sig = get(kappa, n);
  return Sig ? &(Sig->Sigma) : nullptr;
}

//==============================================================================
double NewSigma::getLambda(int kappa, int n) const {
  const auto Sig = get(kappa, n);
  return Sig ? (Sig->lambda) : 1.0;
}

//==============================================================================
//! returns Spinor: Sigma|Fv>
//! @details If Sigma for kappa_v doesn't exist, returns |0>.
DiracSpinor NewSigma::SigmaFv(const DiracSpinor &Fv) const {
  auto Sv = get(Fv.kappa(), Fv.n());
  return Sv ? Sv->lambda * (Sv->Sigma * Fv) : 0.0 * Fv;
}

//==============================================================================
//! Stores scaling factors, lambda, for each kappa (Sigma -> lamda*Sigma)
void NewSigma::scale_Sigma(const std::vector<double> &lambdas) {
  for (std::size_t i = 0; i < m_Sigmas.size() && i < lambdas.size(); ++i) {
    m_Sigmas.at(i).lambda = lambdas.at(i);
  }
}

//==============================================================================
// if n=0, scales _all_
void NewSigma::scale_Sigma(double lambda, int kappa, int n) {
  for (auto &Sig : m_Sigmas) {
    if (Sig.kappa == kappa && (Sig.n == n || n <= 0)) {
      Sig.lambda = lambda;
    }
  }
}

//==============================================================================
//! Prints the scaling factors to screen
void NewSigma::print_scaling() const {
  std::cout << "Scaling factors: lambda = ";
  for (const auto &Sig : m_Sigmas) {
    std::cout << Sig.lambda << ", ";
  }
  std::cout << "\n";
}

//==============================================================================
//! Prints the sub-grid parameters to screen
void NewSigma::print_subGrid() const {
  printf("Sigma sub-grid: r=(%.1e, %.1f)aB with %i points. [i0=%i, "
         "stride=%i]\n",
         m_HF->grid().r(m_i0), m_HF->grid().r(m_i0 + m_stride * (m_size - 1)),
         int(m_size), int(m_i0), int(m_stride));
}

//==============================================================================
// void NewSigma::print_info() const {
//   print_subGrid();
//   if (!m_Sigmas.empty())
//     std::cout << "Have Sigma for:\n";
//   for (const auto &Sig : m_Sigmas) {
//     std::cout << Sig.n << " " << AtomData::kappa_symbol(Sig.kappa)
//               << " en=" << Sig.en;
//     if (Sig.lambda != 1.0) {
//       std::cout << " (lambda=" << Sig.lambda << ")";
//     }
//     std::cout << "\n";
//   }
// }

//==============================================================================
bool NewSigma::read_write(const std::string &fname, IO::FRW::RoW rw) {

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

  // XXXX Add method and basis data, fk etc.
  // // read/write basis config:
  // std::string basis_config =
  //     (rw == IO::FRW::write) ? DiracSpinor::state_config(m_excited) : "";
  // rw_binary(iofs, rw, basis_config);
  // // XXX Add info on method! (or, different filename!)

  // Sub-grid:
  rw_binary(iofs, rw, m_r0, m_rmax, m_stride, m_i0, m_size, m_includeG);

  // Number of kappas (number of Sigma/G matrices)
  std::size_t num_Sigmas = rw == IO::FRW::write ? m_Sigmas.size() : 0;
  rw_binary(iofs, rw, num_Sigmas);

  for (std::size_t iS = 0; iS < num_Sigmas; ++iS) {

    if (rw == IO::FRW::read) {
      m_Sigmas.push_back(
          {0, 0.0,
           GMatrix{m_i0, m_stride, m_size, m_includeG, m_HF->grid_sptr()}, 0,
           1.0}); // don't read/write lamba
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
    for (const auto &Sig : m_Sigmas) {
      fmt::print("kappa = {}, ev = {:.5f}", Sig.kappa, Sig.en);
      if (Sig.n > 0) {
        fmt::print(" (n = {})", Sig.n);
      }
      std::cout << "\n";
    }
  }
  return true;
}

} // namespace MBPT