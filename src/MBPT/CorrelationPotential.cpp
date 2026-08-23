#include "CorrelationPotential.hpp"
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "MBPT/Ladder.hpp"
#include "MBPT/Sigma2.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"
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
  const std::vector<double> &etak, const std::string &ladder_file,
  bool form_derivative, bool fk_both_lines, bool feynman_exchange)
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
    m_fk_both_lines(fk_both_lines),
    m_feynman_exchange(feynman_exchange),
    m_fname(fname),
    m_ladder_file(ladder_file),
    m_form_derivative(form_derivative) {

  std::cout << "\nConstruct Correlation Potential\n";

  // attempt to read in Sigma file:
  // (Just contains Sigma matrix, nothing else)
  const bool read_ok = read_write(fname, IO::FRW::read);

  // Read (separate) ladder file, Sigma_L, if given. Produced by the Ladder{}
  // block; stored separately from the base Sigma (independent of read_ok)
  if (!m_ladder_file.empty()) {
    auto SLs = read_SigmaL(m_ladder_file, m_HF->grid_sptr());
    for (auto &sl : SLs) {
      fmt::print("Sigma_L: kappa = {:>2}, en = {:+.5f} (n = {})", sl.kappa,
                 sl.en, sl.n);
      // de = <v|Sigma_L|v> (basis version of state) - visual check for users
      const auto pFv =
        std::find_if(m_basis.cbegin(), m_basis.cend(), [&sl](const auto &F) {
          return F.n() == sl.n && F.kappa() == sl.kappa;
        });
      if (pFv != m_basis.cend()) {
        const auto de = *pFv * (sl.SL * *pFv);
        fmt::print(", de = {:+.5e}", de);
      }
      std::cout << "\n";
      m_Sigma_L.push_back({sl.kappa, sl.en, std::move(sl.SL), sl.n, 1.0});
    }
    if (m_Sigma_L.empty()) {
      std::cout << "WARNING: no ladder Sigma_L read from: " << m_ladder_file
                << " - ladder will not be included\n";
    }
  }

  if (!read_ok) {

    if (m_method == SigmaMethod::Feynman && m_feynman_exchange) {
      std::cout << "Using Feynman method for direct and exchange diagrams\n";
    } else if (m_method == SigmaMethod::Feynman) {
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
      if (m_fk_both_lines) {
        std::cout << "fk applied to both Coulomb lines in exchange\n";
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
        if (m_fk_both_lines) {
          std::cout << "fk applied to both Coulomb lines in exchange\n";
        }
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
    std::find_if(m_Sigmas.begin(), m_Sigmas.end(), [kappa, n](const auto &s) {
      return s.kappa == kappa && (s.n == n || n <= 0);
    });
  if (it != m_Sigmas.end()) {
    // have sigma already!
    // print deets!
    auto de = Fv ? *Fv * (it->Sigma * *Fv) : 0.0;
    fmt::print("Have Sigma: kappa = {:>2}, en = {:+.5f}, de = {:+.5e}\n",
               it->kappa, it->en, de);
    // May still need the derivative (e.g., Sigma read from an older file)
    if (m_form_derivative && !get_derivative(kappa, it->n)) {
      const auto fk_state = state_fk(it->en, Fv);
      // Backfill stored fk (e.g., Sigma read from an older file)
      if (fk_state && it->fk.empty()) {
        it->fk = *fk_state;
      }
      form_derivative(kappa, it->en, it->n, Fv, it->Sigma,
                      fk_state ? &*fk_state : nullptr);
    }
    return;
  }
  if (Fv) {
    assert(Fv->kappa() == kappa);
    fmt::print("Form Σ for {} at e = {:.4f} au = {:.2f} /cm\n",
               Fv->shortSymbol(), ev, ev * PhysConst::Hartree_invcm);
  } else {
    fmt::print("Form Σ for kappa={} at e = {:.4f}\n", kappa, ev);
  }

  // fk screening factors for this state: calculated once, here; the base
  // Sigma and its derivative use the same values
  const auto fk_state = state_fk(ev, Fv);
  const auto fk_pointer = fk_state ? &*fk_state : nullptr;

  auto S = m_method == SigmaMethod::Feynman ?
             formSigma_F(kappa, ev, Fv, fk_pointer) :
             formSigma_G(kappa, ev, Fv);

  m_Sigmas.push_back({kappa, ev, std::move(S), n, 1.0});

  // Record the fk actually used for this Sigma (calculated, or manual m_fk):
  // stored with Sigma for re-use (e.g., Sigma_2 screening in CI)
  if (fk_state) {
    m_Sigmas.back().fk = *fk_state;
  } else if (!m_fk.empty()) {
    m_Sigmas.back().fk = m_fk;
  }

  if (m_form_derivative && !get_derivative(kappa, n)) {
    form_derivative(kappa, ev, n, Fv, m_Sigmas.back().Sigma, fk_pointer);
  }
}

//==============================================================================
std::optional<std::vector<double>>
CorrelationPotential::state_fk(double ev, const DiracSpinor *Fv) {
  const bool screening = m_Foptions.screening == Screening::include;
  if (m_method != SigmaMethod::Feynman || !m_calculate_fk || !screening ||
      Fv == nullptr) {
    return std::nullopt;
  }
  if (!m_Fy) {
    setup_Feynman();
  }
  auto fk = calculate_fk(ev, *Fv);
  std::cout << "  fk   = {";
  for (const auto &each_fk : fk) {
    printf("%.3f, ", each_fk);
  }
  std::cout << "}\n";
  // If not stored, store first screening factors
  if (m_fk.empty()) {
    m_fk = fk;
  }
  return fk;
}

//==============================================================================
void CorrelationPotential::form_derivative(
  int kappa, double ev, int n, const DiracSpinor *Fv, const GMatrix &Sigma0,
  const std::vector<double> *given_fk) {
  fmt::print("  Form dSigma/dE for kappa={} at e = {:.4f}\n", kappa, ev);
  std::cout << std::flush;
  // One extra Sigma evaluation, re-using the base Sigma and its fk.
  // Print nothing for it: only the Sigma we actually use is reported
  auto dS = m_method == SigmaMethod::Feynman ?
              formSigma_F(kappa, ev + m_delta_en, Fv, given_fk, false) :
              formSigma_G(kappa, ev + m_delta_en, Fv, false);
  dS -= Sigma0;
  dS *= (1.0 / m_delta_en);
  // if (Fv) {
  //   const auto dde = *Fv * (dS * *Fv);
  //   fmt::print("  d(de)/dE({}) = {:+.5f}\n", Fv->shortSymbol(),
  //              dde * PhysConst::Hartree_invcm);
  // }
  m_dSigma.push_back({kappa, ev, std::move(dS), n, 1.0});
}

//==============================================================================
GMatrix CorrelationPotential::formSigma_F(int kappa, double ev,
                                          const DiracSpinor *Fv,
                                          const std::vector<double> *given_fk,
                                          bool print) {

  if (!m_Fy) {
    setup_Feynman();
  }

  std::vector<double> vfk;

  if (given_fk != nullptr) {
    // Screening factors already calculated for this state (see state_fk)
    vfk = *given_fk;
  } else if (m_calculate_fk && m_Fy->screening()) {
    assert(Fv != nullptr && "Cannot calculate fk without Fv");
    vfk = calculate_fk(ev, *Fv);
    if (print) {
      std::cout << "  fk   = {";
      for (auto &e : vfk) {
        printf("%.3f, ", e);
      }
      std::cout << "}\n";
    }
    // If not stored, store first screening factors
    if (m_fk.empty()) {
      m_fk = vfk;
    }
  } else {
    vfk = m_fk;
  }

  if (Fv && print) {
    fmt::print("  de({}) = ", Fv->shortSymbol());
    std::cout << std::flush;
  }

  auto Sd = m_Fy->Sigma_direct(kappa, ev);

  double deD{0.0};
  if (Fv && print) {
    deD = (*Fv) * (Sd * *Fv);
    fmt::print("{:.2f} + ", deD * PhysConst::Hartree_invcm);
    std::cout << std::flush;
  }

  const auto Sx = m_feynman_exchange ?
                    m_Fy->Sigma_exchange(kappa, ev) :
                    m_Gold->Sigma_exchange(kappa, ev, vfk, m_fk_both_lines);

  if (Fv && print) {
    const auto deX = (*Fv) * (Sx * *Fv);
    fmt::print("{:.2f} = {:.2f}\n", deX * PhysConst::Hartree_invcm,
               (deD + deX) * PhysConst::Hartree_invcm);
    std::cout << std::flush;
  }

  if (m_includeBreit_b2) {
    // nb: do some extra work to calculate it seperately (Qk and Pk)..
    // But, since selection rules are different, it's better this way
    if (Fv && print) {
      fmt::print("  de[B2]  = ");
      std::cout << std::flush;
    }
    const auto dS =
      m_Gold->dSigma_Breit2(kappa, ev, m_fk, m_etak, 99, m_n_max_breit);
    if (Fv && print) {
      const auto deB2 = (*Fv) * (dS * *Fv);
      fmt::print("{:.2f}\n", deB2 * PhysConst::Hartree_invcm);
      std::cout << std::flush;
    }
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

  // Sigma_d for every k in one pass (Green's fns shared across k):
  // much faster than separate per-k Sigma_direct calls
  const auto Sd0_k = m_Fy0->Sigma_direct_each_k(v.kappa(), ev);
  const auto SdX_k = m_FyX->Sigma_direct_each_k(v.kappa(), ev);

  std::vector<double> vfk;
  for (auto k = 0ul; k <= 9 && k < Sd0_k.size(); ++k) {
    const auto de0 = v * (Sd0_k[k] * v);
    const auto deX = v * (SdX_k[k] * v);
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
  // Sigma_d for every k in one pass (see calculate_fk)
  const auto Sd0_k = m_FyX->Sigma_direct_each_k(v.kappa(), ev);
  const auto SdH_k = m_Fy->Sigma_direct_each_k(v.kappa(), ev);
  std::vector<double> vetak;
  for (auto k = 0ul; k <= 6 && k < Sd0_k.size(); ++k) {
    // Include screening when calc eta:
    const auto de0 = v * (Sd0_k[k] * v);
    const auto deH = v * (SdH_k[k] * v);
    const auto etak = de0 != 0.0 ? deH / de0 : 1.0;
    vetak.push_back(etak);
    if (std::abs(etak - 1.0) < 0.01 || de0 == 0.0)
      break;
  }
  return vetak;
}

//==============================================================================
GMatrix CorrelationPotential::formSigma_G(int kappa, double ev,
                                          const DiracSpinor *Fv, bool print) {

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
              m_Gold->Sigma_both(kappa, ev, m_fk, m_etak, 99, m_fk_both_lines);

  double deD{0.0};
  if (Fv && print) {
    deD = (*Fv) * (Sd * *Fv);
    fmt::print("  de({}) = {:.2f} ", Fv->shortSymbol(),
               deD * PhysConst::Hartree_invcm);
  }

  if (exchange_seperately) {
    const auto Sx = m_Gold->Sigma_exchange(kappa, ev, m_fk, m_fk_both_lines);
    if (Fv && print) {
      const auto deX = (*Fv) * (Sx * *Fv);
      fmt::print("+ {:.2f} = {:.2f}\n", deX * PhysConst::Hartree_invcm,
                 (deD + deX) * PhysConst::Hartree_invcm);
    }
    Sd += Sx;
  } else if (print) {
    std::cout << "\n" << std::flush;
  }

  if (m_includeBreit_b2) {
    // nb: do some extra work to calculate it seperately (Qk and Pk)..
    // But, since selection rules are different, it's better this way
    if (Fv && print) {
      fmt::print("  de[B2]  = ");
      std::cout << std::flush;
    }
    const auto dS =
      m_Gold->dSigma_Breit2(kappa, ev, m_fk, m_etak, 99, m_n_max_breit);
    if (Fv && print) {
      const auto deB2 = (*Fv) * (dS * *Fv);
      fmt::print("{:.2f}\n", deB2 * PhysConst::Hartree_invcm);
      std::cout << std::flush;
    }
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

  if (m_calculate_fk && !m_Fy0 && m_Fy->screening()) {

    // Fy with no screening (+no hp)
    auto t_Foptions0{m_Foptions};
    t_Foptions0.screening = Screening::exclude;
    t_Foptions0.hole_particle = HoleParticle::exclude;
    m_Fy0 = Feynman(m_HF, m_i0, m_stride, m_size, t_Foptions0, m_n_min_core,
                    m_includeG, false, m_fname, false);

    // Fy with screening (but no hp)
    auto t_FoptionsX{m_Foptions};
    t_FoptionsX.screening = Screening::include;
    t_FoptionsX.hole_particle = HoleParticle::exclude;
    m_FyX = m_Fy->hole_particle() ?
              Feynman(m_HF, m_i0, m_stride, m_size, t_FoptionsX, m_n_min_core,
                      m_includeG, false, m_fname, false) :
              m_Fy;

    // Fy0 and FyX differ only in screening, so the polarisation loop (the
    // expensive step) is formed once and handed from one to the other
    std::optional<std::vector<std::vector<ComplexRMatrix>>> pi_wk{};
    for (auto *Fy : {&*m_Fy0, &*m_FyX}) {
      if (Fy->has_qpiq() || Fy->read_qpiq(m_fname))
        continue;
      if (!pi_wk) {
        pi_wk = Fy->polarisation_wk();
      }
      Fy->form_qpiq(*pi_wk);
      Fy->write_qpiq(m_fname);
    }
  }
}

//==============================================================================
namespace {
// Finds Sigma in list for given kappa (and n); shared by get()/get_ladder()
const SigmaData *find_Sigma(const std::vector<SigmaData> &Sigmas, int kappa,
                            int n) {
  if (n <= 0) {
    // returns FIRST sigma that has correct kappa, order matters!
    const auto it =
      std::find_if(Sigmas.cbegin(), Sigmas.cend(),
                   [kappa](const auto &s) { return s.kappa == kappa; });
    return it != Sigmas.cend() ? &(*it) : nullptr;
  } else {
    // Find first Sigma that matches kappa _and_ n
    const auto it =
      std::find_if(Sigmas.cbegin(), Sigmas.cend(), [kappa, n](const auto &s) {
        return s.kappa == kappa && s.n == n;
      });
    // If not found, look (recursively) for next lowest n
    return it != Sigmas.cend() ? &(*it) : find_Sigma(Sigmas, kappa, n - 1);
  }
}
} // namespace

const SigmaData *CorrelationPotential::get(int kappa, int n) const {
  return find_Sigma(m_Sigmas, kappa, n);
}

const SigmaData *CorrelationPotential::get_ladder(int kappa, int n) const {
  return find_Sigma(m_Sigma_L, kappa, n);
}

const SigmaData *CorrelationPotential::get_derivative(int kappa, int n) const {
  return find_Sigma(m_dSigma, kappa, n);
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
DiracSpinor CorrelationPotential::SigmaFv(const DiracSpinor &Fv) const {
  const auto Sv = get(Fv.kappa(), Fv.n());
  const auto Sl = get_ladder(Fv.kappa(), Fv.n());
  if (!Sv && !Sl)
    return 0.0 * Fv;
  auto SF = Sv ? Sv->Sigma * Fv : 0.0 * Fv;
  if (Sl) {
    SF += Sl->Sigma * Fv;
  }
  // lambda (from base Sigma) scales both the base and ladder parts
  const auto lambda = Sv ? Sv->lambda : 1.0;
  return lambda * SF;
}

//==============================================================================
void CorrelationPotential::print_de(const std::vector<DiracSpinor> &valence) {

  if (valence.empty() || m_basis.empty()) {
    return;
  }

  if (!m_Gold) {
    m_Gold =
      Goldstone(m_basis, m_HF->core(), m_i0, m_stride, m_size, m_n_min_core,
                m_includeG, m_includeBreit_b2 ? m_HF->vBreit() : nullptr);
  }
  const auto &[core, excited] = m_Gold->basis();

  std::cout << "\nMBPT(2) (/cm):\n";
  fmt::print("{:5s} {:>11s} {:>11s} {:>11s} {:>12s} {:>9s}\n", "state",
             "direct", "exchange", "total", "<v|Sigma|v>", "eps");
  for (const auto &v : valence) {
    const auto [de_direct, de_exchange] =
      MBPT::Sigma_vw_direct_exchange(v, v, m_Gold->Yeh(), core, excited);
    const auto de_2 = de_direct + de_exchange;
    const auto de_Sigma = v * SigmaFv(v);
    const auto eps = de_Sigma / de_2 - 1.0;
    fmt::print("{:5s} {:>11.2f} {:>11.2f} {:>11.2f} {:>12.2f} {:>9.1e}\n",
               v.shortSymbol(), de_direct * PhysConst::Hartree_invcm,
               de_exchange * PhysConst::Hartree_invcm,
               de_2 * PhysConst::Hartree_invcm,
               de_Sigma * PhysConst::Hartree_invcm, eps);
  }
  std::cout << std::flush;
}

//==============================================================================
DiracSpinor CorrelationPotential::dSigmaFv(const DiracSpinor &Fv) const {
  const auto dSv = get_derivative(Fv.kappa(), Fv.n());
  if (!dSv) {
    return 0.0 * Fv;
  }
  // same lambda as the base Sigma (nb: ladder has no derivative)
  const auto lambda = getLambda(Fv.kappa(), Fv.n());
  return lambda * (dSv->Sigma * Fv);
}

//==============================================================================
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

  // dSigma/dE matrices. This block was appended to the format later: older
  // files simply end here, and hold no derivatives
  const bool have_deriv_block =
    rw == IO::FRW::write || iofs.peek() != std::fstream::traits_type::eof();
  if (have_deriv_block) {
    std::size_t num_dSigma = rw == IO::FRW::write ? m_dSigma.size() : 0;
    rw_binary(iofs, rw, num_dSigma);
    for (std::size_t iS = 0; iS < num_dSigma; ++iS) {
      if (rw == IO::FRW::read) {
        m_dSigma.push_back(
          {0, 0.0,
           GMatrix{m_i0, m_stride, m_size, m_includeG, m_HF->grid_sptr()}, 0,
           1.0});
      }
      auto &Sig = m_dSigma.at(iS);
      rw_binary(iofs, rw, Sig.kappa, Sig.en, Sig.n);
      auto &Gk = Sig.Sigma;
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
  }

  // Per-Sigma fk screening factors, in m_Sigmas order (appended to the
  // format later still: older files end before this block)
  const bool have_fk_block =
    rw == IO::FRW::write || iofs.peek() != std::fstream::traits_type::eof();
  if (have_fk_block) {
    std::size_t num_fk = rw == IO::FRW::write ? m_Sigmas.size() : 0;
    rw_binary(iofs, rw, num_fk);
    for (std::size_t iS = 0; iS < num_fk && iS < m_Sigmas.size(); ++iS) {
      auto &Sig = m_Sigmas.at(iS);
      int kappa = Sig.kappa;
      rw_binary(iofs, rw, kappa, Sig.fk);
      assert(kappa == Sig.kappa && "fk block out of sync with Sigmas");
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
std::string CorrelationPotential::method_string() const {
  std::string out = m_method == SigmaMethod::Feynman ? "Feynman" : "Goldstone";
  if (m_method == SigmaMethod::Feynman) {
    const auto scr = m_Foptions.screening == Screening::include;
    const auto hp = m_Foptions.hole_particle == HoleParticle::include;
    if (scr && hp) {
      out += ", all-order";
    } else {
      if (scr) {
        out += "+scr";
      }
      if (m_Foptions.hole_particle != HoleParticle::exclude) {
        out += hp ? "+hp" : "+hp0";
      }
    }
  }
  if (!m_ladder_file.empty()) {
    out += ", ladder";
  }
  return out;
}

//==============================================================================
std::vector<double> CorrelationPotential::average_fk(int l_max) const {
  // Weighted by l, not kappa: average each l's fine-structure pair first,
  // then take the mean over the ls (s counts the same as p, d, ...)
  std::vector<std::vector<double>> per_l;
  for (int l = 0; l <= l_max; ++l) {

    std::vector<const std::vector<double> *> stored;
    for (const auto kappa : {l, -(l + 1)}) {
      if (kappa == 0)
        continue;
      const auto Sigma_kappa = get(kappa);
      if (Sigma_kappa != nullptr && !Sigma_kappa->fk.empty()) {
        stored.push_back(&Sigma_kappa->fk);
      }
    }
    if (stored.empty())
      continue;

    auto size = stored.front()->size();
    for (const auto &each : stored) {
      size = std::min(size, each->size());
    }
    std::vector<double> fk_l(size, 0.0);
    for (const auto &each : stored) {
      for (std::size_t k = 0; k < size; ++k) {
        fk_l[k] += (*each)[k] / double(stored.size());
      }
    }
    per_l.push_back(std::move(fk_l));
  }

  if (per_l.empty()) {
    return {};
  }
  auto size = per_l.front().size();
  for (const auto &each : per_l) {
    size = std::min(size, each.size());
  }
  std::vector<double> average(size, 0.0);
  for (const auto &each : per_l) {
    for (std::size_t k = 0; k < size; ++k) {
      average[k] += each[k] / double(per_l.size());
    }
  }
  return average;
}

//==============================================================================
std::string CorrelationPotential::lambda_string() const {
  // Rounded values, so stable across runs
  std::string lambdas;
  for (const auto &Sig : m_Sigmas) {
    if (std::abs(Sig.lambda - 1.0) > 1.0e-8) {
      lambdas += fmt::format("{}={:.4f},", Sig.kappa, Sig.lambda);
    }
  }
  return lambdas;
}

//==============================================================================
void CorrelationPotential::print_info() const {
  for (const auto &Sig : m_Sigmas) {
    fmt::print("kappa = {:>2}, ev = {:+.5f}", Sig.kappa, Sig.en);
    if (std::abs(Sig.lambda - 1.0) > 1.0e-8) {
      fmt::print(" : scaled with λ = {:.5f}", Sig.lambda);
    }
    if (get_derivative(Sig.kappa, Sig.n)) {
      fmt::print(" (with dSigma/dE)");
    }
    std::cout << "\n";
  }
}

} // namespace MBPT