#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracODE/include.hpp"
#include "HF/HartreeFock.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <vector>
//==============================================================================
ContinuumOrbitals::ContinuumOrbitals(const Wavefunction &wf)
  : p_rgrid(wf.grid_sptr()), p_hf(wf.vHF()), m_alpha(wf.alpha()) {}

ContinuumOrbitals::ContinuumOrbitals(const HF::HartreeFock *hf)
  : p_rgrid(hf->grid_sptr()), p_hf(hf), m_alpha(hf->alpha()) {}

//==============================================================================
double ContinuumOrbitals::check_orthog(bool print) const {
  double worst = 0.0;
  if (p_hf == nullptr)
    return worst;
  for (const auto &Fc : orbitals) {
    for (const auto &Fn : p_hf->core()) {
      if (Fn.kappa() != Fc.kappa())
        continue;
      const auto eps = Fc * Fn;
      if (std::abs(eps) > std::abs(worst))
        worst = eps;
      if (print) {
        std::cout << "<" << Fc.shortSymbol() << "|" << Fn.shortSymbol()
                  << "> = ";
        printf("%.1e\n", eps);
      }
    }
  }
  return worst;
}

//==============================================================================
int ContinuumOrbitals::solveContinuumHF(double ec, int min_l, int max_l,
                                        const DiracSpinor *Fi,
                                        bool force_rescale, bool subtract_self,
                                        bool force_orthog_Fi) {

  // include Hartree here? Probably shouldn't, since we do "core Hartree"
  const auto self_consistant = (p_hf->method() == HF::Method::HartreeFock ||
                                p_hf->method() == HF::Method::Hartree ||
                                p_hf->method() == HF::Method::ApproxHF);

  // Also orthogonalise against entire core: (makes no difference)
  const bool orthog_core = force_orthog_Fi;

  using namespace qip::overloads;
  auto vc = p_hf->vlocal();
  // vdir_0 = y^0_{Fi,Fi}: direct self-interaction of ionised orbital.
  // Captured here so it can also be passed to IncludeExchange for the
  // Hermitian projection correction (1-Pc)V_0(1-Pc).
  std::vector<double> vdir_0;
  if ((Fi != nullptr) && subtract_self && self_consistant) {
    // Subtract off the self-interaction direct part: V-self(r) = y^0_ii(r)
    vdir_0 = Coulomb::yk_ab(0, *Fi, *Fi);
    vc -= vdir_0;
  }

  // We may wish to do this to test things, but not for final calculations:
  // if (force_rescale && subtract_self) {
  //   fmt2::styled_print(fg(fmt::color::orange), "\nWARNING 65: ");
  //   fmt::print("Should not subtract self interaction _and_ "
  //              "rescale V(r): do one or the other\n");
  // }

  //Optionally, re-scale large-r part of local potential,
  // so goes like -1/r large r
  // Note: doesn't inclue exchange..
  // This also kills orthogonality for HF...
  if (force_rescale && self_consistant) {
    const auto Zion = std::max(1.0, -vc.back() * p_rgrid->r().back());
    // nb: this agrees best with Dzuba, but bad for orthog
    for (std::size_t i = 0; i < p_rgrid->num_points(); ++i) {
      if (vc[i] > -Zion / p_rgrid->r(i)) {
        vc[i] = -Zion / p_rgrid->r(i);
      }
    }
  }

  // loop through each kappa state
  for (std::size_t k_i = 0; true; ++k_i) {
    const auto kappa = Angular::kindex_to_kappa(k_i);
    const auto l = Angular::l_k(kappa);
    if (l < min_l)
      continue;
    if (l > max_l)
      break;

    auto &Fc = orbitals.emplace_back(0, kappa, p_rgrid);
    Fc.en() = ec;
    // solve initial, without exchange term
    DiracODE::solveContinuum(Fc, ec, vc, m_alpha);
    // Then, include exchange correction:
    if (p_hf != nullptr && !p_hf->is_localQ()) {
      IncludeExchange(&Fc, Fi, vc, vdir_0);
    }
  }

  // Orthogonalise against entire core:
  if (orthog_core) {
    for (auto &Fc : orbitals) {
      for (const auto &Fa : p_hf->core()) {
        Fc.orthog(Fa);
        // if (Fa.kappa() == Fc.kappa())
        // Fc -= (Fc * Fa) * Fa;
      }
    }
  }

  // Force orthogonality to specific state
  if (force_orthog_Fi && Fi != nullptr) {
    for (auto &Fc : orbitals) {
      Fc.orthog(*Fi);
      // if (Fi->kappa() == Fc.kappa()) {
      //   Fc -= (*Fi * Fc) * *Fi;
      // }
    }
  }

  return 0;
}

//******************************************************************************
// V_hp|Fa> = (1/(2j_i+1)) vex(Fa,Fi) + y^0_{ii}·Fa
DiracSpinor ContinuumOrbitals::hp_apply_V0(const DiracSpinor &Fa,
                                           const DiracSpinor *Fi,
                                           const std::vector<double> &vdir_0) {
  auto Vhp =
    (1.0 / Fi->twojp1()) * HF::vexFa(Fa, std::vector<DiracSpinor>{*Fi});
  if (!vdir_0.empty())
    Vhp += vdir_0 * Fa;
  return Vhp;
}

//******************************************************************************
// V_hp_core[a] = V_hp|a>,  Vhp_ba[b,a] = <b|V_hp|a>  for same-kappa core pairs.
void ContinuumOrbitals::hp_precompute(const std::vector<DiracSpinor> &core,
                                      int kappa, const DiracSpinor *Fi,
                                      const std::vector<double> &vdir_0,
                                      std::vector<DiracSpinor> *V_hp_core,
                                      std::vector<double> *Vhp_ba) {
  for (const auto &Fa : core) {
    if (Fa.kappa() == kappa) {
      V_hp_core->push_back(hp_apply_V0(Fa, Fi, vdir_0));
    }
  }
  const auto NK = V_hp_core->size();
  Vhp_ba->assign(NK * NK, 0.0);
  // ib_kc: index into the same-kappa subspace (matches V_hp_core ordering)
  std::size_t ib_kc = 0;
  for (const auto &Fb : core) {
    if (Fb.kappa() != kappa) {
      continue;
    }
    for (std::size_t ia = 0; ia < NK; ++ia) {
      (*Vhp_ba)[ib_kc * NK + ia] = Fb * (*V_hp_core)[ia];
    }
    ++ib_kc;
  }
}

//******************************************************************************
// Returns Pc V_hp|Fc> + V_hp Pc|Fc> - Pc V_hp Pc|Fc>
// where Pc projects onto same-kappa core states.
DiracSpinor ContinuumOrbitals::add_hp_projection(
  const DiracSpinor &Fc, const DiracSpinor *Fi,
  const std::vector<double> &vdir_0, const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &V_hp_core,
  const std::vector<double> &Vhp_ba) {

  const auto NK = V_hp_core.size();
  if (NK == 0)
    return 0.0 * Fc;

  const auto V_hp_Fc = hp_apply_V0(Fc, Fi, vdir_0);
  auto proj = 0.0 * V_hp_Fc;

  // ia_kc, ib_kc: indices into the same-kappa subspace (match V_hp_core / M ordering)
  std::size_t ia_kc = 0;
  for (const auto &Fa : core) {
    if (Fa.kappa() != Fc.kappa()) {
      continue;
    }
    // <a|c>
    const double ov_ac = Fa * Fc;
    // <a|V_hp|c>
    const double ov_aVc = Fa * V_hp_Fc;
    // Pc V_hp |c>
    proj += ov_aVc * Fa;
    // V_hp Pc |c>
    proj += ov_ac * V_hp_core[ia_kc];
    // ib_kc: index for core states that have kappa=kappa_c only
    std::size_t ib_kc = 0;
    for (const auto &Fb : core) {
      if (Fb.kappa() != Fc.kappa()) {
        continue;
      }
      // -Pc V_hp Pc |c>
      proj -= (ov_ac * Vhp_ba[ib_kc * NK + ia_kc]) * Fb;
      ++ib_kc;
    }
    ++ia_kc;
  }
  return proj;
}

//******************************************************************************
void ContinuumOrbitals::IncludeExchange(DiracSpinor *Fc, const DiracSpinor *Fi,
                                        const std::vector<double> &vc,
                                        const std::vector<double> &vdir_0) {
  // Solves the HF exchange problem for the continuum state Fc, using the
  // Hermitian V^{N-1} potential:
  //   V = V^N - (1-Pc)V_hp(1-Pc)
  // where V_hp is the hole-particle potential (direct + exchange of the ionised
  // orbital Fi), and Pc projects onto the same-kappa core states.
  // The symmetric form guarantees that core eigenstates remain eigenstates of
  // H, so continuum states are naturally orthogonal to the core.

  const int max_its = 50;
  const double conv_target = 1.0e-6;

  const bool do_proj =
    (Fi != nullptr) && (p_hf->method() == HF::Method::HartreeFock);

  // Precompute V_hp|a> and M[b,a] = <b|V_hp|a> for same-kappa core pairs.
  // -V_hp already in: vc (direct) and vexFa subtract_one (exchange).
  // Each iteration adds: Pc V_hp|Fc> + V_hp Pc|Fc> - Pc V_hp Pc|Fc>
  // upgrading V^N - V_hp → V^N - (1-Pc) V_hp (1-Pc).
  std::vector<DiracSpinor> V_hp_core;
  std::vector<double> Vhp_ba;
  if (do_proj)
    hp_precompute(p_hf->core(), Fc->kappa(), Fi, vdir_0, &V_hp_core, &Vhp_ba);

  for (int it = 0; it <= max_its; ++it) {
    const auto vx0 = HF::vex_approx(*Fc, p_hf->core());
    const auto vl = qip::add(vc, vx0);
    const auto Fc0 = *Fc;

    if (p_hf->method() == HF::Method::HartreeFock) {
      auto VxFc = HF::vexFa(*Fc, p_hf->core(), 99, Fi) - vx0 * *Fc;
      if (do_proj)
        VxFc +=
          add_hp_projection(*Fc, Fi, vdir_0, p_hf->core(), V_hp_core, Vhp_ba);
      DiracODE::solveContinuum(*Fc, Fc->en(), vl, m_alpha, &VxFc, &Fc0);
    } else {
      DiracODE::solveContinuum(*Fc, Fc->en(), vl, m_alpha);
    }

    const auto delta = Fc0 - *Fc;
    const auto eps = (delta * delta) / (*Fc * *Fc);
    if (eps < conv_target || it == max_its)
      break;

    *Fc += Fc0;
    *Fc *= 0.5;
  }
}

//******************************************************************************
int ContinuumOrbitals::solveContinuumZeff(double ec, int min_l, int max_l,
                                          double Z_eff, const DiracSpinor *Fi,
                                          bool force_orthog) {
  // Solves Dirac equation for H-like potential (Zeff model)
  // Same Zeff as used by DarkARC (eqn B35 of arxiv:1912.08204):
  // Zeff = sqrt{I_{njl} eV / 13.6 eV} * n
  // au: Zeff = sqrt{2 * I_{njl}} * n

  // Also orthogonalise against entire core: (make no difference)
  const bool orthog_core = force_orthog;

  // Zeff potential (pointlike nucleus, spherical with Rn=0):
  const auto vc = Nuclear::sphericalNuclearPotential(Z_eff, 0.0, p_rgrid->r());

  // loop through each kappa state
  for (auto k_i = 0ul; true; ++k_i) {
    const auto kappa = Angular::kindex_to_kappa(k_i);
    const auto l = Angular::l_k(kappa);
    if (l < min_l)
      continue;
    if (l > max_l)
      break;

    auto &Fc = orbitals.emplace_back(0, kappa, p_rgrid);
    Fc.en() = ec;
    // solve initial, without exchange term
    DiracODE::solveContinuum(Fc, ec, vc, m_alpha);

  } // kappa

  // Orthogonalise against entire core:
  if (orthog_core) {
    for (auto &Fc : orbitals) {
      for (const auto &Fa : p_hf->core()) {
        if (Fa.kappa() == Fc.kappa())
          Fc -= (Fc * Fa) * Fa;
      }
    }
  }

  // Forcing orthogonality between continuum states and current core state
  // (have to do this _after_ orthog_core, since that slightly breaks this)
  if (force_orthog) {
    for (auto &Fc : orbitals) {
      if (Fi != nullptr && Fi->kappa() == Fc.kappa()) {
        Fc -= (*Fi * Fc) * *Fi;
      }
    }
  }

  return 0;
}

//==============================================================================
void ContinuumOrbitals::clear() { orbitals.clear(); }
