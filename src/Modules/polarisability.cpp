#include "polarisability.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/TDHF.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/StructureRad.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"
#include <iostream>
#include <string>
#include <vector>

namespace Module {
using namespace Polarisability;

void polarisability(const IO::InputBlock &input, const Wavefunction &wf) {

  std::cout << "\nDipole polarisability:\n";

  input.check(
      {{"rpa", "true/false; include RPA"},
       {"omega", "frequency (for single w)"},
       {"omega_max", "if >0, will run as function of w"},
       {"omega_steps", "number of w steps (from 0) if wmax>0"},
       {"StrucRadNorm", "[only for 'single' frequency]"},
       {"transition", "e.g., '6s,7s' or '6s+,7s+' - not fully checked"}});

  const auto rpaQ = input.get("rpa", true);
  const auto omega = input.get("omega", 0.0);

  const auto srQ = input.get("StrucRadNorm", false);

  // Generare E1 operator and TDHF object:
  const auto he1 = DiracOperator::E1(*wf.rgrid);
  auto dVE1 = ExternalField::TDHF(&he1, wf.getHF());

  // Solve TDHF for core, is doing RPA.
  // nb: even if not doing RPA, need TDHF object for tdhf method
  if (rpaQ) {
    dVE1.solve_core(omega);
  }

  // =================================================================
  // "Static" core + valence dipole polarisabilities:
  // (static if omega = 0), but single omega
  if (omega > 0.0) {
    std::cout << "At single frequency w=" << omega << "\n";
  } else {
    std::cout << "Static\n";
  }
  const auto ac_tdhf =
      alpha_core_tdhf(wf.core, he1, dVE1, omega, wf.getSigma());
  const auto ac_sos = alpha_core_sos(wf.core, wf.spectrum, he1, dVE1, omega);

  std::cout << "           TDHF        SOS\n";
  printf("Core: %9.3f  %9.3f\n", ac_tdhf, ac_sos);

  for (const auto &Fv : wf.valence) {
    const auto av_tdhf =
        alpha_valence_tdhf(Fv, Fv, he1, omega, dVE1, wf.getSigma());
    const auto av_sos = alpha_valence_sos(Fv, wf.spectrum, he1, dVE1, omega);
    printf("%4s: %9.3f  %9.3f", Fv.shortSymbol().c_str(), av_tdhf, av_sos);
    auto eps = std::abs((ac_tdhf + av_tdhf) / (ac_sos + av_sos) - 1.0);
    printf("  =  %9.3f  %9.3f   eps=%.0e\n", ac_tdhf + av_tdhf, ac_sos + av_sos,
           eps);
  }

  // =================================================================
  if (srQ) {
    std::cout << "\nIncluding Structure Radiation + Normalisation\n"
              << std::flush;
    for (const auto &Fv : wf.valence) {
      std::cout << Fv.symbol() << "\n";
      const auto delta_n_max_sum = 2;
      alpha_v_SRN(Fv, wf.spectrum, delta_n_max_sum, wf.basis,
                  wf.en_coreval_gap(), he1, dVE1, omega);
    }
  }

  // =================================================================
  // Transition polarisability, alpha and beta:
  const auto ab_vec = input.get<std::vector<std::string>>("transition", {});
  std::cout << ab_vec.size() << "\n";
  const auto ab_ok = (ab_vec.size() >= 2);
  const auto [na, ka] =
      ab_ok ? AtomData::parse_symbol(ab_vec[0]) : std::pair{0, 0};
  const auto [nb, kb] =
      ab_ok ? AtomData::parse_symbol(ab_vec[1]) : std::pair{0, 0};
  const auto Fa = wf.getState(na, ka);
  const auto Fb = wf.getState(nb, kb);
  if (Fa && Fb) {
    std::cout << "\nScalar dipole transition polarisability, alpha\n";

    const auto a_tdhf =
        alpha_valence_tdhf(*Fb, *Fa, he1, omega, dVE1, wf.getSigma());
    const auto a_sos =
        alpha_valence_sos(*Fb, *Fa, wf.spectrum, he1, dVE1, omega);

    std::cout << Fa->symbol() << " - " << Fb->symbol();
    printf(" : %9.3f  %9.3f\n", a_tdhf, a_sos);

    std::cout << "\nVector dipole transition polarisability, beta (sos)\n";
    const auto [b1, b2] = beta_sos(*Fb, *Fa, wf.spectrum, he1, dVE1, omega);
    std::cout << Fa->symbol() << " - " << Fb->symbol();
    printf(" : %9.3f + %9.3f = %9.3f\n", b1, b2, b1 + b2);
  }

  // =================================================================
  // Dynamic dipole polarisability, alpha:
  const auto omega_max = input.get("omega_max", 0.0);
  const auto omega_steps = std::abs(input.get("omega_steps", 30));
  const auto *const pFv = wf.valence.empty() ? nullptr : &wf.valence.front();
  if (omega_max > 0.0) {
    std::cout << "\nDynamic polarisability:\n";
    // std::cout << "ww    a(w)[TDHF] error  a(w)[SOS]  error\n";
    std::cout << "ww      eps    ac(w) tdhf         sos ";
    if (pFv)
      std::cout << "| ac+v(w) tdhf          sos";
    std::cout << "\n";

    const auto omega_list = qip::uniform_range(0.0, omega_max, omega_steps);
    for (const auto ww : omega_list) {
      if (dVE1.get_eps() > 1.0e-2) {
        // if tdhf didn't converge well last time, start from scratch
        // (otherwise, start from where we left off, since much faster)
        dVE1.clear();
      }
      if (rpaQ)
        dVE1.solve_core(ww, 30, false);
      const auto ac_w_tdhf = alpha_core_tdhf(wf.core, he1, dVE1, ww);
      const auto ac_w_sos = alpha_core_sos(wf.core, wf.spectrum, he1, dVE1, ww);
      printf("%6.4f  %.0e %11.4e %11.4e   ", ww, dVE1.get_eps(), ac_w_tdhf,
             ac_w_sos);
      if (pFv) {
        const auto av_w_tdhf =
            alpha_valence_tdhf(*pFv, *pFv, he1, ww, dVE1, wf.getSigma());
        const auto av_w_sos =
            alpha_valence_sos(*pFv, wf.spectrum, he1, dVE1, ww);

        printf(" %11.4e  %11.4e", ac_w_tdhf + av_w_tdhf, ac_w_sos + av_w_sos);
      }
      std::cout << "\n";
    }
  }
}

//******************************************************************************

namespace Polarisability {

//------------------------------------------------------------------------------
double alpha_core_tdhf(const std::vector<DiracSpinor> &core,
                       const DiracOperator::E1 &he1, ExternalField::TDHF &dVE1,
                       double omega,
                       const MBPT::CorrelationPotential *const Sigma) {
  auto alpha_core = 0.0;
  const auto f = (-1.0 / 3.0);
  for (const auto &Fb : core) {
    // this will include dV
    const auto Xb =
        dVE1.solve_dPsis(Fb, omega, ExternalField::dPsiType::X, Sigma);
    const auto Yb =
        dVE1.solve_dPsis(Fb, omega, ExternalField::dPsiType::Y, Sigma);
    for (const auto &Xbeta : Xb) {
      // no dV here (for closed-shell core)
      alpha_core += he1.reducedME(Xbeta, Fb);
      // alpha_core += he1.reducedME(Fb, Xbeta);
    }
    for (const auto &Ybeta : Yb) {
      alpha_core += he1.reducedME(Ybeta, Fb);
    }
  }
  return f * alpha_core;
}

//------------------------------------------------------------------------------
double alpha_valence_tdhf(const DiracSpinor &Fa, const DiracSpinor &Fb,
                          const DiracOperator::E1 &he1, double omega,
                          ExternalField::TDHF &dVE1,
                          const MBPT::CorrelationPotential *const Sigma) {
  double alpha_v = 0.0;
  // s-p only??
  const auto f = (-1.0 / 3.0) / (Fa.twoj() + 1);
  const auto Xa =
      dVE1.solve_dPsis(Fa, omega, ExternalField::dPsiType::X, Sigma);
  const auto Yb =
      dVE1.solve_dPsis(Fb, omega, ExternalField::dPsiType::Y, Sigma);
  for (const auto &Xalpha : Xa) {
    alpha_v += he1.reducedME(Xalpha, Fb) + dVE1.dV(Xalpha, Fb);
  }
  for (const auto &Ybeta : Yb) {
    alpha_v += he1.reducedME(Ybeta, Fa) + dVE1.dV(Ybeta, Fa);
  }
  return f * alpha_v;
}

//------------------------------------------------------------------------------
double alpha_core_sos(const std::vector<DiracSpinor> &core,
                      const std::vector<DiracSpinor> &basis,
                      const DiracOperator::E1 &he1, ExternalField::TDHF &dVE1,
                      double omega) {

  auto alpha_core = 0.0;
  const auto f = (-2.0 / 3.0);

  // core part: sum_{n,c} |<n|d|c>|^2, n excited states, c core states
  for (const auto &Fb : core) {
    for (const auto &Fn : basis) {
      // if (wf.isInCore(Fn.n, Fn.k))
      //   continue; // if core(HF) = core(basis), these cancel excactly
      if (he1.isZero(Fb.k, Fn.k))
        continue;
      const auto d1 = he1.reducedME(Fn, Fb);
      const auto d2 = d1 + dVE1.dV(Fn, Fb);
      const auto de = Fb.en() - Fn.en();
      alpha_core += std::abs(d1 * d2) * de / (de * de - omega * omega);
    }
  }

  return f * alpha_core;
}

//------------------------------------------------------------------------------
double alpha_valence_sos(const DiracSpinor &Fv,
                         const std::vector<DiracSpinor> &basis,
                         const DiracOperator::E1 &he1,
                         ExternalField::TDHF &dVE1, double omega) {

  auto alpha_v = 0.0;
  const auto f = (-2.0 / 3.0) / (Fv.twoj() + 1);

  for (const auto &Fn : basis) {
    // if (wf.isInCore(Fn.n, Fn.k))
    //   continue; // if core(HF) = core(basis), these cancel excactly
    if (he1.isZero(Fv.k, Fn.k))
      continue;
    const auto d2 = he1.reducedME(Fn, Fv) + dVE1.dV(Fn, Fv);
    // both have dV valence
    const auto de = Fv.en() - Fn.en();
    alpha_v += std::abs(d2 * d2) * de / (de * de - omega * omega);
  }

  return f * alpha_v;
}

//------------------------------------------------------------------------------
double alpha_v_SRN(const DiracSpinor &Fv,
                   const std::vector<DiracSpinor> &spectrum,
                   int delta_n_max_sum,
                   const std::vector<DiracSpinor> &hf_basis,
                   const double en_core, const DiracOperator::E1 &he1,
                   ExternalField::TDHF &dVE1, double omega) {

  // XXX Basis should be HF basis, NOT spectrum
  // "spectrum" may be valence or spectrum

  auto alpha_v0 = 0.0;
  auto alpha_v1 = 0.0;
  auto alpha_v2 = 0.0;
  const auto f = (-2.0 / 3.0) / (Fv.twoj() + 1);

  // std::unique_ptr<MBPT::StructureRad> sr(nullptr);
  // XXX nb: {3,30} no good for < Cs...
  auto sr = MBPT::StructureRad(hf_basis, en_core, {1, 99});

  std::cout << "       d0        dEn      |  a0(n)     da(HF)    da(RPA)  |  "
               "%(HF)     %(RPA)\n";
  for (const auto &Fn : spectrum) {
    if (Fn.en() < en_core)
      continue;
    // Only do for terms with small delta_n
    if (std::abs(Fn.n - Fv.n) > delta_n_max_sum)
      continue;
    if (he1.isZero(Fv.k, Fn.k))
      continue;
    const auto d0 = he1.reducedME(Fn, Fv) + dVE1.dV(Fn, Fv);

    const auto [tb, tbx] = sr.srTB(&he1, Fn, Fv, 0.0, &dVE1);
    const auto [c, cx] = sr.srC(&he1, Fn, Fv, &dVE1);
    const auto [n, nx] = sr.norm(&he1, Fn, Fv, &dVE1);
    const auto d1 = d0 + (tb + c + n);
    const auto d2 = d0 + (tbx + cx + nx);

    const auto de = Fv.en() - Fn.en();
    const auto da_v0 = f * std::abs(d0 * d0) * de / (de * de - omega * omega);
    const auto da_v1 = f * std::abs(d1 * d1) * de / (de * de - omega * omega);
    const auto da_v2 = f * std::abs(d2 * d2) * de / (de * de - omega * omega);
    alpha_v0 += da_v0;
    alpha_v1 += da_v1;
    alpha_v2 += da_v2;

    printf(" %4s %9.2e %9.2e | %9.2e %9.2e %9.2e | %8.1e%% %8.1e%%\n",
           Fn.shortSymbol().c_str(), d0, de, da_v0, da_v1 - da_v0,
           da_v2 - da_v0, (da_v1 - da_v0) / da_v0 * 100,
           (da_v2 - da_v0) / da_v0 * 100);
    std::cout << std::flush;
  }

  const auto srn = (alpha_v1 - alpha_v0);
  const auto srn_x = (alpha_v2 - alpha_v0);

  std::cout << "a(main): " << alpha_v0 << "\n";
  std::cout << "StrucRad+Norm:\n";
  std::cout << "No RPA:  " << srn << " " << srn / alpha_v0 << "\n";
  std::cout << "w/ RPA:  " << srn_x << " " << srn_x / alpha_v0 << "\n";
  std::cout << std::flush;

  return srn;
}

//------------------------------------------------------------------------------
double alpha_valence_sos(const DiracSpinor &Fv, const DiracSpinor &Fw,
                         const std::vector<DiracSpinor> &basis,
                         const DiracOperator::E1 &he1,
                         ExternalField::TDHF &dVE1, double omega) {

  (void)omega; // XXX Where omega go??

  std::cout << "\nWARNING 312: Missing omega in formula??\n";

  auto alpha_v = 0.0;
  const auto f = (-1.0 / 3.0) / (Fv.twoj() + 1);

  // core part: sum_{n,c} |<n|d|c>|^2, n excited states, c core states
  for (const auto &Fn : basis) {
    // if (wf.isInCore(Fn.n, Fn.k))
    //   continue; // if core(HF) = core(basis), these cancel excactly
    if (he1.isZero(Fv.k, Fn.k))
      continue;
    const auto d1 = he1.reducedME(Fn, Fv) + dVE1.dV(Fn, Fv);
    const auto d2 = Fv == Fn ? d1 : he1.reducedME(Fn, Fw) + dVE1.dV(Fn, Fw);
    // both have dV valence
    const auto dev = Fv.en() - Fn.en();
    const auto dew = Fw.en() - Fn.en();
    const auto da = f * d1 * d2 * (1.0 / dew + 1.0 / dev);
    alpha_v += da;
    // std::cout << Fn.symbol() << d1 << " " << d2 << " " << da << "\n";
  }

  return alpha_v;
}

//------------------------------------------------------------------------------
std::pair<double, double> beta_sos(const DiracSpinor &Fv, const DiracSpinor &Fw,
                                   const std::vector<DiracSpinor> &basis,
                                   const DiracOperator::E1 &he1,
                                   ExternalField::TDHF &dVE1, double omega) {

  (void)omega; // XXX Where omega go??

  // XXX ONly for s-s beta!

  auto beta_1 = 0.0;
  auto beta_2 = 0.0;
  const auto f = (-1.0 / 3.0) / (Fv.twoj() + 1);

  // core part: sum_{n,c} |<n|d|c>|^2, n excited states, c core states
  for (const auto &Fn : basis) {
    // if (wf.isInCore(Fn.n, Fn.k))
    //   continue; // if core(HF) = core(basis), these cancel excactly
    if (he1.isZero(Fv.k, Fn.k))
      continue;

    const auto d1 = he1.reducedME(Fn, Fv) + dVE1.dV(Fn, Fv);
    const auto d2 = Fv == Fn ? d1 : he1.reducedME(Fn, Fw) + dVE1.dV(Fn, Fw);
    // both have dV valence
    const auto dev = Fv.en() - Fn.en();
    const auto dew = Fw.en() - Fn.en();
    const auto dB = f * d1 * d2 * (1.0 / dew + 1.0 / dev);
    if (Fn.twoj() == Fv.twoj())
      beta_1 += dB;
    else
      beta_2 += -0.5 * dB;
  }

  return {beta_1, beta_2};
}

} // namespace Polarisability
} // namespace Module
