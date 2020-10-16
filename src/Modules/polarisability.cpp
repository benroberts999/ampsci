#include "polarisability.hpp"
#include "DiracOperator/Operators.hpp"
#include "HF/ExternalField.hpp"
#include "IO/UserInput.hpp"
// #include "MBPT/CorrelationPotential.hpp"
#include "Wavefunction/Wavefunction.hpp"
// #include <algorithm>
// #include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace Module {
using namespace Polarisability;

void polarisability(const IO::UserInputBlock &input, const Wavefunction &wf) {

  std::cout << "\nDipole polarisability:\n";

  input.checkBlock({"rpa", "omega", "transition", "omega_max", "omega_steps"});

  const auto rpaQ = input.get("rpa", true);
  const auto omega = input.get("omega", 0.0);

  // Generare E1 operator and TDHF object:
  const auto he1 = DiracOperator::E1(*wf.rgrid);
  auto dVE1 = HF::ExternalField(&he1, wf.getHF());

  // Solve TDHF for core, is doing RPA.
  // nb: even if not doing RPA, need TDHF object for tdhf method
  if (rpaQ)
    dVE1.solve_TDHFcore(omega);

  // **************
  // "Static" core + valence dipole polarisabilities:
  // (static if omega = 0), but single omega
  const auto ac_tdhf = alpha_core_tdhf(wf.core, he1, dVE1, omega);
  const auto ac_sos = alpha_core_sos(wf.core, wf.spectrum, he1, dVE1, omega);

  std::cout << "           TDHF        SOS\n";
  printf("Core: %9.3f  %9.3f\n", ac_tdhf, ac_sos);

  for (const auto Fv : wf.valence) {
    const auto av_tdhf =
        alpha_valence_tdhf(Fv, Fv, he1, omega, dVE1, wf.getSigma());
    const auto av_sos = alpha_valence_sos(Fv, wf.spectrum, he1, dVE1, omega);
    printf("%4s: %9.3f  %9.3f", Fv.shortSymbol().c_str(), av_tdhf, av_sos);
    printf("  =  %9.3f  %9.3f\n", ac_tdhf + av_tdhf, ac_sos + av_sos);
  }

  // **************
  // Transition polarisability, alpha and beta:
  const auto ab_vec = input.get_list<std::string>("transition", {});
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

  // **************
  // Dynamic dipole polarisability, alpha and beta:
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
    const auto dw = omega_max / omega_steps;
    for (auto ww = 0.0; ww < omega_max; ww += dw) {
      if (dVE1.get_eps() > 1.0e-2) {
        // if tdhf didn't converge well last time, start from scratch
        // (otherwise, start from where we left off, since much faster)
        dVE1.clear_dPsi();
      }
      if (rpaQ)
        dVE1.solve_TDHFcore(ww, 30, false);
      const auto ac_w_tdhf = alpha_core_tdhf(wf.core, he1, dVE1, ww);
      const auto ac_w_sos = alpha_core_sos(wf.core, wf.spectrum, he1, dVE1, ww);
      printf("%6.4f  %.0e %11.4e %11.4e   ", ww, dVE1.get_eps(), ac_w_tdhf,
             ac_w_sos);
      if (pFv) {
        const auto av_w_tdhf =
            alpha_valence_tdhf(*pFv, *pFv, he1, ww, dVE1, wf.getSigma());
        const auto av_w_sos =
            alpha_valence_sos(*pFv, wf.spectrum, he1, dVE1, ww);

        printf(" %11.4e  %11.4e", ac_w_tdhf + av_w_tdhf, av_w_sos);
      }
      std::cout << "\n";
    }
  }
}

//******************************************************************************

namespace Polarisability {

//------------------------------------------------------------------------------
double alpha_core_tdhf(const std::vector<DiracSpinor> &core,
                       const DiracOperator::E1 &he1, HF::ExternalField &dVE1,
                       double omega) {
  auto alpha_core = 0.0;
  const auto f = (-1.0 / 3.0);
  for (const auto &Fb : core) {
    // this will include dV
    const auto Xb = dVE1.solve_dPsis(Fb, omega, HF::dPsiType::X);
    const auto Yb = dVE1.solve_dPsis(Fb, omega, HF::dPsiType::Y);
    for (const auto &Xbeta : Xb) {
      // no dV here (for closed-shell core)
      alpha_core += he1.reducedME(Xbeta, Fb);
    }
    for (const auto &Ybeta : Yb) {
      alpha_core += he1.reducedME(Fb, Ybeta);
    }
  }
  return f * alpha_core;
}

//------------------------------------------------------------------------------
double alpha_valence_tdhf(const DiracSpinor &Fa, const DiracSpinor &Fb,
                          const DiracOperator::E1 &he1, double omega,
                          HF::ExternalField &dVE1,
                          const MBPT::CorrelationPotential *const Sigma) {
  double alpha_v = 0.0;
  // s-p only??
  const auto f = (-1.0 / 3.0) / (Fa.twoj() + 1);
  const auto Xa = dVE1.solve_dPsis(Fa, omega, HF::dPsiType::X, Sigma);
  const auto Yb = dVE1.solve_dPsis(Fb, omega, HF::dPsiType::Y, Sigma);
  for (const auto &Xalpha : Xa) {
    alpha_v += he1.reducedME(Xalpha, Fb) + dVE1.dV(Xalpha, Fb);
  }
  for (const auto &Ybeta : Yb) {
    alpha_v += he1.reducedME(Fa, Ybeta) + dVE1.dV(Fa, Ybeta);
  }
  return f * alpha_v;
}

//------------------------------------------------------------------------------
double alpha_core_sos(const std::vector<DiracSpinor> &core,
                      const std::vector<DiracSpinor> &basis,
                      const DiracOperator::E1 &he1, HF::ExternalField &dVE1,
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
      const auto de = Fb.en - Fn.en;
      alpha_core += std::abs(d1 * d2) * de / (de * de - omega * omega);
    }
  }

  return f * alpha_core;
}

//------------------------------------------------------------------------------
double alpha_valence_sos(const DiracSpinor &Fv,
                         const std::vector<DiracSpinor> &basis,
                         const DiracOperator::E1 &he1, HF::ExternalField &dVE1,
                         double omega) {

  auto alpha_v = 0.0;
  const auto f = (-2.0 / 3.0) / (Fv.twoj() + 1);

  // core part: sum_{n,c} |<n|d|c>|^2, n excited states, c core states

  for (const auto &Fn : basis) {
    // if (wf.isInCore(Fn.n, Fn.k))
    //   continue; // if core(HF) = core(basis), these cancel excactly
    if (he1.isZero(Fv.k, Fn.k))
      continue;
    const auto d2 = he1.reducedME(Fn, Fv) + dVE1.dV(Fn, Fv);
    // both have dV valence
    const auto de = Fv.en - Fn.en;
    alpha_v += std::abs(d2 * d2) * de / (de * de - omega * omega);
  }

  return f * alpha_v;
}

//------------------------------------------------------------------------------
double alpha_valence_sos(const DiracSpinor &Fv, const DiracSpinor &Fw,
                         const std::vector<DiracSpinor> &basis,
                         const DiracOperator::E1 &he1, HF::ExternalField &dVE1,
                         double omega) {

  (void)omega; // XXX Where omega go??

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
    const auto dev = Fv.en - Fn.en;
    const auto dew = Fw.en - Fn.en;
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
                                   HF::ExternalField &dVE1, double omega) {

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
    const auto dev = Fv.en - Fn.en;
    const auto dew = Fw.en - Fn.en;
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
