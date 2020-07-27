#include "Modules/Module_HFAnomaly.hpp"
#include "DiracOperator/Operators.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/UserInput.hpp"
#include "MBPT/DiagramRPA.hpp"
#include "Modules/Module_matrixElements.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace Module {

void HFAnomaly(const IO::UserInputBlock &input, const Wavefunction &wf) {

  const auto rpa = input.get("rpa", true);

  const auto sub_input = input.subBlock("MatrixElements::hfs",
                                        {"mu", "I", "rrms", "parity", "l", "gl",
                                         "mu1", "gl1", "l1", "l2", "I1", "I2"});
  const auto point_in = sub_input.copy_with("F(r)=pointlike");

  const auto ball_in = sub_input.copy_with("F(r)=ball");
  const auto doubly_odd = (wf.Anuc() % 2 == 0);
  const auto sp_in = doubly_odd ? sub_input.copy_with("F(r)=doublyOddBW")
                                : sub_input.copy_with("F(r)=VolotkaBW");

  auto hpt = generateOperator(point_in, wf, false);
  auto hbl = generateOperator(ball_in, wf, false);
  auto hsp = generateOperator(sp_in, wf, true);

  // calc mu, I, gI
  struct ThreeAs {
    ThreeAs(double a, double b, double c) : pt(a), bl(b), sp(c){};
    double pt, bl, sp;
  };

  auto isotope0 = Nuclear::findIsotopeData(wf.Znuc(), wf.Anuc());
  auto mu0 = input.get("mu", isotope0.mu);
  auto I0 = input.get("I", isotope0.I_N);
  auto r_rms0 = input.get("rrms", isotope0.r_rms);
  auto gI0 = mu0 / I0;

  // input.get
  std::unique_ptr<MBPT::DiagramRPA> rpap{nullptr}, rpab{nullptr}, rpas{nullptr};
  if (rpa) {
    std::cout << "\nIncluding RPA (diagram method) - must have basis\n";
    rpap = std::make_unique<MBPT::DiagramRPA>(hpt.get(), wf.basis, wf.core,
                                              wf.identity());
    rpab = std::make_unique<MBPT::DiagramRPA>(hbl.get(), rpap.get());
    rpas = std::make_unique<MBPT::DiagramRPA>(hsp.get(), rpap.get());
    std::cout << "Solving RPA core for point, ball, SP:\n";
    rpap->rpa_core(0.0, true);
    rpab->rpa_core(0.0, true);
    rpas->rpa_core(0.0, true);
  }

  std::vector<ThreeAs> As;
  std::cout << "\nReference: A = " << wf.Anuc() << ", gI = " << gI0
            << ", r_rms = " << r_rms0 << "\n\n";
  std::cout << "      A(point)     A(ball)      A(SP)       | e(ball)  e(SP)\n";
  for (const auto &Fv : wf.valence) {
    auto point = DiracOperator::Hyperfine::hfsA(hpt.get(), Fv);
    auto ball = DiracOperator::Hyperfine::hfsA(hbl.get(), Fv);
    auto sp = DiracOperator::Hyperfine::hfsA(hsp.get(), Fv);
    if (rpa) {
      auto a = DiracOperator::Hyperfine::convertRMEtoA(Fv, Fv);
      point += a * rpap->dV(Fv, Fv);
      ball += a * rpab->dV(Fv, Fv);
      sp += a * rpas->dV(Fv, Fv);
    }
    As.emplace_back(point, ball, sp);
    const auto eps_ball = 100.0 * (ball - point) / point;
    const auto eps_sp = 100.0 * (sp - point) / point;
    printf("%-4s %12.5e %12.5e %12.5e | %7.4f  %7.4f\n",
           Fv.shortSymbol().c_str(), point, ball, sp, eps_ball, eps_sp);
  }

  const auto isotopes = Nuclear::findIsotopeList(wf.Znuc());
  for (const auto &nuc : isotopes) {
    // only do odd isotopes (even must be reference!)
    if (nuc.A % 2 == 0)
      continue;
    auto wfA = Wavefunction(wf.rgrid->params(), {nuc.Z, nuc.A},
                            wf.alpha / PhysConst::alpha);
    const auto gI = nuc.mu / nuc.I_N;
    std::cout << "\nA = " << nuc.A << ", gI = " << gI
              << ", r_rms = " << nuc.r_rms << "\n";
    wfA.hartreeFockCore("HartreeFock", 0.0, wf.coreConfiguration_nice());
    wfA.hartreeFockValence(DiracSpinor::state_config(wf.valence));
    wfA.basis = wf.basis; // OK??

    auto in1 = IO::UserInputBlock("MatrixElements::hfs", {"F(r)=pointlike"});
    auto in2 = IO::UserInputBlock("MatrixElements::hfs", {"F(r)=ball"});
    auto in3 = IO::UserInputBlock("MatrixElements::hfs", {"F(r)=VolotkaBW"});
    auto hpt2 = generateOperator(in1, wfA, false);
    auto hbl2 = generateOperator(in2, wfA, false);
    auto hsp2 = generateOperator(in3, wfA, false);

    std::unique_ptr<MBPT::DiagramRPA> rpap2{nullptr}, rpab2{nullptr},
        rpas2{nullptr};
    if (rpa) {
      std::cout << "\nIncluding RPA (diagram method) - must have basis\n";
      // OK ? or need run with new basis!?
      rpap2 = std::make_unique<MBPT::DiagramRPA>(hpt2.get(), wfA.basis,
                                                 wfA.core, wfA.identity());
      rpab2 = std::make_unique<MBPT::DiagramRPA>(hbl2.get(), rpap.get());
      rpas2 = std::make_unique<MBPT::DiagramRPA>(hsp2.get(), rpap.get());
      std::cout << "Solving RPA core for point, ball, SP:\n";
      rpap2->rpa_core(0.0, true);
      rpab2->rpa_core(0.0, true);
      rpas2->rpa_core(0.0, true);
    }

    // for (const auto &Fv : wf.valence)
    std::cout << "\n      A(point)      e(ball)   e(SP)   | 1D2(ball) "
                 "1D2(SP)  [%]\n";
    for (std::size_t i = 0; i < wf.valence.size(); ++i) {
      const auto &Fv = wf.valence[i];
      const auto [pt0, bl0, sp0] = As[i];
      auto point = DiracOperator::Hyperfine::hfsA(hpt2.get(), Fv);
      auto ball = DiracOperator::Hyperfine::hfsA(hbl2.get(), Fv);
      auto sp = DiracOperator::Hyperfine::hfsA(hsp2.get(), Fv);
      if (rpa) {
        auto a = DiracOperator::Hyperfine::convertRMEtoA(Fv, Fv);
        point += a * rpap2->dV(Fv, Fv);
        ball += a * rpab2->dV(Fv, Fv);
        sp += a * rpas2->dV(Fv, Fv);
      }
      // As.emplace_back(point, ball, sp);
      const auto eps_ball = 100.0 * (ball - point) / point;
      const auto eps_sp = 100.0 * (sp - point) / point;
      auto D12_bl = 100.0 * ((bl0 / ball) * (gI / gI0) - 1.0);
      auto D12_sp = 100.0 * ((sp0 / sp) * (gI / gI0) - 1.0);
      if (std::abs(eps_ball) > 0.001)
        printf("%-4s %12.5e | %8.5f %8.5f | %8.5f  %8.5f\n",
               Fv.shortSymbol().c_str(), point, eps_ball, eps_sp, D12_bl,
               D12_sp);
      else
        printf("%-4s %12.5e | %8.1e %8.1e | %8.1e  %8.1e\n",
               Fv.shortSymbol().c_str(), point, eps_ball, eps_sp, D12_bl,
               D12_sp);
    }
  }
}

} // namespace Module
