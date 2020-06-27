#pragma once
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/Vector.hpp"
#include <string>

namespace UnitTest {
bool HartreeFock(std::ostream &obuff) {
  bool pass = true;

  // Create wavefunction object
  Wavefunction wf({2000, 1.0e-6, 120.0, 0.33 * 120.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.hartreeFockCore("HartreeFock", 0.0, "[Xe]");
  wf.hartreeFockValence("6sp5d");

  Wavefunction wf2({10000, 1.0e-7, 200.0, 0.33 * 150.0, "loglinear", -1.0},
                   {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf2.hartreeFockCore("HartreeFock", 0.0, "[Xe]");
  wf2.hartreeFockValence("6sp5d");

  // test core energies
  const std::vector<double> expected_VD = {
      -1330.118692, -212.5644584, -199.4294701, -186.4365217, -45.96974213,
      -40.44831295, -37.89428811, -28.30949791, -27.77513323, -9.512822,
      -7.44628727,  -6.92099641,  -3.48561693,  -3.39689621,  -1.48980554,
      -0.90789792,  -0.84033918};

  const std::vector<double> expected_VD_v = {
      -0.12736810, -0.08561590, -0.08378548, -0.06441963, -0.06452976};

  auto cmpr = [](const auto &ds, const auto &v) { return (ds.en - v) / v; };
  auto cmpr2 = [](const auto &ds, const auto &ds2) {
    return (ds.en - ds2.en) / ds2.en;
  };

  // Check HF core:
  {
    const auto [eps, at] = qip::compare(wf.core, expected_VD, cmpr);
    const std::string worst = at == wf.core.end() ? "" : at->symbol();
    pass &= qip::check_value(&obuff, "HF core Cs vs. VD " + worst, eps, 0.0,
                             5.0e-6);
  }
  {
    const auto [eps, at] = qip::compare(wf.core, wf2.core, cmpr2);
    const std::string worst = at == wf.core.end() ? "" : at->symbol();
    pass &=
        qip::check_value(&obuff, "HF core Cs grid " + worst, eps, 0.0, 5.0e-6);
  }

  // Check Valence:
  {
    const auto [eps, at] = qip::compare(wf.valence, expected_VD_v, cmpr);
    const std::string worst = at == wf.valence.end() ? "" : at->symbol();
    pass &=
        qip::check_value(&obuff, "HF val Cs vs. VD " + worst, eps, 0.0, 3.0e-6);
  }
  {
    const auto [eps, at] = qip::compare(wf.valence, wf2.valence, cmpr2);
    const std::string worst = at == wf.valence.end() ? "" : at->symbol();
    pass &=
        qip::check_value(&obuff, "HF val Cs grid " + worst, eps, 0.0, 3.0e-6);
  }

  // Check Breit (valence):
  {
    Wavefunction wf_Br({2000, 1.0e-6, 120.0, 0.33 * 120.0, "loglinear", -1.0},
                       {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf_Br.hartreeFockCore("HartreeFock", 1.0, "[Xe]"); // incl Breit
    wf_Br.hartreeFockValence("6sp5d");

    const std::vector<double> VD_vBr = {-0.12735353, -0.08558175, -0.0837724,
                                        -0.06446589, -0.06458303};
    const auto br_vlad = qip::compose([](auto a, auto b) { return a - b; },
                                      VD_vBr, expected_VD_v);

    std::vector<double> me_br;
    for (auto i = 0ul; i < wf.valence.size(); ++i) {
      me_br.push_back(wf_Br.valence[i].en - wf.valence[i].en);
    }

    const auto [eps, at] = qip::compare_eps(me_br, br_vlad);
    pass &= qip::check_value(&obuff, "HF Breit Cs val", eps, 0.0, 5.0e-4);
  }

  //***************************************
  // Yb+
  {
    Wavefunction wf3({2000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear", -1.0},
                     {"Yb", 173, "Fermi", -1.0, -1.0}, 1.0);
    wf3.hartreeFockCore("HartreeFock", 0.0, "[Xe],4f14");
    wf3.hartreeFockValence("6sp5d5f");

    const std::vector<double> dzuba = {-0.41366446, -0.30111301, -0.28830673,
                                       -0.30307172, -0.30088656, -0.12509665,
                                       -0.12507894};
    const auto [eps, at] = qip::compare(wf3.valence, dzuba, cmpr);
    pass &= qip::check_value(&obuff, "HF Yb+ val", eps, 0.0, 7.0e-6);
  }

  return pass;
}

} // namespace UnitTest
