#include "Amplitudes/MatrixElements.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/TDHF.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/StructureRad.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <map>
#include <sstream>
#include <string>
#include <vector>

//==============================================================================
TEST_CASE("Amplitudes: matrix elements", "[Amplitudes][unit]") {

  Wavefunction wf({500, 1.0e-4, 80.0, 20.0, "loglinear", -1.0},
                  {"Li", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[He]");
  wf.solve_valence("2sp");

  auto dE1 = DiracOperator::E1(wf.grid());

  // Single-pair function: t0 is the reduced ME; factor is matel_factor
  for (const auto &a : wf.valence()) {
    for (const auto &b : wf.valence()) {
      if (dE1.isZero(a, b))
        continue;
      const auto me = Amplitudes::matrix_element(a, b, &dE1);
      REQUIRE(me.t0 == Approx(dE1.reducedME(a, b)));
      REQUIRE(me.factor == 1.0);
      REQUIRE(me.dv == 0.0);
      REQUIRE(!me.has_rpa);

      using DiracOperator::MatrixElementType;
      for (const auto type :
           {MatrixElementType::Stretched, MatrixElementType::HFConstant}) {
        const auto me2 =
          Amplitudes::matrix_element(a, b, &dE1, nullptr, nullptr, type);
        REQUIRE(me2.factor == Approx(dE1.matel_factor(type, a, b)));
        REQUIRE(me2.value0() == Approx(me2.factor * me.t0));
      }
    }
  }

  // List driver (no RPA): each element is the reduced ME of its pair, and
  // every allowed pair appears exactly once
  {
    using Amplitudes::Frequency;
    Amplitudes::MEoptions options{Frequency::transition, Frequency::fixed};
    options.print = false;
    const auto mes = Amplitudes::matrix_elements(wf.valence(), &dE1, nullptr,
                                                 nullptr, options);

    REQUIRE(!mes.empty());
    std::map<std::pair<std::string, std::string>, int> count;
    for (const auto &m : mes) {
      const auto *Fa = wf.getState(m.a);
      const auto *Fb = wf.getState(m.b);
      REQUIRE(Fa != nullptr);
      REQUIRE(Fb != nullptr);
      REQUIRE(m.t0 == Approx(dE1.reducedME(*Fa, *Fb)));
      REQUIRE(m.omega == Approx(Fa->en() - Fb->en()));
      REQUIRE(m.dv == 0.0);
      ++count[{m.a, m.b}];
    }
    for (const auto &[ab, n] : count) {
      REQUIRE(n == 1);
    }

    // Two-list overload: every allowed pair, each calculated once
    const auto mes2 = Amplitudes::matrix_elements(
      wf.valence(), wf.valence(), &dE1, nullptr, nullptr, options);
    REQUIRE(mes2.size() >= mes.size());
    for (const auto &m : mes2) {
      const auto *Fa = wf.getState(m.a);
      const auto *Fb = wf.getState(m.b);
      REQUIRE(m.t0 == Approx(dE1.reducedME(*Fa, *Fb)));
    }
  }

  // Frequency-dependent operator (E1v): driver evaluates each pair at its
  // own (signed) transition frequency, using h_minus for negative omega
  {
    auto hv = DiracOperator::E1v(wf.alpha(), 0.0);
    auto hv_minus = DiracOperator::E1v(wf.alpha(), 0.0);

    using Amplitudes::Frequency;
    Amplitudes::MEoptions options{Frequency::transition, Frequency::fixed};
    options.print = false;
    options.calculate_both = true;
    const auto mes = Amplitudes::matrix_elements(wf.valence(), &hv, &hv_minus,
                                                 nullptr, options);

    REQUIRE(!mes.empty());
    for (const auto &m : mes) {
      const auto *Fa = wf.getState(m.a);
      const auto *Fb = wf.getState(m.b);
      REQUIRE(Fa != nullptr);
      REQUIRE(Fb != nullptr);
      auto hv_w = DiracOperator::E1v(wf.alpha(), m.omega);
      REQUIRE(m.t0 == Approx(hv_w.reducedME(*Fa, *Fb)));
    }
  }

  // With RPA: dv is dV, solved by the driver
  {
    auto rpa = ExternalField::TDHF(&dE1, wf.vHF());
    using Amplitudes::Frequency;
    Amplitudes::MEoptions options{Frequency::transition, Frequency::fixed};
    options.print = false;
    rpa.solve_core(0.0, options.rpa_iterations, false);
    const auto mes =
      Amplitudes::matrix_elements(wf.valence(), &dE1, nullptr, &rpa, options);

    auto rpa0 = ExternalField::TDHF(&dE1, wf.vHF());
    rpa0.solve_core(0.0, options.rpa_iterations, false);

    for (const auto &m : mes) {
      const auto *Fa = wf.getState(m.a);
      const auto *Fb = wf.getState(m.b);
      REQUIRE(m.has_rpa);
      REQUIRE(m.dv == Approx(rpa0.dV(*Fa, *Fb)));
    }
  }
}

//==============================================================================
TEST_CASE("Amplitudes: matrix elements with RPA",
          "[Amplitudes][ExternalField][integration]") {

  Wavefunction wf({1000, 1.0e-5, 100.0, 20.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Xe]");
  wf.solve_valence("6sp");

  fmt::print("\nAmplitudes::matrix_elements: E1 with RPA (Cs)\n");
  fmt::print("{:5s} {:5s} {:>13s} {:>13s} {:>9s}\n", "a", "b", "expected",
             "found", "eps");

  // E1, TDHF RPA at fixed w=0: against a separately-solved RPA
  {
    auto dE1 = DiracOperator::E1(wf.grid());

    auto rpa = ExternalField::TDHF(&dE1, wf.vHF());
    using Amplitudes::Frequency;
    Amplitudes::MEoptions options{Frequency::transition, Frequency::fixed};
    options.print = false;
    rpa.solve_core(0.0, options.rpa_iterations, false);
    const auto mes =
      Amplitudes::matrix_elements(wf.valence(), &dE1, nullptr, &rpa, options);

    auto rpa0 = ExternalField::TDHF(&dE1, wf.vHF());
    rpa0.solve_core(0.0, options.rpa_iterations, false);

    REQUIRE(!mes.empty());
    for (const auto &m : mes) {
      const auto *Fa = wf.getState(m.a);
      const auto *Fb = wf.getState(m.b);
      REQUIRE(Fa != nullptr);
      REQUIRE(Fb != nullptr);
      const auto expected = dE1.reducedME(*Fa, *Fb) + rpa0.dV(*Fa, *Fb);
      const auto found = m.value();
      const auto eps = std::abs((found - expected) / expected);
      fmt::print("{:5s} {:5s} {:13.6e} {:13.6e} {:9.1e}\n", m.a, m.b, expected,
                 found, eps);
      REQUIRE(found == Approx(expected).epsilon(1.0e-10));
    }
  }

  // hfs (HFConstant): factor * (t0 + dv) against direct evaluation
  {
    const auto h = DiracOperator::generate("hfs", {"hfs", "print=false;"}, wf);

    auto rpa = ExternalField::TDHF(h.get(), wf.vHF());
    using Amplitudes::Frequency;
    Amplitudes::MEoptions options{Frequency::transition, Frequency::fixed};
    options.print = false;
    options.rpa_iterations = 1; // first-order RPA: quick
    options.type = DiracOperator::MatrixElementType::HFConstant;
    rpa.solve_core(0.0, options.rpa_iterations, false);
    const auto mes = Amplitudes::matrix_elements(wf.valence(), h.get(), nullptr,
                                                 &rpa, options);

    REQUIRE(!mes.empty());
    for (const auto &m : mes) {
      const auto *Fa = wf.getState(m.a);
      const auto *Fb = wf.getState(m.b);
      const auto factor =
        h->matel_factor(DiracOperator::MatrixElementType::HFConstant, *Fa, *Fb);
      const auto expected =
        factor * (h->reducedME(*Fa, *Fb) + rpa.dV(*Fa, *Fb));
      const auto found = m.value();
      const auto eps = std::abs((found - expected) / expected);
      fmt::print("{:5s} {:5s} {:13.6e} {:13.6e} {:9.1e}\n", m.a, m.b, expected,
                 found, eps);
      REQUIRE(found == Approx(expected).epsilon(1.0e-10));
    }
  }

  // Frequency::transition: RPA re-solved at each transition frequency
  {
    auto dE1 = DiracOperator::E1(wf.grid());

    auto rpa = ExternalField::TDHF(&dE1, wf.vHF());
    using Amplitudes::Frequency;
    Amplitudes::MEoptions options{Frequency::transition, Frequency::transition};
    options.print = false;
    const auto mes =
      Amplitudes::matrix_elements(wf.valence(), &dE1, nullptr, &rpa, options);

    auto rpa0 = ExternalField::TDHF(&dE1, wf.vHF());
    for (const auto &m : mes) {
      const auto *Fa = wf.getState(m.a);
      const auto *Fb = wf.getState(m.b);
      rpa0.clear();
      rpa0.solve_core(m.omega, options.rpa_iterations, false);
      const auto expected = rpa0.dV(*Fa, *Fb);
      const auto found = m.dv;
      const auto eps = std::abs((found - expected) / expected);
      fmt::print("{:5s} {:5s} {:13.6e} {:13.6e} {:9.1e}\n", m.a, m.b, expected,
                 found, eps);
      // driver may iterate from previous solution: small differences allowed
      REQUIRE(found == Approx(expected).epsilon(1.0e-5));
    }
  }
}

//==============================================================================
//! sr_matrix_elements against direct evaluation with MBPT::StructureRad
TEST_CASE("Amplitudes: SR matrix elements",
          "[Amplitudes][StrucRad][MBPT][unit]") {

  Wavefunction wf({400, 1.0e-4, 50.0, 33.0, "loglinear"}, {"Na", -1, "Fermi"},
                  1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-5);
  wf.solve_valence("3sp");
  // wf.formBasis({"20spd", 20, 6, {3.4e-4, 1.7e-2}, 30.0});
  wf.formBasis({"20spd", 20, 6, {1.0e-4}, 30.0});

  auto hE1 = DiracOperator::E1(wf.grid());

  MBPT::StructureRad sr(wf.basis(), wf.FermiLevel(), {2, 10}, "", 99, {}, {},
                        false);
  using Amplitudes::Frequency;
  Amplitudes::SRNoptions options{Frequency::transition, Frequency::fixed};
  options.print = false;
  const auto mes =
    Amplitudes::sr_matrix_elements(wf.valence(), &hE1, &sr, nullptr, options);

  // Direct evaluation, separate StructureRad object:
  MBPT::StructureRad sr0(wf.basis(), wf.FermiLevel(), {2, 10}, "", 99, {}, {},
                         false);
  sr0.solve_core(&hE1, nullptr);

  REQUIRE(!mes.empty());
  for (const auto &m : mes) {
    const auto *Fa = wf.getState(m.a);
    const auto *Fb = wf.getState(m.b);
    REQUIRE(Fa != nullptr);
    REQUIRE(Fb != nullptr);
    REQUIRE(m.t0 == Approx(hE1.reducedME(*Fa, *Fb)));
    // Frequency::fixed: SR denominators at the RPA frequency (0), not w_ab
    REQUIRE(m.sr == Approx(sr0.SR(*Fa, *Fb, 0.0)));
    REQUIRE(m.norm == Approx((sr0.f_norm(*Fa) + sr0.f_norm(*Fb)) * m.t0));
    REQUIRE(m.bo == Approx(sr0.BO(*Fa, *Fb)));
    REQUIRE(m.total() ==
            Approx(m.factor * (m.t0 + m.dv + m.sr + m.norm + m.bo)));
  }
}
