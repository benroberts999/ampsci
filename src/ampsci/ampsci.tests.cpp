#include "ampsci.hpp"
#include "DiracOperator/GenerateOperator.hpp"
#include "ExternalField/include.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "catch2/catch.hpp"
#include "qip/Random.hpp"
#include "qip/Vector.hpp"

//============================================================================
TEST_CASE("ampsci - basic unit test",
          "[ampsci][HartreeFock][Correlations][json][unit]") {

  // This is just a test case that failed to converge before an update
  std::cout << " Simple 'unit' test with small parameters only. Not meant to "
               "be accurate\n";

  std::string run_label = "deleteme_" + qip::random_string(5);

  const auto input_string = std::string{R"(
    Atom {
      Z = Cs;
      run_label = )" + run_label + R"(;
      json_out = true;
    }
    HartreeFock {
      core = [Xe];
      valence = 6sp;
    }
    Grid {
      r0 = 1.0e-5;
      rmax = 75.0;
      num_points = 1000;
    }
    Basis {
      number = 30;
      order = 7;
      rmax = 35.0;
      states = 25spdfg;
    }
    Correlations {
      n_min_core = 3;
      each_valence = false;
    }
    Module::MatrixElements {
      operator = hfs;
      rpa = diagram;
      off-diagonal = false;
    }
    Module::MatrixElements {
      operator = E1;
      rpa = TDHF;
      omega = 0.0;
    }
  )"};

  const auto wf = ampsci(IO::InputBlock{"", input_string});

  REQUIRE(wf.Znuc() == 55);
  REQUIRE(wf.identity() == "Cs1_" + run_label);

  // require all valence states converge
  for (const auto &v : wf.valence()) {
    REQUIRE(v.eps() < 1.0e-9);
  }

  // Open and read the JSON output file
  std::ifstream file("Cs1_" + run_label + ".json");
  nlohmann::json json_file;
  file >> json_file;

  // Check what we read in matched: energies etc.
  for (const auto &c : wf.core()) {
    REQUIRE(json_file["wavefunctions"]["core"][c.shortSymbol()]["en"] ==
            Approx(c.en()));
    REQUIRE(json_file["wavefunctions"]["core"][c.shortSymbol()]["kappa"] ==
            Approx(c.kappa()));
  }

  // More complicated: check the read-in wavefunctions and grid match
  const std::vector<double> r = json_file["radial"]["r"];
  const std::vector<double> dr = json_file["radial"]["dr"];
  for (const auto &v : wf.valence()) {
    REQUIRE(json_file["wavefunctions"]["valence"][v.shortSymbol()]["en"] ==
            Approx(v.en()));

    const std::vector<double> f =
      json_file["wavefunctions"]["valence"][v.shortSymbol()]["f"];
    const std::vector<double> g =
      json_file["wavefunctions"]["valence"][v.shortSymbol()]["g"];

    // <r> expectation value:
    const auto r_ev_wf = v * (wf.grid().r() * v);

    using namespace qip::overloads;
    const auto integrand_json = (f * r * f + g * r * g) * dr;
    const auto r_ev_json = NumCalc::integrate(1.0, 0, 0, integrand_json);
    REQUIRE(r_ev_json == Approx(r_ev_wf));
  }
}

//============================================================================
TEST_CASE("ampsci - Integration Tests",
          "[ampsci][HartreeFock][Correlations][integration][slow]") {

  // This is just a test case that failed to converge before an update
  std::cout << "Accuracy test for Cs\n";
  std::cout
    << "Note: to make quick, doesn't do fullest numerics - so is not 100% "
       "best possible\n";
  std::cout << "Does not include Breit/QED, structure radiation etc.\n";

  std::string run_label = "deleteme_" + qip::random_string(5);

  const auto input_string = std::string{R"(
    Atom {
      Z = Cs;
      run_label = )" + run_label + R"(;
    }
    HartreeFock {
      core = [Xe];
      valence = 6sp5d4f;
    }
    Grid {
      r0 = 1.0e-6;
      rmax = 120.0;
      num_points = 1600;
    }
    Basis {
      number = 30;
      order = 7;
      rmax = 35.0;
      states = 30spdfghi;
    }
    Correlations {
      n_min_core = 3;
      all_order = true;
      each_valence = false;
      fk = 0.7, 0.6, 0.83, 0.89, 0.95, 0.97, 0.99;
    }
  )"};

  // run ampsci
  fmt2::styled_print(fg(fmt::color::blue), "\nRun AMPSCI:\n");
  const auto wf = ampsci(IO::InputBlock{"", input_string});

  // --------- Energies -----------

  // Expt. energies:
  std::map<std::string, double> nist{{"6s+", -31406.4677}, {"6p-", -20228.1996},
                                     {"6p+", -19674.1606}, {"5d-", -16907.2109},
                                     {"5d+", -16809.6254}, {"4f-", -6934.2408},
                                     {"4f+", -6934.4222}};

  fmt2::styled_print(fg(fmt::color::blue), "\nCompare energies:\n");
  fmt::print("{:3} {:9} [{:9}] {:5}\n", "v", "En", "NIST", "Eps(%)");
  for (const auto &v : wf.valence()) {
    // 0.1% for s, p; 0.2% for f; 1% for d
    const auto eps_target = v.l() == 2 ? 1.0e-2 : v.l() == 3 ? 2.0e-3 : 1.0e-3;

    const auto nist_val = nist[v.shortSymbol()];
    const auto ampsci_val = v.en() * PhysConst::Hartree_invcm;
    const auto eps = (ampsci_val / nist_val - 1.0);

    fmt::print("{:3} {:9.2f} [{:9.2f}] {:+4.2f}%\n", v.shortSymbol(),
               ampsci_val, nist_val, eps * 100.0);

    CHECK(std::abs(eps) < eps_target);
  }

  // --------- E1 -----------

  std::map<std::string, double> E1_exp{{"6s+6p-", 4.5057}, {"6s+6p+", 6.3398}};

  fmt2::styled_print(fg(fmt::color::blue),
                     "\nCompare E1 matrix elements (No SR):\n");
  const auto e1 = DiracOperator::generate("E1", {}, wf);
  auto dVe1 = ExternalField::TDHF(e1.get(), wf.vHF());
  dVe1.solve_core(0.0); // ignore frequency dependence ?
  fmt::print("\n{:3} {:3}  {:7}  [{:7}]  {:5}\n", "v", "w", "E1", "Expt",
             "Eps(%)");
  for (const auto &v : wf.valence()) {
    for (const auto &w : wf.valence()) {

      auto exp_val = E1_exp[v.shortSymbol() + w.shortSymbol()];
      if (exp_val == 0.0)
        continue;

      const auto ampsci_val = std::abs(e1->reducedME(v, w) + dVe1.dV(v, w));
      const auto eps = (ampsci_val / exp_val - 1.0);

      fmt::print("{:3}-{:3}  {:7.4f}  [{:7.4f}]  {:+4.2f}%\n", v.shortSymbol(),
                 w.shortSymbol(), ampsci_val, exp_val, eps * 100.0);

      CHECK(std::abs(eps) < 0.3e-2); // 0.3%
    }
  }

  // --------- HFS -----------

  fmt2::styled_print(fg(fmt::color::blue),
                     "\nCompare Hyperfine A constants:\n");
  std::map<std::string, double> A_exp{
    {"6s+", 2298.1579425}, {"6p-", 291.9135}, {"6p+", 50.28163},
    // {"5d-", 48.78},
    // {"5d+", -21.24}
  };

  const auto hfs = DiracOperator::generate(
    "hfs", {"hfs_options", "k=1; F=SingleParticle;"}, wf);
  auto dVhfs =
    ExternalField::DiagramRPA(hfs.get(), wf.basis(), wf.vHF(), wf.identity());

  // auto sr = MBPT::StructureRad(wf.basis(), wf.FermiLevel(), {4, 20});

  // fmt::print("{:3} {:7} + {:7} = {:7}  [{:7}] {:5}\n", "v", "A0", "SR+N",
  //  "Final", "Expt.", "Eps(%)");
  fmt::print("{:3} {:7}  [{:7}] {:5}\n", "v", "A0", "SR+N", "Final", "Expt.",
             "Eps(%)");
  for (const auto &v : wf.valence()) {
    auto exp_val = A_exp[v.shortSymbol()];
    if (exp_val == 0.0)
      continue;

    // 7% with no SR (need large SR basis to fix high-l states)
    const auto eps_target = 7.0e-2;

    const auto f =
      hfs->matel_factor(DiracOperator::MatrixElementType::HFConstant, v, v);
    const auto ampsci_val_0 = f * (hfs->reducedME(v, v) + dVhfs.dV(v, v));

    // const auto [d, d2] = sr.srn(hfs.get(), v, v, 0.0);

    // const auto srn = d * f;

    const auto ampsci_val = ampsci_val_0; // + srn;

    const auto eps = (ampsci_val / exp_val - 1.0);

    // fmt::print("{:3} {:7.2f} + {:7.2f} = {:7.2f}  [{:7.2f}] {:+4.2f}%\n",
    //            v.shortSymbol(), ampsci_val_0, 0.0, ampsci_val, exp_val,
    //            eps * 100.0);
    fmt::print("{:3} {:7.2f}  [{:7.2f}] {:+4.2f}%\n", v.shortSymbol(),
               ampsci_val_0, exp_val, eps * 100.0);

    CHECK(std::abs(eps) < eps_target);
  }
}