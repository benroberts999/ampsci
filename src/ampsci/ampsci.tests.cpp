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
  std::cout << "Simple 'unit' test with small parameters only. Not meant to "
               "be accurate\n";

  std::string run_label = "deleteme_" + qip::random_string(5);

  IO::InputBlock input{""};
  {
    IO::InputBlock atom{"Atom"};
    atom.set("Z", std::string{"Cs"});
    atom.set("run_label", run_label);
    atom.set("json_out", true);
    input.set_block("Atom", atom);

    IO::InputBlock hf{"HartreeFock"};
    hf.set("core", std::string{"[Xe]"});
    hf.set("valence", std::string{"6sp"});
    input.set_block("HartreeFock", hf);

    IO::InputBlock grid{"Grid"};
    grid.set("r0", 1.0e-5);
    grid.set("rmax", 75.0);
    grid.set("num_points", 800);
    input.set_block("Grid", grid);

    IO::InputBlock basis{"Basis"};
    basis.set("number", 30);
    basis.set("order", 7);
    basis.set("rmax", 35.0);
    basis.set("states", std::string{"12spdf"});
    input.set_block("Basis", basis);

    IO::InputBlock corr{"Correlations"};
    corr.set("stride", 10);
    corr.set("n_min_core", 4);
    corr.set("each_valence", false);
    corr.set("read", false);
    input.set_block("Correlations", corr);

    IO::InputBlock me_hfs{"MatrixElements"};
    me_hfs.set("operator", std::string{"hfs"});
    me_hfs.set("rpa", std::string{"diagram"});
    me_hfs.set("off-diagonal", false);
    input.add_module(me_hfs);

    IO::InputBlock me_e1{"MatrixElements"};
    me_e1.set("operator", std::string{"E1"});
    me_e1.set("rpa", std::string{"TDHF"});
    me_e1.set("omega", 0.0);
    input.add_module(me_e1);
  }

  const auto wf = ampsci(input);

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
TEST_CASE("ampsci - basic Feynman unit test",
          "[ampsci][Feynman][Correlations][unit]") {

  // This is just a test case that failed to converge before an update
  std::cout << "Simple 'unit' test with small parameters only. Not meant to "
               "be accurate\n";

  std::string run_label = "deleteme_" + qip::random_string(5);

  IO::InputBlock input{""};
  {
    IO::InputBlock atom{"Atom"};
    atom.set("Z", std::string{"Cs"});
    atom.set("run_label", run_label);
    input.set_block("Atom", atom);

    IO::InputBlock hf{"HartreeFock"};
    hf.set("core", std::string{"[Xe]"});
    hf.set("valence", std::string{"6sp"});
    hf.set("QED", true);
    input.set_block("HartreeFock", hf);

    IO::InputBlock nuc{"Nucleus"};
    nuc.set("type", std::string{"Gaussian"});
    input.set_block("Nucleus", nuc);

    IO::InputBlock grid{"Grid"};
    grid.set("r0", 1.0e-5);
    grid.set("rmax", 80.0);
    grid.set("num_points", 1000);
    input.set_block("Grid", grid);

    IO::InputBlock basis{"Basis"};
    basis.set("number", 35);
    basis.set("order", 7);
    basis.set("rmax", 40.0);
    basis.set("states", std::string{"25spdfg"});
    basis.set("type", std::string{"Johnson"});
    input.set_block("Basis", basis);

    IO::InputBlock spec{"Spectrum"};
    spec.set("number", 70);
    spec.set("order", 7);
    spec.set("rmax", 70.0);
    spec.set("states", std::string{"6spd"});
    input.set_block("Spectrum", spec);

    IO::InputBlock corr{"Correlations"};
    corr.set("n_min_core", 5);
    corr.set("all_order", true);
    corr.set("stride", 14);
    corr.set("each_valence", false);
    corr.set("imag_omega", std::vector<double>{0.02, 2.0});
    corr.set("read", false);
    input.set_block("Correlations", corr);
  }

  const auto wf = ampsci(input);

  CHECK(wf.Znuc() == 55);
  CHECK(wf.identity() == "Cs1q_" + run_label);

  // const auto &spectrum = wf.sp

  // require all valence states converge
  for (const auto &v : wf.valence()) {
    REQUIRE(v.eps() < 1.0e-9);
  }

  std::cout << "\nCompare valence and spectrum:\n";
  fmt::print("{:3} {:8} {:8} {:6} {:6}\n", "v", "E_val", "E_spec", "eps",
             "1-<v|s>");
  for (const auto &v : wf.spectrum()) {
    const auto p1 = wf.getState(v.n(), v.kappa());
    const auto p2 = wf.getState(v.shortSymbol());
    REQUIRE(p1 == p2);

    // ensure is valence state:
    if (wf.isInCore(v.n(), v.kappa())) {
      REQUIRE(p1 != nullptr);
      continue;
    }

    if (p1 == nullptr)
      continue;

    REQUIRE(wf.isInValence(v.n(), v.kappa()));

    REQUIRE(v.n() == p1->n());
    REQUIRE(v.kappa() == p1->kappa());

    const auto eps = (p1->en() - v.en()) / (p1->en() + v.en());
    const auto ort = 1.0 - v * (*p1);

    fmt::print("{:3} {:8.5f} {:8.5f} {:6.0e} {:6.0e}\n", v.shortSymbol(),
               p1->en(), v.en(), eps, ort);

    CHECK(v.en() == Approx(p1->en()).epsilon(1.0e-3));
    CHECK(ort == Approx(0.0).margin(1.0e-6));
  }
}

//============================================================================
TEST_CASE("ampsci - basic CI unit test", "[ampsci][CI][unit]") {

  // This is just a test case that failed to converge before an update
  std::cout << "Simple 'unit' test with small parameters only. Not meant to "
               "be accurate\n";

  std::string run_label = "deleteme_" + qip::random_string(5);

  IO::InputBlock input1{""};
  {
    IO::InputBlock atom{"Atom"};
    atom.set("Z", std::string{"Be"});
    atom.set("run_label", run_label);
    input1.set_block("Atom", atom);

    IO::InputBlock hf{"HartreeFock"};
    hf.set("core", std::string{"[He]"});
    input1.set_block("HartreeFock", hf);

    IO::InputBlock grid{"Grid"};
    grid.set("r0", 1.0e-4);
    grid.set("rmax", 50.0);
    grid.set("num_points", 700);
    input1.set_block("Grid", grid);

    IO::InputBlock basis{"Basis"};
    basis.set("number", 25);
    basis.set("order", 7);
    basis.set("rmax", 40.0);
    basis.set("states", std::string{"16spdf"});
    input1.set_block("Basis", basis);

    IO::InputBlock ci{"CI"};
    ci.set("ci_basis", std::string{"12spdf"});
    ci.set("J+", std::vector<int>{0, 1});
    ci.set("J-", std::vector<int>{0, 1});
    ci.set("num_solutions", 2);
    ci.set("sigma1", true);
    ci.set("sigma2", true);
    ci.set("cis2_basis", std::string{"3spd"});
    input1.set_block("CI", ci);

    IO::InputBlock ci_me{"CI_matrixElements"};
    ci_me.set("operator", std::string{"E1"});
    ci_me.set("rpa", false);
    ci_me.set("ci_basis", std::string{"12spdf"});
    ci_me.set("omega", 0.0);
    ci_me.set("J", std::vector<int>{0, 1});
    input1.add_module(ci_me);
  }
  const auto wf = ampsci(input1);

  // Expt. for lowest few states:
  const auto E0 = -146882.86 - 75192.64;
  std::map<std::string, double> nist{{"0+", E0},
                                     {"0-", E0 + 21978.31},
                                     {"1+", E0 + 52080.943},
                                     {"1-", E0 + 21978.925}};

  REQUIRE(wf.Znuc() == 4);
  REQUIRE(wf.identity() == "Be2_" + run_label);

  std::cout << "\nCompare CI energies (rough - very small basis)\n";
  fmt::print("{:3} {:4} {:7} [{:7}] {:6}\n", "Jpi", "conf", "En", "Expt.",
             "eps");
  const auto &ci_wfs = wf.CIwfs();
  for (auto &ci_wf : ci_wfs) {
    const auto twoj = ci_wf.twoJ();
    const auto parity = ci_wf.parity();

    const auto ci_wf2 = wf.CIwf(twoj / 2, parity);
    REQUIRE(ci_wf2 == &ci_wf);

    const auto ci_en = ci_wf.energy(0) * PhysConst::Hartree_invcm;

    const auto key = std::to_string(twoj / 2) + (parity == 1 ? "+" : "-");
    const auto E_exp = nist[key];

    const auto eps = ci_en / E_exp - 1.0;

    fmt::print("{:3} {:4} {:7.0f} [{:7.0f}] {:6.0e}\n", key,
               ci_wf.info(0).config, ci_en, E_exp, eps);
  }

  // run again, by reading in qk files
  // Input almost the same, but we add 'no_new_integrals'
  IO::InputBlock input_v2{""};
  {
    IO::InputBlock atom{"Atom"};
    atom.set("Z", std::string{"Be"});
    atom.set("run_label", run_label);
    input_v2.set_block("Atom", atom);

    IO::InputBlock hf{"HartreeFock"};
    hf.set("core", std::string{"[He]"});
    input_v2.set_block("HartreeFock", hf);

    IO::InputBlock grid{"Grid"};
    grid.set("r0", 1.0e-4);
    grid.set("rmax", 50.0);
    grid.set("num_points", 700);
    input_v2.set_block("Grid", grid);

    IO::InputBlock basis{"Basis"};
    basis.set("number", 25);
    basis.set("order", 7);
    basis.set("rmax", 40.0);
    basis.set("states", std::string{"16spdf"});
    input_v2.set_block("Basis", basis);

    IO::InputBlock ci{"CI"};
    ci.set("no_new_integrals", true);
    ci.set("sort_output", true);
    ci.set("ci_basis", std::string{"12spdf"});
    ci.set("J+", std::vector<int>{0, 1});
    ci.set("J-", std::vector<int>{0, 1});
    ci.set("num_solutions", 2);
    ci.set("sigma1", true);
    ci.set("sigma2", true);
    ci.set("cis2_basis", std::string{"3spd"});
    ci.set("print_details", false);
    input_v2.set_block("CI", ci);
  }

  const auto wf2 = ampsci(input_v2);

  // loop through original wavefunction:
  for (auto &ci_wf : ci_wfs) {
    const auto twoj = ci_wf.twoJ();
    const auto parity = ci_wf.parity();

    // Find corresponding CI solution from second:
    const auto ci_wf2 = wf2.CIwf(twoj / 2, parity);
    REQUIRE(ci_wf2 != nullptr);

    for (std::size_t i = 0; i < ci_wf.num_solutions(); ++i) {

      const auto E1 = ci_wf.energy(i) * PhysConst::Hartree_invcm;
      const auto E2 = ci_wf2->energy(i) * PhysConst::Hartree_invcm;

      const auto key = std::to_string(twoj / 2) + (parity == 1 ? "+" : "-");

      const auto eps = E1 / E2 - 1.0;

      fmt::print("{:3} {:1} {:4} {:7.0f} [{:7.0f}] {:6.0e}\n", key, i,
                 ci_wf.info(i).config, E1, E2, eps);

      REQUIRE(ci_wf.info(i).gJ == Approx(ci_wf2->info(i).gJ));
      REQUIRE(E1 == Approx(E2));
      REQUIRE(ci_wf.info(i).config == ci_wf2->info(i).config);
    }
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

  std::string run_label = "deleteme_" + qip::random_string(3);

  IO::InputBlock input{""};
  {
    IO::InputBlock atom{"Atom"};
    atom.set("Z", std::string{"Cs"});
    atom.set("run_label", run_label);
    input.set_block("Atom", atom);

    IO::InputBlock hf{"HartreeFock"};
    hf.set("core", std::string{"[Xe]"});
    hf.set("valence", std::string{"6sp5d4f"});
    input.set_block("HartreeFock", hf);

    IO::InputBlock grid{"Grid"};
    grid.set("r0", 1.0e-6);
    grid.set("rmax", 120.0);
    grid.set("num_points", 1600);
    input.set_block("Grid", grid);

    IO::InputBlock basis{"Basis"};
    basis.set("number", 30);
    basis.set("order", 7);
    basis.set("rmax", 35.0);
    basis.set("states", std::string{"30spdfghi"});
    input.set_block("Basis", basis);

    IO::InputBlock corr{"Correlations"};
    corr.set("n_min_core", 3);
    corr.set("all_order", true);
    corr.set("each_valence", false);
    corr.set("fk", std::vector<double>{0.7, 0.6, 0.83, 0.89, 0.95, 0.97, 0.99});
    input.set_block("Correlations", corr);
  }

  // run ampsci
  fmt2::styled_print(fg(fmt::color::blue), "\nRun AMPSCI:\n");
  const auto wf = ampsci(input);

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
  const auto e1 = DiracOperator::generate("E1", IO::InputBlock{}, wf);
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
    {"6s+", 2298.1579425}, {"6p-", 291.9135}, {"6p+", 50.28163}};

  IO::InputBlock hfs_opts{"hfs_options"};
  hfs_opts.set("k", 1);
  hfs_opts.set("F", std::string{"SingleParticle"});
  const auto hfs = DiracOperator::generate("hfs", hfs_opts, wf);
  auto dVhfs =
    ExternalField::DiagramRPA(hfs.get(), wf.basis(), wf.vHF(), wf.identity());
  dVhfs.solve_core(0.0);

  auto sr = MBPT::StructureRad(wf.basis(), wf.FermiLevel(), {4, 18});
  sr.solve_core(hfs.get(), &dVhfs);

  fmt::print("\n{:3} {:>7}  {:>7}  [{:>7}] {:5}\n", "v", "A0", "+SR+N", "Expt.",
             "Eps(%)");
  for (const auto &v : wf.valence()) {
    auto exp_val = A_exp[v.shortSymbol()];
    if (exp_val == 0.0)
      continue;

    // 1% for s, 3% for p. (need large SR basis to fix high-l states)
    const auto eps_target = v.l() == 0 ? 1.0e-1 : 3.0e-2;

    const auto f =
      hfs->matel_factor(DiracOperator::MatrixElementType::HFConstant, v, v);
    const auto ampsci_val_0 = f * (hfs->reducedME(v, v) + dVhfs.dV(v, v));

    const auto d = sr.srn(v, v, hfs.get(), &dVhfs);

    const auto srn = d * f;

    const auto ampsci_val = ampsci_val_0 + srn;

    const auto eps = (ampsci_val / exp_val - 1.0);

    fmt::print("{:3} {:7.2f}  {:7.2f}  [{:7.2f}] {:+4.2f}%\n", v.shortSymbol(),
               ampsci_val_0, ampsci_val, exp_val, eps * 100.0);

    CHECK(std::abs(eps) < eps_target);
  }
}