#include "CI_Integrals.hpp"
#include "ConfigurationInteraction.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Random.hpp"

//==============================================================================
TEST_CASE("CI: Configuration Interaction unit tests", "[CI][unit]") {

  std::cout << "CI, unit tests (not meant to be accurate)\n";

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"He", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[]", std::nullopt, 1.0e-5);
  wf.formBasis(SplineBasis::Parameters("6spd", 20, 6, 1.0e-2, 1.0e-2, 20.0));

  std::string qk_filename = "deleteme_" + qip::random_string(3) + ".qk.abf";
  std::string ci_filename = "deleteme_" + qip::random_string(3) + ".ci.abf";

  const IO::InputBlock input{"CI", "ci_basis = 6spd;"
                                   "n_min_core = 1;"
                                   "J + = 0, 1;"
                                   "J - = 0, 1;"
                                   "qk_file = " +
                                     qk_filename +
                                     ";"
                                     "num_solutions = 2;"
                                     "ci_file = " +
                                     ci_filename + ";"};

  // Initial CI calculation (writes to file):
  const auto CIWFs = CI::configuration_interaction(input, wf);

  // Expected values: (not accurate, just simple unit test)
  std::vector J{0, 1, 0, 1};
  std::vector parity{1, 1, -1, -1};
  std::vector energy{std::vector{-2.84834207, -2.13740814},
                     {-2.17164001, -2.0657259},
                     {-2.13100411, -2.05532175},
                     {-2.13100096, -2.12092003}};
  std::vector gj{
    std::vector{0.0, 0.0}, {1.9999, 1.9999}, {0.0, 0.0}, {1.5, 1.0}};

  for (std::size_t i = 0; i < CIWFs.size(); ++i) {
    const auto &ci_wf = CIWFs.at(i);
    REQUIRE(ci_wf.twoJ() == 2 * J[i]);
    REQUIRE(ci_wf.parity() == parity[i]);
    for (std::size_t j = 0ul; j < ci_wf.num_solutions(); ++j) {
      REQUIRE(ci_wf.energy(j) == Approx(energy[i][j]).epsilon(1.0e-2));
      REQUIRE(ci_wf.info(j).gJ == Approx(gj[i][j]).epsilon(1.0e-2));
    }
  }

  //-----------------------------------------------------------------------
  std::cout << "\nRe-run CI, reading in from file:\n";

  // Re-"run" identical CI calculation
  // This time: should read in prev solution from disk
  // Read from ci file; compare energies, gJ, and coefficients
  auto input2 = input;
  input2.add("read_only = true;");
  const auto CIWFs_2 = CI::configuration_interaction(input2, wf);

  REQUIRE(CIWFs.size() == CIWFs_2.size());

  for (std::size_t i = 0; i < CIWFs.size(); ++i) {
    const auto &w = CIWFs.at(i);
    const auto &r = CIWFs_2.at(i);
    REQUIRE(w.twoJ() == r.twoJ());
    REQUIRE(w.parity() == r.parity());
    REQUIRE(w.num_solutions() == r.num_solutions());
    REQUIRE(w.CSFs().size() == r.CSFs().size());
    for (std::size_t s = 0; s < w.num_solutions(); ++s) {
      REQUIRE(w.energy(s) == Approx(r.energy(s)));
      REQUIRE(w.info(s).gJ == Approx(r.info(s).gJ));
      const auto wc = w.coefs(s);
      const auto rc = r.coefs(s);
      REQUIRE(wc.size() == rc.size());
      REQUIRE(wc.size() == w.CSFs().size());
      for (std::size_t k = 0; k < wc.size(); ++k) {
        REQUIRE(wc[k] == Approx(rc[k]));
      }
    }
  }

  //-----------------------------------------------------------------------
  std::cout << "\nRe-run CI, reading in from file: but ask for more states\n";
  // Now, request 3 solutions:
  // ci file has 2, so must re-solve and overwrite
  const IO::InputBlock input3{"CI", "ci_basis = 6spd;"
                                    "n_min_core = 1;"
                                    "J + = 0, 1;"
                                    "J - = 0, 1;"
                                    "qk_file = " +
                                      qk_filename +
                                      ";"
                                      "num_solutions = 3;"
                                      "print_details = false;"
                                      "ci_file = " +
                                      ci_filename + ";"};
  // don't read in
  const IO::InputBlock input4{"CI", "ci_basis = 6spd;"
                                    "n_min_core = 1;"
                                    "J + = 0, 1;"
                                    "J - = 0, 1;"
                                    "qk_file = false;"
                                    "num_solutions = 3;"
                                    "print_details = false;"
                                    "ci_file = " +
                                      ci_filename + ";"};

  // Reads from file (2 stored, 3 required, re-calculates)
  const auto CIWFs_3a = CI::configuration_interaction(input3, wf);
  // calculate 3 from scratch
  const auto CIWFs_3b = CI::configuration_interaction(input3, wf);

  REQUIRE(CIWFs_3a.size() == CIWFs.size());
  REQUIRE(CIWFs_3b.size() == CIWFs.size());
  for (std::size_t i = 0; i < CIWFs_3a.size(); ++i) {
    const auto &a = CIWFs_3a.at(i);
    const auto &b = CIWFs_3b.at(i);
    REQUIRE(a.num_solutions() == 3);
    REQUIRE(b.num_solutions() == 3);
    // First 2 solutions match CIWFs
    for (std::size_t s = 0; s < 2; ++s) {
      REQUIRE(a.energy(s) == Approx(CIWFs.at(i).energy(s)));
      REQUIRE(a.info(s).gJ == Approx(CIWFs.at(i).info(s).gJ));
    }
    // All 3 identical between solve and read
    for (std::size_t s = 0; s < 3; ++s) {
      REQUIRE(a.energy(s) == Approx(b.energy(s)));
      REQUIRE(a.info(s).gJ == Approx(b.info(s).gJ));
      const auto ac = a.coefs(s);
      const auto bc = b.coefs(s);
      REQUIRE(ac.size() == bc.size());
      for (std::size_t k = 0; k < ac.size(); ++k) {
        REQUIRE(ac[k] == Approx(bc[k]));
      }
    }
  }

  //-----------------------------------------------------------------------
  // basic/Misc tests

  // Term(int two_J, int L, int two_S, int parity)

  REQUIRE(CI::Term_Symbol(2, 3, 2, +1) == "3^F_1");
  REQUIRE(CI::Term_Symbol(2, 3, 2, -1) == "3^F°_1");
  REQUIRE(CI::Term_Symbol(1, 3, 2, -1) == "3^F°_1/2");
  REQUIRE(CI::Term_Symbol(1, 2, 1, -1) == "2^D°_1/2");
  REQUIRE(CI::Term_Symbol(6, 0, 0, 1) == "1^S_3");
}
