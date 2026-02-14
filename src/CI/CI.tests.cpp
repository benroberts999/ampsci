#include "CI_Integrals.hpp"
#include "ConfigurationInteraction.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Random.hpp"

//==============================================================================
TEST_CASE("CI: Configuration Interaction, unit tests", "[CI][unit]") {

  std::cout << "\n----------------------------------------\n";
  std::cout << "CI, unit tests (not meant to be accurate)\n";

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"He", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[]", 1.0e-5);
  wf.formBasis(SplineBasis::Parameters("6spd", 20, 6, 1.0e-2, 1.0e-2, 20.0));

  std::string qk_filename = "deleteme_" + qip::random_string(6) + ".qk.abf";

  const IO::InputBlock input{"CI", "ci_basis = 6spd;"
                                   "n_min_core = 1;"
                                   "J + = 0, 1;"
                                   "J - = 0, 1;"
                                   "qk_file = " +
                                     qk_filename +
                                     ";"
                                     "num_solutions = 2;"};

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

  // again, this time, should read in qk file
  const auto CIWFs_2 = CI::configuration_interaction(input, wf);

  for (std::size_t i = 0; i < CIWFs.size(); ++i) {
    const auto &ci_wf = CIWFs.at(i);
    const auto &ci_wf_2 = CIWFs_2.at(i);
    REQUIRE(ci_wf.twoJ() == ci_wf_2.twoJ());
    REQUIRE(ci_wf.parity() == ci_wf_2.parity());
    for (std::size_t j = 0ul; j < ci_wf.num_solutions(); ++j) {
      REQUIRE(ci_wf.energy(j) == Approx(ci_wf_2.energy(j)));
      REQUIRE(ci_wf.info(j).gJ == Approx(ci_wf_2.info(j).gJ));
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
