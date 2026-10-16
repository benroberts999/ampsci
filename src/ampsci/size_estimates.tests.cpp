#include "ampsci/size_estimates.hpp"
#include "CI/CSF.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "ampsci/ampsci.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include "qip/Random.hpp"
#include <filesystem>
#include <string>

//==============================================================================
// The size estimates from an input file, against the tables, CI blocks, and
// ci file that the run of that input actually forms
TEST_CASE("SizeEstimates: input vs run", "[ampsci][CI][unit]") {

  std::cout << "SizeEstimates: input vs run\n";

  const std::string atom =
    "Atom { Z = Li; }"
    "HartreeFock { core = [He]; valence = 2sp; }"
    "Grid { r0 = 1.0e-5; rmax = 30.0; num_points = 800; }"
    "Basis { number = 20; order = 6; r0 = 1.0e-3; "
    "rmax = 30.0; states = 6spd; }";

  // Two electrons, with Sigma_1 (calculated here, since there is no
  // Correlations block) and Sigma_2
  {
    const auto tag = "deleteme_" + qip::random_string(4);
    const auto ci_file = tag + ".ci.abf";
    const IO::InputBlock input(
      "ampsci", atom +
                  "CI { ci_basis = 4sp; J = 0, 1; num_solutions = 2; "
                  "sigma1 = true; sigma2 = true; cis2_basis = 3sp; "
                  "print_details = false; qk_file = " +
                  tag + ".qk.abf; sk_file = " + tag +
                  ".sk.abf; ci_file = " + ci_file + "; }");

    const auto estimate = estimate_sizes(input);
    print_size_estimates(estimate);
    const auto wf = ampsci(input);
    const auto &ints = wf.CI_integrals();

    REQUIRE(estimate.qk_CI_required);
    REQUIRE(estimate.qk_MBPT_required);
    REQUIRE(!estimate.qk_all_required);
    REQUIRE(estimate.qk_CI + estimate.qk_MBPT == ints.qk.count());

    // All Q^k over the full basis (all k), as structure radiation forms
    Coulomb::QkTable qk_all;
    qk_all.fill(wf.basis(), Coulomb::YkTable(wf.basis()), -1, false);
    REQUIRE(estimate.qk_all == qk_all.count());

    // Structure radiation (as a module, or within one) requires the full table
    for (const auto &module :
         {"Module::StructureRadiation{}", "Module::StrucRad{}",
          "Module::CI_matrixElements{ StructureRadiation{} }"}) {
      const IO::InputBlock with_SR("ampsci", atom + module);
      REQUIRE(estimate_sizes(with_SR).qk_all_required);
    }

    REQUIRE(estimate.ci_blocks.size() == wf.CIwfs().size());
    double file_GB = 0.0;
    for (std::size_t i = 0; i < estimate.ci_blocks.size(); ++i) {
      const auto &block = estimate.ci_blocks[i];
      const auto &psi = wf.CIwfs()[i];
      REQUIRE(block.twoJ == psi.twoJ());
      REQUIRE(block.parity == psi.parity());
      REQUIRE(block.num_CSFs == psi.num_CSFs());
      file_GB += block.file_GB;
    }
    REQUIRE(file_GB * 1.0e9 ==
            Approx(double(std::filesystem::file_size(ci_file))));
  }

  // Three electrons, no MBPT
  {
    const auto tag = "deleteme_" + qip::random_string(4);
    const auto ci_file = tag + ".ci.abf";
    const IO::InputBlock input("ampsci",
                               atom +
                                 "CI { num_electrons = 3; ci_basis = 4sp; "
                                 "J = 1/2, 3/2; num_solutions = 2; "
                                 "print_details = false; qk_file = " +
                                 tag + ".qk.abf; ci_file = " + ci_file + "; }");

    const auto estimate = estimate_sizes(input);
    print_size_estimates(estimate);
    const auto wf = ampsci(input);

    REQUIRE(estimate.num_electrons == 3);
    // MBPT table counted, but not formed by this run
    REQUIRE(estimate.qk_MBPT > 0);
    REQUIRE(!estimate.qk_MBPT_required);
    REQUIRE(estimate.qk_CI == wf.CI_integrals().qk.count());
    REQUIRE(estimate.ci_blocks.size() == wf.CIwfs().size());
    double file_GB = 0.0;
    for (std::size_t i = 0; i < estimate.ci_blocks.size(); ++i) {
      REQUIRE(estimate.ci_blocks[i].num_CSFs == wf.CIwfs()[i].num_CSFs());
      file_GB += estimate.ci_blocks[i].file_GB;
    }
    REQUIRE(file_GB * 1.0e9 ==
            Approx(double(std::filesystem::file_size(ci_file))));
  }
}
