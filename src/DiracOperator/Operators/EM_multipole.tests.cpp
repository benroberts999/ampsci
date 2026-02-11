#include "DiracOperator/GenerateOperator.hpp"
#include "Maths/Grid.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "include.hpp"
#include "qip/Vector.hpp"
#include <iostream>
#include <string>
#include <vector>

TEST_CASE("EM_multipole operators", "[DiracOperator][unit][EM_multipole][jL]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "EM_multipole\n";

  // Construct a wavefunction similar to other tests
  Wavefunction wf({100, 1.0e-4, 50.0, 1.0, "loglinear"}, {"Cs", 133, "Fermi"});

  // Populate simple valence orbitals (H-like) so getState(...) works
  auto &orbs = wf.valence();
  for (int ik = 0; ik <= Angular::indexFromKappa(7); ik++) {
    int kappa = Angular::kappaFromIndex(ik);
    int n = Angular::l_k(kappa) + 1;
    orbs.push_back(DiracSpinor::exactHlike(n, kappa, wf.grid_sptr(), 1.0));
  }

  const auto &rvec = wf.grid().r();

  // We will test these omegas. To ensure the JL_table lookups are exact we
  // build the JL_table q-grid from the *exact* q=alpha*omega values we will
  // request (plus a very small omega used in the small-q approximation).
  const auto omega0 = 1.0e-6;
  const std::vector<double> omegas{omega0, 1.0e-2, 1.0, 10.0};

  using namespace qip::overloads;
  std::vector<double> qvec = omegas * PhysConst::alpha;

  const int max_L = 3;
  SphericalBessel::JL_table jl_table(max_L + 1, qvec, rvec);

  // Operators to test: name and whether VEk_Len should be skipped for gamma5
  const std::vector<std::string> op_names{"VEk",  "VMk", "VLk",
                                          "Vk_w", "Sk",  "VEk_Len"};

  // lambda to generate operator using "helper" functions
  auto use_helper_function =
      [&](const std::string &name, int k, bool gamma5,
          SphericalBessel::JL_table *jl =
              nullptr) -> std::unique_ptr<DiracOperator::TensorOperator> {
    using namespace DiracOperator;
    // Use multipole helpers when they map naturally to the operator.
    if (name == "VEk")
      return multipole::V_sigma_K(wf.grid(), +1, k, gamma5, jl);
    if (name == "VMk")
      return multipole::V_sigma_K(wf.grid(), 0, k, gamma5, jl);
    if (name == "VLk")
      return multipole::V_sigma_K(wf.grid(), -1, k, gamma5, jl);
    if (name == "Vk_w")
      return multipole::Phi_K(wf.grid(), k, gamma5, jl);
    if (name == "Sk")
      return multipole::S_K(wf.grid(), k, gamma5, jl);
    if (name == "VEk_Len" && !gamma5)
      return std::make_unique<DiracOperator::VEk_Len>(wf.grid(), k, 1.0e-6, jl);
    return std::unique_ptr<DiracOperator::TensorOperator>(nullptr);
  };

  // 1) For each operator, check that radial_rhs reduced product equals radialIntegral
  for (const auto &name : op_names) {
    std::cout << name << "\n";
    for (int k = 0; k <= max_L; ++k) {
      for (double omega : omegas) {
        for (bool gamma5 : {false, true}) {
          if (name == "VEk_Len" && gamma5)
            continue;

          // Version one: using 'generate'
          std::string opts = std::string("k=") + std::to_string(k) +
                             "; omega=" + std::to_string(omega) + "; ";
          if (gamma5)
            opts += std::string("gamma5=true;");
          auto h = DiracOperator::generate(name, {"", opts}, wf);

          //using "helper" functions
          auto h_3 = use_helper_function(name, k, gamma5);
          h_3->updateFrequency(omega);

          // Using Bessel loolup table
          auto h_4 = use_helper_function(name, k, gamma5, &jl_table);
          h_4->updateFrequency(omega);

          for (const auto &a : orbs) {
            for (const auto &b : orbs) {

              const auto rme = h->reducedME(a, b);
              const auto rad_int = h->radialIntegral(a, b);
              const auto C_ang = h->angularF(a.kappa(), b.kappa());
              const auto rme_2 = a * h->reduced_rhs(a.kappa(), b);

              const auto rme_3 = h_3->reducedME(a, b);
              const auto rme_4 = h_4->reducedME(a, b);

              // 1. test redME vs rad
              REQUIRE(rad_int * C_ang == Approx(rme));

              // 2. test red_lhs vs. redME
              if (rme != 0.0)
                REQUIRE(rme_2 == Approx(rme));
              else
                REQUIRE(rme_2 == Approx(rme).margin(1.0e-17));

              // 3. test 'directly constructed version (helper function)
              // (includes test of update omega)
              REQUIRE(rme_3 == Approx(rme));

              // 4. Test the Bessel table version
              REQUIRE(rme_4 == Approx(rme));
            }
          }
        }
      }
    }
  }

  // small qr limit: electric and magnetic only
  // Add others when we can
  for (const auto omega : {1.0e-6, 1.0e-5}) {
    DiracOperator::E1 E1(wf.grid());
    DiracOperator::E1v E1v(PhysConst::alpha, omega);
    DiracOperator::M1 M1(wf.grid(), PhysConst::alpha, omega);

    // converts "transition" operators to "moment" form
    const auto ff = DiracOperator::multipole::moment_factor(1, omega);

    auto E1_w = use_helper_function("VEk_Len", 1, false);
    E1_w->updateFrequency(omega);

    auto E1v_w = use_helper_function("VEk", 1, false);
    E1v_w->updateFrequency(omega);

    auto M1_w = use_helper_function("VMk", 1, false);
    M1_w->updateFrequency(omega);

    const auto eps = omega * 1.0e-3;

    for (const auto &a : orbs) {
      for (const auto &b : orbs) {
        // Compare reduced matrix elements with table-backed versions
        REQUIRE(E1.reducedME(a, b) ==
                Approx(ff * E1_w->reducedME(a, b)).epsilon(eps));

        REQUIRE(E1v.reducedME(a, b) ==
                Approx(ff * E1v_w->reducedME(a, b)).epsilon(eps));

        // negative? factored i out?
        REQUIRE(-1 * M1.reducedME(a, b) * PhysConst::muB_CGS ==
                Approx(ff * M1_w->reducedME(a, b)).epsilon(eps));
      }
    }
  }
}