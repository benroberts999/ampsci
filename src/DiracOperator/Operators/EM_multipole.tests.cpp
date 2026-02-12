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
  const std::vector<std::string> op_names{
      "VE_Len", "VE", "VM", "VL", "VT", "AE", "AM", "AL", "AT", "S", "P"};

  auto use_helper_function =
      [&](const std::string &name, int k,
          SphericalBessel::JL_table *jl =
              nullptr) -> std::unique_ptr<DiracOperator::TensorOperator> {
    using namespace DiracOperator;

    // Vector
    if (name == "VE")
      return std::make_unique<VEk>(wf.grid(), k, 1.0e-4, jl);
    if (name == "VM")
      return std::make_unique<VMk>(wf.grid(), k, 1.0e-4, jl);
    if (name == "VL")
      return std::make_unique<VLk>(wf.grid(), k, 1.0e-4, jl);
    if (name == "VT")
      return std::make_unique<Phik>(wf.grid(), k, 1.0e-4, jl);
    if (name == "VE_Len")
      return std::make_unique<VEk_Len>(wf.grid(), k, 1.0e-6, jl);

    // Axial
    if (name == "AE")
      return std::make_unique<AEk>(wf.grid(), k, 1.0e-4, jl);
    if (name == "AM")
      return std::make_unique<AMk>(wf.grid(), k, 1.0e-4, jl);
    if (name == "AL")
      return std::make_unique<ALk>(wf.grid(), k, 1.0e-4, jl);
    if (name == "AT")
      return std::make_unique<Phi5k>(wf.grid(), k, 1.0e-4, jl);

    // Scalar / pseudoscalar
    if (name == "S")
      return std::make_unique<Sk>(wf.grid(), k, 1.0e-4, jl);
    if (name == "P")
      return std::make_unique<S5k>(wf.grid(), k, 1.0e-4, jl);

    std::cout << "Error: unknown operator name '" << name << "'\n";
    return std::make_unique<NullOperator>();
  };

  // 1) For each operator, check that radial_rhs reduced product equals radialIntegral
  for (const auto &name : op_names) {
    std::cout << name << "\n";
    for (int k = 0; k <= max_L; ++k) {
      for (double omega : omegas) {

        // Version one: using 'generate'
        std::string opts = std::string("k=") + std::to_string(k) +
                           "; omega=" + std::to_string(omega) + "; ";

        const auto op_type = std::string(1, name.at(0));
        const auto component =
            name.size() > 1 ? std::string(1, name.at(1)) : "";
        const auto form = name.size() > 3 ? std::string(1, name.at(3)) : "";

        opts += "type=" + op_type + ";";
        opts += "component=" + component + ";";
        opts += "form=" + form + ";";

        auto h = DiracOperator::generate("Multipole", {"", opts}, wf);

        //using "helper" functions
        auto h_3 = use_helper_function(name, k);
        h_3->updateFrequency(omega);

        // Using Bessel loolup table
        auto h_4 = use_helper_function(name, k, &jl_table);
        h_4->updateFrequency(omega);

        for (const auto &a : orbs) {
          for (const auto &b : orbs) {

            const auto rme = h->reducedME(a, b);
            const auto rad_int = h->radialIntegral(a, b);
            const auto C_ang = h->angularF(a.kappa(), b.kappa());
            const auto rme_2 = a * h->reduced_rhs(a.kappa(), b);

            const auto rme_3 = h_3->reducedME(a, b);
            const auto rme_4 = h_4->reducedME(a, b);

            if (h->isZero(a, b)) {
              REQUIRE(rme == 0.0);
              REQUIRE(rad_int == 0.0);
              continue;
            }

            // if (!h->isZero(a, b) && k > 0) {
            // REQUIRE(rme != 0.0);
            // No, sometimes zero if (ka-kb)=0 or similar..
            // }

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

  // small qr limit: electric and magnetic only
  // Add others when we can
  for (const auto omega : {1.0e-6, 1.0e-5}) {
    DiracOperator::E1 E1(wf.grid());
    DiracOperator::E1v E1v(PhysConst::alpha, omega);
    DiracOperator::M1 M1(wf.grid(), PhysConst::alpha, omega);

    // converts "transition" operators to "moment" form
    const auto ff = DiracOperator::multipole::moment_factor(1, omega);

    auto E1_w = DiracOperator::VEk_Len(wf.grid(), 1, 0.0, nullptr);
    E1_w.updateFrequency(omega);

    auto E1v_w = DiracOperator::VEk(wf.grid(), 1, 0.0, nullptr);
    E1v_w.updateFrequency(omega);

    auto M1_w = DiracOperator::VMk(wf.grid(), 1, 0.0, nullptr);
    M1_w.updateFrequency(omega);

    const auto eps = omega * 1.0e-3;

    for (const auto &a : orbs) {
      for (const auto &b : orbs) {
        // Compare reduced matrix elements with table-backed versions

        // nb: extra sign, because E1 defined as -|e|r, rather than |e|r
        REQUIRE(E1.reducedME(a, b) ==
                Approx(-ff * E1_w.reducedME(a, b)).epsilon(eps));

        REQUIRE(E1v.reducedME(a, b) ==
                Approx(-ff * E1v_w.reducedME(a, b)).epsilon(eps));

        REQUIRE(M1.reducedME(a, b) * PhysConst::muB_CGS ==
                Approx(ff * M1_w.reducedME(a, b)).epsilon(eps));
      }
    }
  }
}