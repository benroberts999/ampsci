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

  // Construct a wavefunction similar to other tests
  Wavefunction wf({75, 1.0e-4, 50.0, 1.0, "loglinear"}, {"Cs", 133, "Fermi"});

  // Populate simple valence orbitals (H-like) so getState(...) works
  auto &orbs = wf.valence();
  for (auto ik = 0ul; ik <= 7; ik++) {
    int kappa = Angular::kindex_to_kappa(ik);
    int n = Angular::l_k(kappa) + 1;
    orbs.push_back(DiracSpinor::exactHlike(n, kappa, wf.grid_sptr(), 1.0));
  }

  const auto &rvec = wf.grid().r();

  // We will test these omegas. To ensure the JL_table lookups are exact we
  // build the JL_table q-grid from the *exact* q=alpha*omega values we will
  // request (plus a very small omega used in the small-q approximation).
  const auto omega0 = 1.0e-6;
  const std::vector<double> omegas{omega0, 0.5, 10.0};

  using namespace qip::overloads;
  std::vector<double> qvec = omegas * PhysConst::alpha;

  const int max_L = 3;
  SphericalBessel::JL_table jl_table(max_L + 1, qvec, rvec);

  // Operators to test: name and whether VEk_Len should be skipped for gamma5
  const std::vector<std::string> op_names{
    "VE_Len", "VE", "VM", "VL", "VT", "AE", "AM", "AL", "AT", "S", "P"};

  auto use_helper_function =
    [&](const std::string &name, bool low_q, int k,
        SphericalBessel::JL_table *jl =
          nullptr) -> std::unique_ptr<DiracOperator::TensorOperator> {
    using namespace DiracOperator;

    if (low_q) {
      // Vector
      if (name == "VE")
        return std::make_unique<VEk_lowq>(wf.grid(), k, 1.0e-4);
      if (name == "VM")
        return std::make_unique<VMk_lowq>(wf.grid(), k, 1.0e-4);
      if (name == "VL")
        return std::make_unique<VLk_lowq>(wf.grid(), k, 1.0e-4);
      if (name == "VT")
        return std::make_unique<Phik_lowq>(wf.grid(), k, 1.0e-4);

      // Axial
      if (name == "AE")
        return std::make_unique<AEk_lowq>(wf.grid(), k, 1.0e-4);
      if (name == "AM")
        return std::make_unique<AMk_lowq>(wf.grid(), k, 1.0e-4);
      if (name == "AL")
        return std::make_unique<ALk_lowq>(wf.grid(), k, 1.0e-4);
      if (name == "AT")
        return std::make_unique<Phi5k_lowq>(wf.grid(), k, 1.0e-4);

      // Scalar / pseudoscalar
      if (name == "S")
        return std::make_unique<Sk_lowq>(wf.grid(), k, 1.0e-4);
      if (name == "P")
        return std::make_unique<S5k_lowq>(wf.grid(), k, 1.0e-4);
    }

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
    for (const auto &low_q : {false, true}) {
      for (int k = 0; k <= max_L; ++k) {
        for (double omega : omegas) {

          if (low_q && name == "VE_Len")
            continue;

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
          if (low_q) {
            opts += "low_q=true;";
          }

          auto h = DiracOperator::generate("Multipole", {"", opts}, wf);

          //using "helper" functions
          auto h_3 = use_helper_function(name, low_q, k);
          h_3->updateFrequency(omega);

          // Using Bessel loolup table
          std::unique_ptr<DiracOperator::TensorOperator> h_4{nullptr};
          if (!low_q) {
            h_4 = use_helper_function(name, low_q, k, &jl_table);
            h_4->updateFrequency(omega);
          }

          for (const auto &a : orbs) {
            for (const auto &b : orbs) {

              // nb: These ops depend on _frequency_, not momentum..
              if (low_q && k == 0 && (name == "P" || name == "AT")) {
                const auto omega_ab = a.en() - b.en();
                h->updateFrequency(omega_ab);
                h_3->updateFrequency(omega_ab);
              }

              const auto rme = h->reducedME(a, b);
              const auto rad_int = h->radialIntegral(a, b);
              const auto C_ang = h->angularF(a.kappa(), b.kappa());
              const auto rme_2 = a * h->reduced_rhs(a.kappa(), b);

              const auto rme_3 = h_3->reducedME(a, b);
              const auto rme_4 = h_4 ? h_4->reducedME(a, b) : 0.0;

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
              if (h_4)
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

//==============================================================================
TEST_CASE("EM_multipole updateRank", "[DiracOperator][unit][EM_multipole]") {

  Wavefunction wf({75, 1.0e-4, 50.0, 1.0, "loglinear"}, {"Cs", 133, "Fermi"});

  auto &orbs = wf.valence();
  for (auto ik = 0ul; ik <= 5; ik++) {
    int kappa = Angular::kindex_to_kappa(ik);
    int n = Angular::l_k(kappa) + 1;
    orbs.push_back(DiracSpinor::exactHlike(n, kappa, wf.grid_sptr(), 1.0));
  }

  // type/component pairs covering all Lorentz structures
  struct OpSpec {
    char type;
    char comp;
  };
  const std::vector<OpSpec> ops = {
    {'V', 'E'}, {'V', 'M'}, {'V', 'L'}, {'V', 'T'}, {'A', 'E'},
    {'A', 'M'}, {'A', 'L'}, {'A', 'T'}, {'S', 'T'}, {'P', '_'},
  };

  const std::vector<int> Ks = {1, 2, 5};
  const std::vector<double> omegas = {1.0e-3, 1.0, 10.0};

  for (const auto &spec : ops) {

    // Single polymorphic operator constructed once at K=Ks[0], omega=omegas[0];
    // updated via updateRank + updateFrequency for each subsequent (k, omega).
    auto op_update = DiracOperator::MultipoleOperator(
      wf.grid(), 0, 0.0, spec.type, spec.comp, false);
    std::cout << spec.type << " " << spec.comp << std::endl;
    // std::cout << op_update->name() << "\n";
    REQUIRE(op_update != nullptr);

    for (int k : Ks) {
      for (double omega : omegas) {

        // Fresh operator constructed directly for this (k, omega)
        const auto op_fresh = DiracOperator::MultipoleOperator(
          wf.grid(), k, omega, spec.type, spec.comp, false);
        REQUIRE(op_fresh != nullptr);

        // Update the persistent operator to the same (k, omega)
        op_update->updateRank(k);
        op_update->updateFrequency(omega);

        for (const auto &a : orbs) {
          for (const auto &b : orbs) {

            if (!op_fresh->isZero(a, b)) {
              REQUIRE(!op_update->isZero(a, b));
            }
            if (op_fresh->isZero(a, b)) {
              REQUIRE(op_update->isZero(a, b));
              continue;
            }

            const auto rme_fresh = op_fresh->reducedME(a, b);
            const auto rme_update = op_update->reducedME(a, b);

            REQUIRE(rme_update == Approx(rme_fresh));
          }
        }
      }
    }
  }
}

//==============================================================================
TEST_CASE("EM_multipole clone", "[DiracOperator][unit][EM_multipole]") {

  Wavefunction wf({75, 1.0e-4, 50.0, 1.0, "loglinear"}, {"Cs", 133, "Fermi"});

  auto &orbs = wf.valence();
  for (auto ik = 0ul; ik <= 4; ik++) {
    int kappa = Angular::kindex_to_kappa(ik);
    int n = Angular::l_k(kappa) + 1;
    orbs.push_back(DiracSpinor::exactHlike(n, kappa, wf.grid_sptr(), 1.0));
  }

  const std::vector<double> omegas = {1.0e-3, 1.0, 10.0};
  const std::vector<int> Ks = {1, 2, 7}; // must be in order

  struct OpSpec {
    char type;
    char comp;
  };
  const std::vector<OpSpec> ops = {
    {'V', 'E'}, {'V', 'M'}, {'V', 'L'}, {'V', 'T'}, {'A', 'E'},
    {'A', 'M'}, {'A', 'L'}, {'A', 'T'}, {'S', '_'}, {'P', 'P'},
  };

  // Build a JL_table covering the omegas we will use
  using namespace qip::overloads;
  std::vector<double> qvec = omegas * PhysConst::alpha;
  const int max_L = Ks.back(); // must be in order
  SphericalBessel::JL_table jl_table(max_L + 1, qvec, wf.grid().r());

  for (const auto use_jl : {false, true}) {
    const SphericalBessel::JL_table *jl = use_jl ? &jl_table : nullptr;

    for (const auto &spec : ops) {
      for (int k : Ks) {
        const double omega0 = omegas[0];
        const double omega1 = omegas[1];

        // Construct original at (k, omega0)
        auto orig = DiracOperator::MultipoleOperator(
          wf.grid(), k, omega0, spec.type, spec.comp, false, jl);
        REQUIRE(orig != nullptr);

        // Clone at the same state
        auto cloned = orig->clone();
        REQUIRE(cloned != nullptr);

        // If jl table was provided, the clone should share the same pointer
        if (use_jl) {
          auto *orig_em =
            dynamic_cast<DiracOperator::EM_multipole *>(orig.get());
          auto *clone_em =
            dynamic_cast<DiracOperator::EM_multipole *>(cloned.get());
          REQUIRE(orig_em != nullptr);
          REQUIRE(clone_em != nullptr);
          REQUIRE(orig_em->jl() == clone_em->jl());
        }

        // Update original to omega1 — clone must still be at omega0
        orig->updateFrequency(omega1);

        for (const auto &a : orbs) {
          for (const auto &b : orbs) {
            if (orig->isZero(a, b))
              continue;

            // Clone at omega0 should be stable (unchanged by update to orig)
            const auto rme_clone = cloned->reducedME(a, b);
            const auto rme_orig = orig->reducedME(a, b);
            // same K, so if one is zero, both zero
            if (rme_orig == 0.0) {
              REQUIRE(rme_clone == rme_orig);
            } else {
              REQUIRE(rme_clone != rme_orig);
            }
            REQUIRE(rme_clone == Approx(cloned->reducedME(a, b)));
          }
        }

        // Now update clone to omega1 as well — results must match original
        cloned->updateFrequency(omega1);

        for (const auto &a : orbs) {
          for (const auto &b : orbs) {
            if (orig->isZero(a, b)) {
              REQUIRE(cloned->isZero(a, b));
              continue;
            }
            REQUIRE(cloned->reducedME(a, b) == Approx(orig->reducedME(a, b)));
          }
        }
      }
    }
  }
}