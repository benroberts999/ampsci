#include "DiracOperator/GenerateOperator.hpp"
#include "Maths/Grid.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include "include.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
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

  // Worst violation, per operator, of the hermiticity relation between
  // <a||h||b> and <b||h||a> (see below)
  std::map<std::string, double> worst_symmetry;

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

              // 5. Hermiticity: <b||h||a> = symm_sign <a||h||b>, where the
              // sign depends on whether the stored reduced matrix elements
              // are the real or the imaginary part (Realness). Only the
              // exact operators are checked: the low-q forms at K = 0 and
              // the length form VEk_Len hold on shell (omega = e_a - e_b)
              // only, so at fixed omega they mix both swap symmetries (as
              // E1v); VEk_Len is checked against E1 below
              if (!low_q && name != "VE_Len") {
                const auto rme_ba = h->reducedME(b, a);
                const auto expected = h->symm_sign(a, b) * rme;
                const auto scale = std::max(std::abs(rme), 1.0e-30);
                auto &worst = worst_symmetry[name];
                worst = std::max(worst, std::abs(rme_ba - expected) / scale);
              }
            }
          }
        }
      }
    }
  }

  // Hermiticity (5): every operator must obey its own symmetry rule
  std::cout << "Hermiticity of the reduced matrix elements:\n";
  for (const auto &[name, worst] : worst_symmetry) {
    fmt::print(
      "  {:>6}: worst |<b||h||a> - s <a||h||b>| / |<a||h||b>| = {:.1e}\n", name,
      worst);
  }
  for (const auto &[name, worst] : worst_symmetry) {
    INFO(name);
    REQUIRE(worst < 1.0e-8);
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

        // With the Hermitian (uniform i^{K+1}) convention T^(+1)_1 ->
        // +(sqrt2/3)<alpha>, matching E1 = -|e| r and E1v directly
        REQUIRE(E1.reducedME(a, b) ==
                Approx(ff * E1_w.reducedME(a, b)).epsilon(eps));

        REQUIRE(E1v.reducedME(a, b) ==
                Approx(ff * E1v_w.reducedME(a, b)).epsilon(eps));

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

        // Update original to omega1 -- clone must still be at omega0
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

        // Now update clone to omega1 as well -- results must match original
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

//==============================================================================
TEST_CASE("EM_multipole low-q vs full-q",
          "[DiracOperator][unit][EM_multipole]") {
  // The low-q (long-wavelength) forms must reproduce the full operators as
  // q -> 0, including the sign: this pins the relative phase conventions of
  // the two implementations (the 2026-09 audit found VEk_lowq, VLk_lowq and
  // S5k with the wrong sign, undetected by any test).
  // Exact H-like (Z=1) Dirac orbitals; a second orbital of the same kappa is
  // included for the first three kappas so that the orthogonal same-kappa
  // pairs used by the kappa_a = kappa_b branches (AE, AL, K=0 S) are
  // exercised. Diagonal elements are skipped where the low-q form assumes
  // <a|b> = 0 (AE, AL, S). Not tested: K=0 vector temporal (an O(q^2) term,
  // needs exact orthogonality) and K=0 vector longitudinal (no low-q form
  // is implemented; the paper's B36 has one).
  using namespace DiracOperator;

  Wavefunction wf({2000, 1.0e-5, 60.0, 1.0, "loglinear"}, {"Cs", 133, "Fermi"});
  auto &orbs = wf.valence();
  for (auto ik = 0ul; ik <= 5; ik++) {
    const int kappa = Angular::kindex_to_kappa(ik);
    const int l = Angular::l_k(kappa);
    orbs.push_back(DiracSpinor::exactHlike(l + 1, kappa, wf.grid_sptr(), 1.0));
    if (ik <= 2)
      orbs.push_back(
        DiracSpinor::exactHlike(l + 2, kappa, wf.grid_sptr(), 1.0));
  }

  // q = alpha*omega; (qr)^2 < 1e-9 over the grid, so the low-q forms hold
  // to that order (the K=0 AT and P forms use the transition frequency
  // internally, via the exact gamma^5 relations, exact for H-like states)
  const double omega = 1.0e-5;
  const auto &gr = wf.grid();

  const auto make = [&](const std::string &name, int k,
                        bool low_q) -> std::unique_ptr<TensorOperator> {
    if (low_q) {
      if (name == "VE")
        return std::make_unique<VEk_lowq>(gr, k, omega);
      if (name == "VL")
        return std::make_unique<VLk_lowq>(gr, k, omega);
      if (name == "VM")
        return std::make_unique<VMk_lowq>(gr, k, omega);
      if (name == "VT")
        return std::make_unique<Phik_lowq>(gr, k, omega);
      if (name == "AE")
        return std::make_unique<AEk_lowq>(gr, k, omega);
      if (name == "AL")
        return std::make_unique<ALk_lowq>(gr, k, omega);
      if (name == "AM")
        return std::make_unique<AMk_lowq>(gr, k, omega);
      if (name == "AT")
        return std::make_unique<Phi5k_lowq>(gr, k, omega);
      if (name == "S")
        return std::make_unique<Sk_lowq>(gr, k, omega);
      if (name == "P")
        return std::make_unique<S5k_lowq>(gr, k, omega);
    } else {
      if (name == "VE")
        return std::make_unique<VEk>(gr, k, omega);
      if (name == "VL")
        return std::make_unique<VLk>(gr, k, omega);
      if (name == "VM")
        return std::make_unique<VMk>(gr, k, omega);
      if (name == "VT")
        return std::make_unique<Phik>(gr, k, omega);
      if (name == "AE")
        return std::make_unique<AEk>(gr, k, omega);
      if (name == "AL")
        return std::make_unique<ALk>(gr, k, omega);
      if (name == "AM")
        return std::make_unique<AMk>(gr, k, omega);
      if (name == "AT")
        return std::make_unique<Phi5k>(gr, k, omega);
      if (name == "S")
        return std::make_unique<Sk>(gr, k, omega);
      if (name == "P")
        return std::make_unique<S5k>(gr, k, omega);
    }
    return nullptr;
  };

  // The K=0 pseudoscalar elements are alpha^3 w^2 suppressed and come from
  // a near-cancellation in the full operator (P^(+)[j_0]), so the full-q
  // values are only good to ~1% on this grid for the smallest w
  struct LowqCase {
    std::string name;
    int k;
    bool skip_diagonal;
    double tol;
  };
  const std::vector<LowqCase> cases{
    {"VE", 1, false, 1.0e-4}, {"VL", 1, false, 1.0e-4},
    {"VM", 1, false, 1.0e-4}, {"VT", 1, false, 1.0e-4},
    {"AE", 1, true, 1.0e-4},  {"AL", 1, true, 1.0e-4},
    {"AM", 1, false, 1.0e-4}, {"AT", 1, false, 1.0e-4},
    {"S", 1, false, 1.0e-4},  {"P", 1, false, 1.0e-4},
    {"AL", 0, false, 1.0e-4}, {"AT", 0, false, 1.0e-4},
    {"P", 0, false, 3.0e-2},  {"S", 0, true, 1.0e-4}};

  std::cout
    << "Low-q vs full-q multipoles (worst |low-q - full| / max|full|):\n";
  for (const auto &[name, k, skip_diagonal, tol] : cases) {
    const auto h_full = make(name, k, false);
    const auto h_lowq = make(name, k, true);
    REQUIRE(h_full != nullptr);
    REQUIRE(h_lowq != nullptr);

    std::vector<std::tuple<std::string, double, double>> pairs;
    double max_abs = 0.0;
    for (const auto &a : orbs) {
      for (const auto &b : orbs) {
        const bool diagonal = (a.n() == b.n() && a.kappa() == b.kappa());
        if (h_full->isZero(a, b) || (skip_diagonal && diagonal))
          continue;
        const auto full = h_full->reducedME(a, b);
        const auto lowq = h_lowq->reducedME(a, b);
        pairs.emplace_back(a.shortSymbol() + " " + b.shortSymbol(), full, lowq);
        max_abs = std::max(max_abs, std::abs(full));
      }
    }
    REQUIRE(!pairs.empty());
    REQUIRE(max_abs > 0.0);

    double worst = 0.0;
    for (const auto &[label, full, lowq] : pairs) {
      INFO(name << " K=" << k << " <" << label << ">: full = " << full
                << ", low-q = " << lowq);
      REQUIRE(lowq == Approx(full).epsilon(tol).margin(tol * max_abs));
      worst = std::max(worst, std::abs(lowq - full) / max_abs);
    }
    fmt::print("  {:>2} K={}: {:.1e}\n", name, k, worst);
  }
}
