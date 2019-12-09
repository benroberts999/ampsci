#include "Dirac/Wavefunction.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/UserInput.hpp"
#include "Maths/BSplines.hpp"
#include "Modules/Module_runModules.hpp"
#include <iostream>
#include <string>
//
#include "Dirac/Operators.hpp"
#include "Maths/Matrix_linalg.hpp"
#include "testSplines.hpp"

int main(int argc, char *argv[]) {
  ChronoTimer timer("\nhartreeFock");
  std::string input_file = (argc > 1) ? argv[1] : "hartreeFock.in";
  std::cout << "Reading input from: " << input_file << "\n";

  // Input options
  UserInput input(input_file);

  // Get + setup atom parameters
  auto input_ok = input.check("Atom", {"Z", "A", "varAlpha2"});
  auto atom_Z = input.get<std::string>("Atom", "Z");
  auto atom_A = input.get("Atom", "A", -1);
  auto var_alpha = [&]() {
    auto varAlpha2 = input.get("Atom", "varAlpha2", 1.0);
    return (varAlpha2 > 0) ? std::sqrt(varAlpha2) : 1.0e-25;
  }();

  // Get + setup Grid parameters
  input_ok = input_ok && input.check("Grid", {"r0", "rmax", "num_points",
                                              "type", "b", "fixed_du"});
  auto r0 = input.get("Grid", "r0", 1.0e-5);
  auto rmax = input.get("Grid", "rmax", 150.0);
  auto num_points = input.get("Grid", "num_points", 1600ul);
  auto du_tmp = input.get("Grid", "fixed_du", -1.0); // >0 means calc num_points
  if (du_tmp > 0)
    num_points = 0;
  auto b = input.get("Grid", "b", 0.33 * rmax);
  auto grid_type = input.get<std::string>("Grid", "type", "loglinear");
  if (b <= r0 || b >= rmax)
    grid_type = "logarithmic";

  // Get + setup nuclear parameters
  input_ok =
      input_ok && input.check("Nucleus", {"A", "rrms", "skin_t", "type"});
  atom_A = input.get("Nucleus", "A", atom_A); // over-writes "atom" A
  auto nuc_type = input.get<std::string>("Nucleus", "type", "Fermi");
  auto rrms = input.get("Nucleus", "rrms", -1.0); /*<0 means lookup default*/
  auto skint = input.get("Nucleus", "skin_t", -1.0);

  // create wavefunction object
  Wavefunction wf(atom_Z, {num_points, r0, rmax, b, grid_type, du_tmp},
                  {atom_Z, atom_A, nuc_type, rrms, skint}, var_alpha);

  // Parse input for HF method
  input_ok =
      input_ok &&
      input.check("HartreeFock", {"core", "valence", "convergence", "method",
                                  "Green_H", "Green_d", "Tietz_g", "Tietz_t",
                                  "orthonormaliseValence", "sortOutput"});
  auto str_core = input.get<std::string>("HartreeFock", "core", "[]");
  auto eps_HF = input.get("HartreeFock", "convergence", 1.0e-12);
  auto HF_method = HartreeFock::parseMethod(
      input.get<std::string>("HartreeFock", "method", "HartreeFock"));

  if (!input_ok)
    return 1;

  std::cout << "\nRunning for " << wf.atom() << "\n"
            << wf.nuclearParams() << "\n"
            << wf.rgrid.gridParameters() << "\n"
            << "********************************************************\n";
  // For when using Hartree, or a parametric potential:
  double H_d = 0.0, g_t = 0.0;
  if (HF_method == HFMethod::GreenPRM) {
    H_d = input.get("HartreeFock", "Green_H", 0.0);
    g_t = input.get("HartreeFock", "Green_d", 0.0);
    std::cout << "Using Greens Parametric Potential\n";
  } else if (HF_method == HFMethod::TietzPRM) {
    H_d = input.get("HartreeFock", "Tietz_g", 0.0);
    g_t = input.get("HartreeFock", "Tietz_t", 0.0);
    std::cout << "Using Tietz Parametric Potential\n";
  } else if (HF_method == HFMethod::Hartree) {
    std::cout << "Using Hartree Method (no Exchange)\n";
  }

  { // Solve Hartree equations for the core:
    ChronoTimer t(" core");
    wf.hartreeFockCore(HF_method, str_core, eps_HF, H_d, g_t);
  }

  // Adds effective polarision potential to direct potential
  // (After HF core, before HF valence)
  auto a_eff = input.get("dV", "a_eff", 0.0);
  if (a_eff > 0) { // a=0.61 works well for Cs ns, n=6-18
    auto r_cut = input.get("dV", "r_cut", 1.0);
    auto dV = [=](double x) { return -0.5 * a_eff / (x * x * x * x + r_cut); };
    for (auto i = 0u; i < wf.rgrid.num_points; ++i) {
      wf.vdir[i] += dV(wf.rgrid.r[i]);
    }
  }

  // Solve for the valence states:
  auto valence_list = (wf.Ncore() < wf.Znuc())
                          ? input.get<std::string>("HartreeFock", "valence", "")
                          : "";
  if (valence_list != "") {
    // 'if' is only for output format, nothing bad happens if below are called
    ChronoTimer t("  val");
    wf.hartreeFockValence(valence_list);
    if (input.get("HartreeFock", "orthonormaliseValence", false))
      wf.orthonormaliseOrbitals(wf.valence_orbitals, 2);
  }

  // Output results:
  std::cout << "\nHartree Fock: " << wf.atom() << "\n";
  auto sorted = input.get("HartreeFock", "sortOutput", true);
  wf.printCore(sorted);
  wf.printValence(sorted);

  // run each of the modules
  Module::runModules(input, wf);
  // return 1;

  // Get + setup nuclear parameters
  // auto basis_ok =
  input.check("Basis", {"number", "order", "r0", "rmax", "print"});
  auto n_spl = input.get("Basis", "number", 0ul);
  auto k_spl = input.get("Basis", "order", 0ul);
  auto r0_spl = input.get("Basis", "r0", 0.0);
  auto rmax_spl = input.get("Basis", "rmax", 0.0);
  // auto print_spl = input.get("Basis", "print", false);

  // if (k_spl > 0 && n_spl >= k_spl && basis_ok) {
  // BSplines bspl(n_spl, k_spl, wf.rgrid, r0_spl, rmax_spl);

  auto basis = test_splines(-1, n_spl, k_spl, r0_spl, rmax_spl, wf.rgrid);
  // std::cin.get();
  // auto basis = wf.core_orbitals;

  auto Hd = DirectHamiltonian(wf.vnuc, wf.vdir, wf.get_alpha());

  std::cout << "\n\n";
  LinAlg::SqMatrix Aij((int)basis.size());
  LinAlg::SqMatrix Sij((int)basis.size());
#pragma omp parallel for
  for (auto i = 0; i < (int)basis.size(); i++) {
    const auto &si = basis[i];
    for (auto j = 0; j < (int)basis.size(); j++) {
      const auto &sj = basis[j];
      auto VexPsi_j = HartreeFock::vex_psia_any(sj, wf.core_orbitals);
      auto VexPsi_i = HartreeFock::vex_psia_any(si, wf.core_orbitals);

      // Vex seems to make incredibly small contribution??
      // auto aij = Hd.matrixEl(si, sj) + 0.5 * (si * VexPsi_j + (sj *
      // VexPsi_i));
      auto aij = Hd.matrixEl(si, sj) + (si * VexPsi_j);
      // auto aji = Hd.matrixEl(sj, si) + (sj * VexPsi_i);
      Aij[i][j] = aij; // / std::sqrt((si * si) * (sj * sj));
      Sij[i][j] = si * sj;
    }
    // std::cout << "\n";
  }
  std::cout << "\nFilled Matrix\n\n";

  // Aij.clip_low(1.0e-9); //?
  // Sij.clip_low(1.0e-8); //?
  // Aij.make_symmetric(); //?

  std::cout << "Worst A:" << Aij.check_symmetric() << "\n";
  std::cout << "Worst S:" << Sij.check_symmetric() << "\n";
  // std::cin.get();

  LinAlg::test3(Aij, Sij);

  //*********************************************************
  //               TESTS
  //*********************************************************

  // // needs: #include "Physics/AtomData.hpp" (for AtomData::listOfStates_nk)
  // bool test_hf_basis = false;
  // if (test_hf_basis) {
  //   auto basis_lst = AtomData::listOfStates_nk("9spd8f");
  //   std::vector<DiracSpinor> basis = wf.core_orbitals;
  //   HartreeFock hfbasis(wf, basis, 1.0e-6);
  //   hfbasis.verbose = false;
  //   for (const auto &nk : basis_lst) {
  //     if (wf.isInCore(nk.n, nk.k))
  //       continue;
  //     basis.emplace_back(DiracSpinor(nk.n, nk.k, wf.rgrid));
  //     auto tmp_vex = std::vector<double>{};
  //     hfbasis.hf_valence_approx(basis.back(), tmp_vex);
  //   }
  //   wf.orthonormaliseOrbitals(basis, 2);
  //   wf.printValence(false, basis);
  //   std::cout << "\n Total time: " << timer.reading_str() << "\n";
  // }

  return 0;
}

//******************************************************************************
