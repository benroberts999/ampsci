#include "Modules/exampleModule.hpp"
#include "DiracOperator/DiracOperator.hpp" //For E1 operator
#include "HF/Breit.hpp"
#include "IO/InputBlock.hpp"
#include "LinAlg/LinAlg.hpp"
#include "MBPT/Feynman.hpp"
#include "MBPT/RDMatrix.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <vector>

namespace Module {

// want this function to construct the VBr coordinate operator
std::vector<MBPT::GMatrix> form_VBr(const IO::InputBlock &input,
                                    const Wavefunction &wf) {
  // returns pointer to HF class that will have core and valence orbitals for the atom I give to the input file
  auto vHF = wf.vHF();

  const auto r0 = input.get({"Grid"}, "r0", 1.0e-6);
  const auto rmax = input.get({"Grid"}, "rmax", 120.0);

  const auto default_stride = [&]() {
    // By default, choose stride such that there is 150 points over [1e-4,30]
    const auto dstride =
        int(wf.grid().getIndex(30.0) - wf.grid().getIndex(1.0e-4)) / 150;
    return (dstride <= 2) ? 2 : dstride;
  }();
  const auto stride = input.get({"Correlations"}, "stride", default_stride);
  const auto i0 = vHF->grid().getIndex(r0);
  const auto size = (vHF->grid().getIndex(rmax) - i0) / stride + 1;

  auto max_ki = 2 * DiracSpinor::max_l(wf.core());

  // ---- CLASS ACTUALLY STARTS DOING STUFF HERE ---- //
  auto Vx_kappa = std::vector<MBPT::GMatrix>(
      std::size_t(max_ki + 1), {i0, stride, size, false, vHF->grid_sptr()});

  // initialises the hydrogen wave function object
  Wavefunction wfH(vHF->grid_sptr(), {"1", 0, "Ball"});

  wfH.set_HF();
  wfH.solve_core(false);
  wfH.formBasis(SplineBasis::Parameters(
      "90spdfghi", 90, 9, 1e-4, 0.0, 90.0,
      true)); // forms the hydrogen wavefunction basis; NOTE: the last boolean option is for including the negative energy states in the basis

  // constructs the Breit matrix as VBr(r1,r2)=\sum_i[VBrF_i](r1)F_i^†(r2)
  for (const auto &Fn : wfH.basis()) {
    auto VxFn = vHF->VBr(Fn);
    if (Fn.k_index() < max_ki + 1) {
      Vx_kappa[Fn.k_index()].add(
          VxFn, Fn,
          1.0); // multiplies VBrFa(r1) by Fa^†(r2) from the right
    }
  }

  // includes integration measures
  for (auto &Vx : Vx_kappa) {
    Vx.dri_in_place();
    Vx.drj_in_place();
  }
  return Vx_kappa;
}

/*
// want this function to calculate matrix elements of the VBr operator with respect to core and valence HF states for the atom i give into the input file
void exampleModule(const IO::InputBlock &input, const Wavefunction &wf) {

  // returns pointer to HF class that will have core and valence orbitals for the atom I give to the input file
  auto vHF = wf.vHF();

  const auto r0 = input.get({"Grid"}, "r0", 1.0e-6);
  const auto rmax = input.get({"Grid"}, "rmax", 120.0);

  const auto default_stride = [&]() {
    // By default, choose stride such that there is 150 points over [1e-4,30]
    const auto dstride =
        int(wf.grid().getIndex(30.0) - wf.grid().getIndex(1.0e-4)) / 150;
    return (dstride <= 2) ? 2 : dstride;
  }();
  const auto stride = input.get({"Correlations"}, "stride", default_stride);
  const auto i0 = vHF->grid().getIndex(r0);
  const auto size = (vHF->grid().getIndex(rmax) - i0) / stride + 1;

  auto max_ki = 2 * DiracSpinor::max_l(wf.core());

  // construct coordinate VBr matrix
  std::vector<MBPT::GMatrix> VBr = form_VBr(input, wf);

  std::cout << "State         Regular      Coordinate matrix" << std::endl;

  // loop over all core orbitals of the HF core for the atom I give to the input file, not the hydrogen core
  for (const auto &Fi : wf.valence()) {
    auto expval_gold = Fi * vHF->VBr(Fi);
    auto expval_feyn =
        Fi * (VBr.at(std::size_t(Angular::indexFromKappa(Fi.kappa()))) * Fi);

    std::cout << Fi.symbol() << "      " << expval_gold << "       "
              << expval_feyn << std::endl;
  }
}
*/

// this is my attempt at writing the code for forming the BCHF basis using Derevianko's method of constructing the B-spline basis
void exampleModule(const IO::InputBlock &input, const Wavefunction &wf) {

  // extracting only the s orbital basis states
  auto s_basis =
      qip::select_if(wf.basis(), [](const auto &f) { return f.kappa() == -1; });

  int N = s_basis.size(); // s-orbital basis size

  // allows us to calculate [VBrFa](r) without BCHF basis from input file
  HF::Breit Br(1.0);

  std::cout << "Basis size is " << N << std::endl;
  auto vHF = wf.vHF();  // pointer to HF
  auto basis = s_basis; // vector of spline basis states

  // extracts energy of highest core electron, Fermi energy
  //auto Fermi_en = vHF->core()[vHF->core().size() - 1].en();

  //std::cout << "Fermi energy is " << Fermi_en << std::endl;

  // initialise matrix to store (dU)ik + ei where dU = U_BCHF - U_CHF = VBr
  LinAlg::Matrix<double> dUik_plus_ei{N, N};

#pragma omp parallel for collapse(2)
  for (int i = 0; i < N; i++) {
    for (int k = 0; k < N; k++) {

      // extracting the ith and kth basis state
      const auto &Fi = basis[i];
      const auto &Fk = basis[k];

      auto dU_ik = Fi * Br.VbrFa(Fk, wf.core()); // (VBr)_ik

      dUik_plus_ei[i][k] = dU_ik;

      if (i == k) {
        dUik_plus_ei[i][k] += Fi.en();
      }
    }
  }
  // solves eigenvalue equation for BCHF eigenvalue varepsilonj
  auto [e_values, e_vectors] = LinAlg::symmhEigensystem(dUik_plus_ei);

  std::cout << std::endl;

  for (int j = 0; j < s_basis.size(); j++) {
    // only print energies of states with energies above Fermi level?
    //if (e_values.size() < Fermi_en) {
    std::cout << e_values[j] * PhysConst::Hartree_invcm << " "
              << s_basis.at(j).en() * PhysConst::Hartree_invcm << std::endl;
    //}
  }
}

// =================================================================================================================== //
/*
void exampleModule(const IO::InputBlock &input, const Wavefunction &wf) {
  // This is an example module, designed to help you write new modules

  // In this example, we will solve a new wavefunction, assuming a different
  // nuclear charge distribution, and see the difference in the energies and E1
  // matrix elements this produces.

  // Read in some optional input options: A, rms, and type
  // "checkBlock" is optional, but helps prevent errors on input options:
  input.check({{"A", "Atomic mass number [0]"},
               {"rrms", "root-mean-square nuclear radii [-1]"},
               {"type", "Fermi/ball/pointlike [Fermi]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  // If no A is specified, use 0 (i.e., poinlike)
  auto A = input.get("A", 0);
  // if rrms<0, the default value (for the given A) will be looked up
  auto rrms = input.get("rrms", -1.0);
  // Will assume a Fermi nucleus, unless A or r =0
  auto nuc_type = input.get<std::string>("type", "Fermi");
  // Set nuc. type explicitely to 'pointlike' if A=0, or r_rms = 0.0
  if (A == 0 || rrms == 0.0)
    nuc_type = "pointlike";
  if (nuc_type == "pointlike") {
    A = 0;
    rrms = 0.0;
  }

  // Use the same Grid and alpha, but different nuclear parameters (except Z)
  Wavefunction wfpt(wf.grid().params(),
                    {wf.Znuc(), A, nuc_type, rrms, Nuclear::default_t},
                    wf.alpha() / PhysConst::alpha);

  std::cout << "\n";
  IO::print_line();
  std::cout << "Calculating finite-nuclear size corrections for \n"
            << wf.atom() << ", " << wf.nucleus() << "\nvs.\n"
            << wfpt.atom() << ", " << wfpt.nucleus() << "\n\n";

  // Note: Assume only Hartree-Fock approximation, no Breit, QED, Sigma etc.

  // Solve Hartree-Fock core for new wavefuinction:
  wfpt.solve_core("HartreeFock", 0.0, wf.coreConfiguration());
  // print the new core energies:
  wfpt.printCore();

  // Solve Hartree-Fock valence
  wfpt.solve_valence(DiracSpinor::state_config(wf.valence()));
  wfpt.printValence();

  // Calculate the energy shifts in atomic units, and GHz
  std::cout << "\nFinite nuclear charge energy shifts:\n";
  std::cout << "  state            au           GHz\n";
  for (auto i = 0ul; i < wf.valence().size(); ++i) {
    const auto del_e = wf.valence()[i].en() - wfpt.valence()[i].en();
    printf("%7s  %12.5e  %12.5e\n", wf.valence()[i].symbol().c_str(), del_e,
           del_e * PhysConst::Hartree_GHz);
  }

  // Now, we calculate E1 matrix elements.
  std::cout << "\nFinite nuclear charge shifts to E1 reduced MEs (no RPA):\n";

  std::cout << "                      A=" << wf.Anuc()
            << "         A=" << wfpt.Anuc() << "         Shift\n";

  // 1) Create E1 operator:
  const auto hE1 = DiracOperator::E1(wf.grid());

  // 2) Loop through each pair of valence states, calc E1 matrix elements:
  for (auto a = 0ul; a < wf.valence().size(); ++a) {
    const auto &Fa = wf.valence()[a];
    const auto &F0a = wfpt.valence()[a]; // pointlike orbital
    for (auto b = 0ul; b < a; ++b) {
      const auto &Fb = wf.valence()[b];
      const auto &F0b = wfpt.valence()[b]; // pointlike orbital

      // Skip the MEs which are zero due to selection rules:
      if (hE1.isZero(Fa.kappa(), Fb.kappa()))
        continue;

      const auto e1 = hE1.reducedME(Fa, Fb);
      const auto e10 = hE1.reducedME(F0a, F0b);
      const auto delta = e1 - e10;

      std::cout << Fa << " - " << Fb << ": ";
      printf("%12.5e  %12.5e  %12.5e\n", e1, e10, delta);
    }
  }
}
*/

} // namespace Module
