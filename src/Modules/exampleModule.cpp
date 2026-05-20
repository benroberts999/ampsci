#include "DiracOperator/GenerateOperator.hpp"
#include "DiracOperator/Operators/RadialF.hpp"
#include "DiracOperator/include.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/Modules.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include <iostream>

// exampleModule.cpp
//
// Template for writing a custom ampsci module.
//
// HOW TO ADD A NEW MODULE:
//   1. Copy this file. Rename it (e.g., myModule.cpp).
//   2. Rename `exampleModule` to your module name throughout this file
//      (function, registrar string, and function pointer).
//   3a. Internal: drop the file anywhere under src/ and recompile.
//   3b. External: place the file anywhere outside of ampsci, add to your Makefile:
//         EXTERNAL_MODULES = path/to/myModule.cpp
//       then recompile.
//   No other files need to be edited.

namespace Module {

// Declare out module function
// (This is optional; but decalring first allows us to register at the
//  top of the file. Purely preference.)
void exampleModule(const IO::InputBlock &input, const Wavefunction &wf);

// Register module:
// The first string is the name used in input files: Module::exampleModule{}.
// Note: must be inside an anonymous namespace. Just copy this example
namespace {
const Register r_exampleModule{
  "exampleModule",
  "Example module: Simple example of things you can do do in a module",
  &exampleModule};
} // namespace

//==============================================================================
//==============================================================================

// Now, define the function body of the module:

void exampleModule(const IO::InputBlock &input, const Wavefunction &wf) {
  // This is an example module, designed to help you write new modules.
  //
  // It demonstrates:
  //  - Parsing user input options
  //  - Constructing a DiracOperator directly, and via generate()
  //  - Looping over core and valence states
  //  - Constructing a second Wavefunction and comparing to the external one

  // Declare and validate input options.
  // Users see these descriptions when they run: ./ampsci -m exampleModule
  input.check(
    {{"operator",
      "Operator for matrix elements. "
      "Run 'ampsci -o' for a list of available operators. [default: E1]"},
     {"options{}",
      "Optional options for the operator (see ampsci -o 'operator')"},
     {"var_alpha", "Fractional variation in alpha for comparison wf: "
                   "alpha = var_alpha * alpha_0 [0.01]"}});
  if (input.has_option("help"))
    return;

  const auto oper = input.get<std::string>("operator", "E1");
  // Use user-supplied options block if present, otherwise empty (defaults)
  const auto h_options =
    input.getBlock("options").value_or(IO::InputBlock(oper, {}));

  const auto var_alpha = input.get("var_alpha", 0.01);

  //----------------------------------------------------------------------------
  // 1. Construct a scalar |r| operator directly (not via generate()).
  //    RadialF with power=1 gives the scalar operator with radial factor r,
  //    so <a|r|b> = int[ (fa*fb + ga*gb) r dr ] \delta_{kappa_a,kappa_b}
  const DiracOperator::RadialF hr(wf.grid(), 1.0);  // r
  const DiracOperator::RadialF hr2(wf.grid(), 2.0); // r^2

  std::cout << "\nCore orbitals: energy and <r> expectation value\n";
  std::cout << "\n  state    energy (au)       <r> (aB)     r_rms (aB)\n";

  for (const auto &Fc : wf.core()) {
    // For a scalar operator, radialIntegral(a,a) = <a|r|a>
    const double r_exp = hr.radialIntegral(Fc, Fc);
    const double r2_exp = hr2.radialIntegral(Fc, Fc);
    const auto r_rms = std::sqrt(r2_exp);
    fmt::print("  {:5s}  {:13.6f}  {:13.6e}  {:13.6e}\n", Fc.shortSymbol(),
               Fc.en(), r_exp, r_rms);
  }
  std::cout << "\n";

  //----------------------------------------------------------------------------
  // 2. Construct the operator via generate() -- this is the runtime approach,
  //    which reads operator name and options from user input.
  const auto h = DiracOperator::generate(oper, h_options, wf);

  std::cout << "Valence reduced matrix elements: " << h->name() << "\n";
  std::cout << "  " << std::string(36, '-') << "\n";
  fmt::print("  {:5s}  {:5s}  {:>13s}\n", "a", "b", "<a||h||b>");
  std::cout << "  " << std::string(36, '-') << "\n";
  for (const auto &Fa : wf.valence()) {
    for (const auto &Fb : wf.valence()) {
      if (Fb.en() >= Fa.en())
        continue; // upper triangle only
      if (h->isZero(Fa, Fb))
        continue;
      fmt::print("  {:5s}  {:5s}  {:+13.6e}\n", Fa.shortSymbol(),
                 Fb.shortSymbol(), h->reducedME(Fa, Fb));
    }
  }
  std::cout << "\n";

  //----------------------------------------------------------------------------
  // 3. Construct a new Wavefunction with a different alpha (or method),
  //    solve HF, then compare matrix elements state-by-state.
  //    getState(n, kappa) looks up a state by quantum numbers -- returns
  //    nullptr if not found. This is safe even if the two wfs differ.

  fmt::print("Constructing comparison wavefunction:\n");
  fmt::print("  var_alpha = {}\n\n", var_alpha);

  // Same grid and nucleus, different alpha
  Wavefunction wf2(wf.grid().params(), wf.nucleus(), var_alpha);
  wf2.solve_core(wf.vHF()->method(), wf.coreConfiguration());
  wf2.solve_valence(DiracSpinor::state_config(wf.valence()));

  const auto h2 = DiracOperator::generate(oper, h_options, wf2);

  fmt::print("Comparison: {} matrix elements\n", h->name());
  fmt::print("  (alpha_ref = {:.6f}, alpha_comp = {:.6f})\n",
             wf.alpha() / PhysConst::alpha, wf2.alpha() / PhysConst::alpha);
  std::cout << "  " << std::string(58, '-') << "\n";
  fmt::print("  {:5s}  {:5s}  {:>13s}  {:>13s}  {:>13s}\n", "a", "b",
             "ME (ref)", "ME (comp)", "delta");
  std::cout << "  " << std::string(58, '-') << "\n";

  for (const auto &Fa : wf.valence()) {
    for (const auto &Fb : wf.valence()) {
      if (Fb.en() >= Fa.en())
        continue;
      if (h->isZero(Fa, Fb))
        continue;

      // Look up matching states by quantum numbers; skip if missing in wf2
      const auto *Fa2 = wf2.getState(Fa.n(), Fa.kappa());
      const auto *Fb2 = wf2.getState(Fb.n(), Fb.kappa());
      if (!Fa2 || !Fb2)
        continue;

      const double me = h->reducedME(Fa, Fb);
      const double me2 = h2->reducedME(*Fa2, *Fb2);
      fmt::print("  {:5s}  {:5s}  {:+13.6e}  {:+13.6e}  {:+13.6e}\n",
                 Fa.shortSymbol(), Fb.shortSymbol(), me, me2, me2 - me);
    }
  }
}

} // namespace Module
