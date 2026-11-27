#include "Modules/StructRad.hpp"
#include "DiracOperator/GenerateOperator.hpp"
#include "ExternalField/CorePolarisation.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/FeynmanSR.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Modules/Modules.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace Module {

// Declare, register, then define below.
namespace {
const Register r_StructRad{"StructRad", "description", &StructRad};
} // namespace

void StructRad(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check({
    {"operator", "e.g., E1, hfs [E1]"},
    {"options{}", "options specific to operator; blank by dflt"},
    {"rpa", "RPA for the internal operator line: true(=TDHF), false, TDHF, "
            "basis, diagram [false]"},
    {"rpa_omega", "Frequency the RPA is solved at [0.0]"},
    {"stride", "stride [such that ~150 points between r0 and rmax]"},
    {"n_min_core", "Minimum core state [1]"},
    {"rmax", "Maximum r value [40]"},
    {"r0", "Minimum r value, [1e-4]"},
    {"includeG", "include lower g component [false]"},
    {"exchange", "include the exchange diagrams [false]"},
    {"screening", "Include all-orders screening [false]"},
    {"hole_particle", "Include all-orders hole-particle interaction [false]"},
    {"both_orders", "Also calculate each pair with w and v swapped "
                    "(numerical test) [false]"},
  });

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  IO::ChronoTimer timer("StructRad");

  const auto oper = input.get<std::string>("operator", "E1");
  const auto h_options =
    input.getBlock("options").value_or(IO::InputBlock(oper, {}));
  const auto h = DiracOperator::generate(oper, h_options, wf);

  // Core polarisation of the external field, solved once at a fixed frequency
  const auto rpa_method = input.get<std::string>("rpa", "false");
  const auto rpa_omega = input.get("rpa_omega", 0.0);
  const auto rpa = ExternalField::make_rpa(rpa_method, h.get(), wf.vHF(), true,
                                           wf.basis(), wf.identity());
  if (rpa) {
    rpa->solve_core(rpa_omega);
  }

  auto r0 = input.get("r0", 1e-4);
  auto rmax = input.get("rmax", 40.0);
  const auto i0 = wf.grid().getIndex(r0);
  const auto imax = wf.grid().getIndex(rmax);
  // By default, choose stride such that there are ~150 points between r0
  // and rmax
  const auto default_stride = std::max((imax - i0) / 150, 1ul);
  const std::size_t stride = input.get("stride", default_stride);
  const auto n_min_core = input.get("n_min_core", 1);
  const bool includeG = input.get("includeG", false);
  // Sub-grid as for the Correlations Feynman (so Q*Pi*Q cache can be shared)
  const auto gsize = (imax - i0) / stride + 1;
  const auto screening = input.get("screening", false);
  const auto hole_particle = input.get("hole_particle", false);
  const auto both_orders = input.get("both_orders", false);

  // Real part of the frequency contour, chosen as for the Correlations
  // Feynman: with hole-particle, the loop poles lie closer to zero
  const auto hp_excitation =
    hole_particle ?
      std::optional{MBPT::lowest_hp_excitation(*wf.vHF(), n_min_core)} :
      std::nullopt;
  const auto omre =
    -std::abs(MBPT::best_omre(wf.core(), wf.valence(), true, hp_excitation));

  MBPT::FeynmanOptions options{};
  options.omre = omre;
  if (screening) {
    options.screening = MBPT::Screening::include;
  }
  if (hole_particle) {
    options.hole_particle = MBPT::HoleParticle::include;
  }

  fmt::print("Sigma sub-grid: r=({:.1e}, {:.1f})aB with {} points. "
             "[stride={}]\n",
             wf.grid().r(i0), wf.grid().r(i0 + stride * (gsize - 1)), gsize,
             stride);

  // Same disk-cache prefix as the Correlations Feynman: reads its Q*Pi*Q
  // when every option and the grids match
  auto Fyn = MBPT::FeynmanSR(wf.vHF(), i0, stride, gsize, options, n_min_core,
                             includeG, true, wf.identity());

  std::cout << "\nStructure radiation, broken up by Feynman diagram: "
            << h->name() << "\n";
  std::cout << "Reduced matrix elements <w||h||v> (" << h->units() << ")\n";

  // Values are printed as +d.[n_digits]e-xx
  const int n_digits = 3;
  const auto width = n_digits + 7;

  // (b_x) is not calculated: its column is left blank
  std::cout << "\n";
  fmt::print("{:5s} {:5s}", "w", "v");
  for (const auto label : {"a", "b", "c", "ax", "bx", "cx", "total"}) {
    fmt::print(" {:>{}s}", label, width);
  }
  fmt::print("\n");
  // Each value is printed as soon as it is known
  const auto print_value = [&](double value) {
    fmt::print(" {:+{}.{}e}", value, width, n_digits);
    std::cout << std::flush;
  };
  const auto &valence = wf.valence();
  for (std::size_t iv = 0; iv < valence.size(); ++iv) {
    // Each pair once (w >= v), unless both orders are asked for
    for (std::size_t iw = both_orders ? 0 : iv; iw < valence.size(); ++iw) {
      const auto &v = valence[iv];
      const auto &w = valence[iw];
      if (h->isZero(w.kappa(), v.kappa()))
        continue;
      const auto omega = w.en() - v.en();

      fmt::print("{:5s} {:5s}", w.shortSymbol(), v.shortSymbol());
      std::cout << std::flush;

      const auto a =
        w *
        (Fyn.dSigma_a(*h, w.kappa(), v.kappa(), v.en(), omega, rpa.get()) * v);
      print_value(a);

      const auto [dSigma_b, dSigma_c] =
        Fyn.dSigma_bc(*h, w.kappa(), v.kappa(), v.en(), omega, rpa.get());
      const auto b = w * (dSigma_b * v);
      const auto c = w * (dSigma_c * v);
      print_value(b);
      print_value(c);

      auto total = a + b + c;

      const auto [dSigma_ax, dSigma_cx] =
        Fyn.dSigma_ax_cx(*h, w.kappa(), v.kappa(), v.en(), omega, rpa.get());
      const auto ax = w * (dSigma_ax * v);
      const auto cx = w * (dSigma_cx * v);
      print_value(ax);
      fmt::print(" {:>{}s}", "", width);
      print_value(cx);
      total += ax + cx;

      print_value(total);
      fmt::print("\n");
      std::cout << std::flush;
    }
  }

  // // 2) The direct diagrams, from the functions for the E1 operator only
  // std::cout << "\nFunctions for the E1 operator only:\n";
  // fmt::print("{:5s} {:5s}", "w", "v");
  // for (const auto label : {"a", "b+c", "total"}) {
  //   fmt::print(" {:>{}s}", label, width);
  // }
  // fmt::print("\n");
  // for (const auto &w : wf.valence()) {
  //   for (const auto &v : wf.valence()) {
  //     if (h->isZero(w.kappa(), v.kappa()))
  //       continue;
  //     const auto omega = w.en() - v.en();

  //     const auto a =
  //       w * (Fyn.Sigma_SR_direct(w.kappa(), v.en(), v.kappa(), omega) * v);
  //     const auto bc =
  //       w * (Fyn.Sigma_SR_direct2(w.kappa(), v.en(), v.kappa(), omega) * v);

  //     fmt::print("{:5s} {:5s}", w.shortSymbol(), v.shortSymbol());
  //     for (const auto value : {a, bc, a + bc}) {
  //       fmt::print(" {:+{}.{}e}", value, width, n_digits);
  //     }
  //     fmt::print("\n");
  //     std::cout << std::flush;
  //   }
  // }
  // std::cout << "\n";
}

} // namespace Module
