#include "GreenQED.hpp"
#include "DiracOperator/include.hpp" //For E1 operator
#include "IO/InputBlock.hpp"
#include "MBPT/Feynman.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include <complex>

namespace Module {

void GreenQED(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
    {{"r0", "Minimum r to calculate Green's fn [1.0e-3]"}, //
     {"rmax", "Maximum r to calculate Green's fn [50.0]"},
     {"num_points",
      "Number of radial points for Green's function (used for stride). [150]"},
     {"stride", "Explicitely set stride. Will over-ride num_points if set."},
     {"kappa", "Kappa value for initial test [default: -1]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  // Determine "sub" grid:
  const auto r0_target = input.get("r0", 1.0e-3);
  const auto rmax_target = input.get("rmax", 50.0);

  const std::size_t i0 = wf.grid().getIndex(r0_target);
  const std::size_t imax = wf.grid().getIndex(rmax_target);
  assert(imax > i0 && "Require rmax > r0");

  const std::size_t in_size = input.get("num_points", 150ul);
  const std::optional<std::size_t> in_stride = input.get<std::size_t>("stride");

  std::size_t size{0}, stride{0};
  if (in_stride.has_value()) {
    assert(*in_stride > 0 && "Cannot have stride of zero");
    stride = *in_stride;
    size = (imax - i0) / stride + 1;
  } else {
    assert(in_size > 1 && "Cannot have num_points <= 1");
    size = in_size;
    stride = stride = std::max(1ul, (imax - i0) / (size - 1));
  }
  assert(size > 1 && stride > 0);

  // actual r0,rmax might be slightly different, due to finite grid, stride
  const auto r0 = wf.grid().r(i0);
  const auto rmax = wf.grid().r(i0 + stride * size);
  fmt::print(
    "Grid for Green's function: {:.1e} - {:.1f} with {} points [stide = {}]\n",
    r0, rmax, size, stride);

  // We don't need QPQ (at least for now)
  bool construct_qpq = false;
  bool verbose = false;
  // I think this only affects Sigma anyway:
  bool include_G = true;

  MBPT::Feynman Fy(wf.vHF(), i0, stride, size, {}, 1, include_G, verbose, "",
                   construct_qpq);

  int kappa = input.get("kappa", -1);
  std::complex<double> en{-0.5, 0.1};

  // regular method:
  const auto g1 = Fy.green(kappa, en);

  // Use basis:
  const auto g2 = Fy.green_basis(kappa, en, wf.basis());

  // Solve Dirac eq. with complex energy directly (never checked if this works)
  const auto g3 = Fy.green_v2(kappa, en);

  // Simple quick test:
  // <v|G(e)|v> = sum_n <v|n><v|n> / (e - en) = 1 / (e - ev)

  // ampsci can do: G|v> which returns a DiracSpinor
  // G does not include integration measure, so we must include it
  // If G is complex, returns a pair: G*F = {re(GF), im(GF)}

  std::cout << "\nState, <v|G|>, expected, error\n\n";
  for (const auto &[name, g_ptr] : {std::pair{"Normal Green's fn", &g1},
                                    {"Basis Green's fn", &g2},
                                    {"Complex Schrodinger Green's fn", &g3}}) {

    std::cout << "For " << name << "\n";
    const auto &gt = *g_ptr;
    for (const auto &v : wf.valence()) {

      if (v.kappa() != kappa)
        continue;

      const auto [re_Gv, im_Gv] = gt.drj() * v;
      const std::complex<double> value = {v * re_Gv, v * im_Gv};
      const auto expected = 1.0 / (en - v.en());

      auto eps = std::abs(value - expected) / std::abs(expected);

      std::cout << v << " " << value << " " << expected << " " << eps << "\n";
    }
    std::cout << "\n";
  }
}

} // namespace Module
