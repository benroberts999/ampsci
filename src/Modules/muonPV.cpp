#include "Modules/muonPV.hpp"
#include "DiracODE/BoundState.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

void muonPV(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"rrms", "root-mean-square nuclear radii for muonic atom. Will "
                "use Atomic rrms by default"},
       {"type", "Nuclear density type: Fermi/ball/Gaussian/pointlike [Fermi]"},
       {"mass", "Mass of muon [206.7682830]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  auto rrms = input.get("rrms", wf.get_rrms());
  auto nuc_type = input.get<std::string>("type", "Fermi");

  int Z = wf.Znuc();
  double t = 2.3;
  double c = Nuclear::c_hdr_formula_rrms_t(rrms, t);

  // Use a different radial grid for muonic atom
  double r0 = 1.0e-7 / Z;
  double rmax = 80.0 / Z;
  std::size_t n_steps = 50000;
  double b = 1.0;
  auto radial_grid = std::make_shared<const Grid>(
      Grid{r0, rmax, n_steps, GridType::loglinear, b});
  std::cout << "Muonic grid: " << radial_grid->gridParameters() << "\n";

  const auto Vnuc = Nuclear::fermiNuclearPotential(Z, t, c, radial_grid->r());

  const double m_muon = input.get("mass", 206.7682830);

  //============================================================================

  std::cout << "\nStep 0: Check pointlike non-relativistic muon solutions:\n";
  std::cout << "(Z = " << Z << ")\n\n";

  std::cout << "        En            <r>           <1/r>\n";
  const auto V0 = Nuclear::sphericalNuclearPotential(Z, 0.0, radial_grid->r());
  double eps = 0.0;
  for (int n = 1; n < 5; ++n) {
    for (int kappa : {-1, 1, -2, 2, 3}) {

      const auto l = Angular::l_k(kappa);
      if (l >= n)
        continue;

      double e_ex = -0.5 * m_muon * Z * Z / (n * n); // "exact"
      double en0 = 0.85 * e_ex;                      // initial guess
      const double alpha0 = 1.0e-16 * PhysConst::alpha;
      const auto F =
          DiracODE::boundState(n, kappa, en0, radial_grid, V0, {}, alpha0,
                               1.0e-14, nullptr, nullptr, Z, m_muon);

      const auto e = F.en();
      const auto r_ev = F * (radial_grid->r() * F);
      const auto invr_ev = F * (radial_grid->rpow(-1) * F);

      const auto r_ex =
          1.0 / (2.0 * Z * m_muon) * (3.0 * n * n - l * (l + 1.0));
      const auto invr_ex = m_muon * Z / (n * n);

      const auto teps =
          std::max({std::abs(e / e_ex - 1.0), std::abs(r_ev / r_ex - 1.0),
                    std::abs(invr_ev / invr_ex - 1.0)});
      eps = std::max(eps, teps);

      fmt::print("{:4s}   {:.5e}   {:.5e}   {:.5e}\n", F.shortSymbol(), e, r_ev,
                 invr_ev);
      fmt::print("exact: {:.5e}   {:.5e}   {:.5e}\n", en0, r_ex, invr_ex);
    }
  }
  std::cout << "Worst epsilon: " << eps << "\n";

  //============================================================================
  std::cout << "\nStep 1: Check pointlike Relativistic muon solutions:\n";
  std::cout << "(Z = " << Z << ")\n\n";

  std::cout << "        En            <r>           <1/r>\n";
  eps = 0.0;
  for (int n = 1; n < 5; ++n) {
    for (int kappa : {-1, 1, -2, 2, 3}) {

      const auto l = Angular::l_k(kappa);
      if (l >= n)
        continue;

      const double alpha = PhysConst::alpha;

      // "exact"
      const auto e_ex = m_muon * AtomData::diracen(Z, n, kappa, alpha);

      // initial guess (force different from exact)
      const auto e0 = 0.85 * e_ex;

      const auto F =
          DiracODE::boundState(n, kappa, e0, radial_grid, V0, {}, alpha,
                               1.0e-14, nullptr, nullptr, Z, m_muon);

      const auto e = F.en();
      const auto r_ev = F * (radial_grid->r() * F);
      const auto invr_ev = F * (radial_grid->rpow(-1) * F);

      // const auto e_ex = m_muon * AtomData::diracen(Z, n, kappa, alpha);

      const auto teps = std::abs(e / e_ex - 1.0);
      eps = std::max(eps, teps);

      fmt::print("{:4s}   {:.5e}   {:.5e}   {:.5e}\n", F.shortSymbol(), e, r_ev,
                 invr_ev);
      fmt::print("exact: {:.5e}\n", e_ex);
    }
  }
  std::cout << "Worst epsilon: " << eps << "\n";

  //============================================================================
  // 2s-2p

  std::cout << "\n---------------------------------------------\n";
  std::cout << "Step 2: Scaling of PNC with Z\n\n";

  fmt::print("{:2s} {:8s} {:8s} {:8s} {:8s} {:8s} {:8s}\n", "Z", "Nrms/fm",
             "Arms/fm", "dE_bn", "Dz_an", "W_nb", "APVz");
  for (int t_Z = 1; t_Z <= 99; ++t_Z) {

    // Get alpha: allow non-relativistic calculations
    const auto alpha = wf.alpha();

    // Get A (for H, demand N=1)
    const int t_A = std::max(2, AtomData::defaultA(t_Z));

    // Look up nuclear rrms from Angeli tables
    double t_rrms = Nuclear::find_rrms(t_Z, t_A);
    // if t_rrms==0, means wasn't found. Use approx formula
    if (t_rrms == 0.0)
      t_rrms = Nuclear::approximate_r_rms(t_A, t_Z);

    // Half-density radius, c, for PNC operator
    const double t_c = Nuclear::c_hdr_formula_rrms_t(t_rrms, t);

    // Use a different radial grid for each muonic atom
    const double tr0 = 1.0e-8 / t_Z;
    const double trmax = 50.0 / t_Z;
    const std::size_t tn_steps = 10000;
    const double tb = trmax / 5.0;
    auto t_grid = std::make_shared<const Grid>(
        Grid{tr0, trmax, tn_steps, GridType::loglinear, tb});

    // Fermi (charge) distribution for Nuclear potential
    const auto t_Vnuc =
        Nuclear::fermiNuclearPotential(t_Z, t, t_c, t_grid->r());

    // Neutron number (for hW units only)
    const int t_N = t_A - t_Z;

    // Weak and dipole operators
    const auto hw = DiracOperator::PNCnsi(t_c, t, *t_grid, t_N, "i(Qw/N)e-11");
    const auto d = DiracOperator::E1(*t_grid);

    // initial energy guess
    // 0.9 roughly to account for finite nuc. size. Doesn't really matter
    const auto ea0 = 0.9 * m_muon * AtomData::diracen(t_Z, 1, -1, alpha);
    const auto es0 = 0.9 * m_muon * AtomData::diracen(t_Z, 2, -1, alpha);
    const auto ep0 = 0.9 * m_muon * AtomData::diracen(t_Z, 2, 1, alpha);

    // Get single-particle muonic wavefunctions
    const auto Fa =
        DiracODE::boundState(1, -1, ea0, t_grid, t_Vnuc, {}, alpha, 1.0e-14,
                             nullptr, nullptr, t_Z, m_muon);

    const auto Fb =
        DiracODE::boundState(2, -1, es0, t_grid, t_Vnuc, {}, alpha, 1.0e-14,
                             nullptr, nullptr, t_Z, m_muon);

    const auto Fp =
        DiracODE::boundState(2, 1, ep0, t_grid, t_Vnuc, {}, alpha, 1.0e-14,
                             nullptr, nullptr, t_Z, m_muon);

    // Check bound state solution is OK
    const auto teps = std::max(Fb.eps(), Fp.eps());
    if (teps > 1.0e-6) {
      std::cout << "Warning: didn't converge: " << teps << "\n";
      continue;
    }

    // Convert RME to z-component
    const auto Angular_z = d.rme3js(Fa.twoj(), Fp.twoj());

    // Matrix elements and energy denominators:
    const auto d_ap = Angular_z * d.reducedME(Fa, Fp);
    const auto h_pb = hw.radialIntegral(Fp, Fb);
    const auto dE_bn = Fb.en() - Fp.en();

    using namespace qip::overloads;

    // Atomic RMS radius (for p-state, just example):
    const auto rev_p = std::sqrt(Fp * (t_grid->r() * t_grid->r() * Fp));

    // PNC amplitude (just main term)
    const auto pnc = d_ap * h_pb / dE_bn;

    fmt::print("{:<2} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e}\n", t_Z, t_rrms,
               rev_p * PhysConst::aB_fm, dE_bn, d_ap, h_pb, pnc);
  }
}

} // namespace Module
