#include "Modules/muonPV.hpp"
#include "DiracODE/BoundState.hpp"
#include "DiracODE/InhomogenousGreens.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/FGRadPot.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"

namespace Module {

void muonPV(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"states", "States to calculate [4sp]"},
       {"t", "t: skin thickness parameter in fm [2.3]"},
       {"type", "Nuclear density type: Fermi/ball/Gaussian/pointlike [Fermi]"},
       {"mass", "Mass of muon (in units of electron mass) [206.7682830]"},
       {"N_steps", "Number of grid points to use for Muon solutions (is "
                   "different from neutral case). [10000]"},
       {"write_wf", "Write wavefunctions (f and g) to file [false]"},
       {"Uehling", "Include Uehling potential [false]"},
       {"Screening",
        "Account for electron screening (direct potental)? [false]"},
       {"Z_scaling", "Perform basic 1s-2s calculations for Z=1-99 [true]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const double m_muon = input.get("mass", 206.7682830);
  const auto rrms = input.get("rrms", wf.get_rrms());
  const double t = input.get("t", 2.3);
  const auto nuc_type = input.get<std::string>("type", "Fermi");
  const auto write_wf = input.get("write_wf", false);
  // nb: Uehling/screening not used in first few tests/checks
  const auto Uehling = input.get("Uehling", false);
  const auto elec_screening = input.get("Screening", false);
  const auto Z_scaling = input.get("Z_scaling", true);

  const auto states_str = input.get("states", std::string{"4sp"});

  const int Z = wf.Znuc();
  const int A = wf.Anuc();

  const auto states = AtomData::listOfStates_nk(states_str);

  // Use a different radial grid for muonic atom
  const double r0 = 1.0e-5 / Z / m_muon;
  const double rmax = 1.0e5 / Z / m_muon;
  const std::size_t n_steps = input.get("N_steps", 10000ul);
  const double b = rmax / 10.0;

  auto radial_grid = std::make_shared<const Grid>(
      Grid{r0, rmax, n_steps, GridType::loglinear, b});
  std::cout << "Muonic grid:\n" << radial_grid->gridParameters() << "\n";

  std::cout << "\n"
            << "Running for 'muon' with mass:\nM = " << m_muon
            << " m_e = " << m_muon * PhysConst::m_e_MeV << " MeV\n";

  //============================================================================

  std::cout << "\n-------------------------\n"
            << "First, do some checks to ensure muonic Dirac solver working\n";

  std::cout << "\nStep 0: Check pointlike non-relativistic muon solutions:\n";
  std::cout << "(Z = " << Z << ")\n";
  std::cout << "(M = " << m_muon << ")\n\n";

  std::cout << "        En            <r>           <1/r>\n";
  const auto V0 = Nuclear::sphericalNuclearPotential(Z, 0.0, radial_grid->r());
  double eps = 0.0;
  for (const auto [n, kappa, x_en] : states) {
    const auto l = Angular::l_k(kappa);

    double e_ex = -0.5 * m_muon * Z * Z / (n * n); // "exact"
    double en0 = 0.85 * e_ex;                      // initial guess
    const double alpha0 = 1.0e-16 * PhysConst::alpha;
    const auto F =
        DiracODE::boundState(n, kappa, en0, radial_grid, V0, {}, alpha0,
                             1.0e-14, nullptr, nullptr, Z, m_muon);
    if (F.eps() > 1.0e-9) {
      std::cout << "# Warning: didn't converge\n";
    }

    const auto e = F.en();
    const auto r_ev = F * (radial_grid->r() * F);
    const auto invr_ev = F * (radial_grid->rpow(-1) * F);

    const auto r_ex = 1.0 / (2.0 * Z * m_muon) * (3.0 * n * n - l * (l + 1.0));
    const auto invr_ex = m_muon * Z / (n * n);

    const auto teps =
        std::max({std::abs(e / e_ex - 1.0), std::abs(r_ev / r_ex - 1.0),
                  std::abs(invr_ev / invr_ex - 1.0)});
    eps = std::max(eps, teps);

    fmt::print("{:4s}   {:.5e}   {:.5e}   {:.5e}  [{:.1e}]\n", F.shortSymbol(),
               e, r_ev, invr_ev, eps);
    fmt::print("exact: {:.5e}   {:.5e}   {:.5e}\n", e_ex, r_ex, invr_ex);
  }
  std::cout << "Worst error: " << eps << "\n";

  //============================================================================
  std::cout << "\nStep 1: Check pointlike Relativistic muon solutions:\n";
  std::cout << "(Z = " << Z << ")\n";
  std::cout << "(M = " << m_muon << ")\n\n";

  std::cout << "        En            <r>           <1/r>\n";
  eps = 0.0;
  for (const auto [n, kappa, x_en] : states) {

    const double alpha = PhysConst::alpha;

    // "exact"
    const auto e_ex = m_muon * AtomData::diracen(Z, n, kappa, alpha);

    // initial guess (force different from exact)
    const auto e0 = 0.85 * e_ex;

    const auto F =
        DiracODE::boundState(n, kappa, e0, radial_grid, V0, {}, alpha, 1.0e-14,
                             nullptr, nullptr, Z, m_muon);

    if (F.eps() > 1.0e-9) {
      std::cout << "# Warning: didn't converge\n";
    }

    const auto e = F.en();
    const auto r_ev = F * (radial_grid->r() * F);
    const auto invr_ev = F * (radial_grid->rpow(-1) * F);

    const auto teps = std::abs(e / e_ex - 1.0);
    eps = std::max(eps, teps);

    fmt::print("{:4s}   {:.5e}   {:.5e}   {:.5e}  [{:.1e}]\n", F.shortSymbol(),
               e, r_ev, invr_ev, eps);
    fmt::print("exact: {:.5e}\n", e_ex);
  }
  std::cout << "Worst error: " << eps << "\n";

  //============================================================================
  {
    using namespace qip::overloads;
    std::cout << "\n---------------------------------------------\n";
    std::cout << "Step 2: Muonic " << AtomData::atomicSymbol(Z) << "\n";
    std::cout << "Z = " << Z << ", A=" << A << "\n";
    std::cout << "Using " << nuc_type << " nuclear potential\n";

    // Nuclear radius, in atomic units:
    const double Rn_au = std::sqrt(5.0 / 3) * rrms / PhysConst::aB_fm;

    std::cout << "Nulcear: rrms = " << rrms
              << " (Rn = " << Rn_au * PhysConst::aB_fm << "), t = " << t
              << "  fm\n";

    if (Uehling) {
      std::cout << "Including Uehling potential (w/ FNS)\n";
    }
    if (elec_screening && Z > 1) {
      std::cout << "Including electron screening potential (V -> V+Vdir)\n";
    }

    std::cout << "M = " << m_muon << "\n\n";

    const auto alpha = wf.alpha();

    auto Vnuc = Nuclear::formPotential(
        Nuclear::Nucleus(Z, A, nuc_type, rrms, t), radial_grid->r());

    if (Uehling) {
      const auto V_ueh = [Z, Rn_au](double r) {
        return FGRP::V_Uehling(Z, r, Rn_au);
      };
      Vnuc += qip::apply_to(V_ueh, radial_grid->r());
    }
    if (elec_screening && Z > 1) {
      // Interpolate neutral direct potential onto muonic grid:
      Vnuc += Interpolator::interpolate(wf.grid().r(), wf.vHF()->vdir(),
                                        radial_grid->r());
    }

    std::cout << "Energies:\n";
    std::cout << "nk    Rinf  eps    R_rms (a0)    E (au)            E "
                 "(keV)\n";

    std::vector<DiracSpinor> muon_Fs;
    for (const auto [n, kappa, x_en] : states) {
      const auto e0 = m_muon * AtomData::diracen(Z, n, kappa, alpha);
      const auto &Fnk = muon_Fs.emplace_back(
          DiracODE::boundState(n, kappa, e0, radial_grid, Vnuc, {}, alpha,
                               1.0e-14, nullptr, nullptr, Z, m_muon));

      const auto R_rms =
          std::sqrt(Fnk * (radial_grid->r() * radial_grid->r() * Fnk));

      fmt::print("{:4s} {:5.2f}  {:5.0e}  {:.5e}  {:.9e}  {:.9e}\n",
                 Fnk.shortSymbol(), Fnk.rinf(), Fnk.eps(), R_rms, Fnk.en(),
                 Fnk.en() * PhysConst::Hartree_eV / 1.0e3);
    }

    // PNC and E1 matrix elements:
    const double c = Nuclear::c_hdr_formula_rrms_t(rrms, t);
    const int N = A - Z;
    const auto hw = DiracOperator::PNCnsi(c, t, *radial_grid, N, "i(Qw/N)e-11");
    const auto d = DiracOperator::E1(*radial_grid);

    std::cout << "\nMatrix elements:\n";
    std::cout << "e    o      E_e-E_o (au)    <e|dz|o> (ea0)   <e|hw|o> "
                 "(i[Qw/N]e-11)\n";
    for (const auto &Fe : muon_Fs) {
      for (const auto &Fo : muon_Fs) {

        if (Fe.parity() != 1)
          continue;
        if (Fo.parity() != -1)
          continue;
        if (d.isZero(Fe, Fo) && hw.isZero(Fe, Fo))
          continue;

        // Convert RME to z-component
        const auto Angular_z = d.rme3js(Fe.twoj(), Fo.twoj());
        const auto dz_eo = Angular_z * d.reducedME(Fe, Fo);

        const auto hw_eo = hw.radialIntegral(Fe, Fo); //scalar operator

        fmt::print("{:4s} {:4s}  {:14.7e}  {:14.7e}    ", Fe.shortSymbol(),
                   Fo.shortSymbol(), Fe.en() - Fo.en(), dz_eo);
        if (hw.isZero(Fe, Fo)) {
          fmt::print(" {:<13}\n", 0);
        } else {
          fmt::print("{:14.7e}\n", hw_eo);
        }
      }
    }

    std::cout << "\n\n----------------------------------------------------\n";
    std::cout << "Calculate full PNC using Mixed States method\n";
    std::cout << "Solve: (H - e_a)δFa = -h*Fa\n";
    std::cout << "Then PNC = <Fa|d|δFb> + <δFa|d|Fb>\n";
    std::cout << "* [units: i(Qw/N)e-11]\n";

    // Calculate PNC using mixed states
    for (const auto &Fa : muon_Fs) {
      // do for s states only
      if (Fa.kappa() != -1)
        continue;
      // Find correspinding (n+1)s state:
      const auto pFb = DiracSpinor::find(Fa.n() + 1, Fa.kappa(), muon_Fs);
      // "main" intermediate state:
      const auto pF2p = DiracSpinor::find(Fa.n() + 1, 1, muon_Fs);

      // ensure it exists:
      if (pFb == nullptr && pF2p == nullptr)
        continue;
      const auto &Fb = *pFb;
      const auto &F2p = *pF2p;

      std::cout << "\n" << Fa << " " << Fb << "\n";

      // Solve inhomogenous Dirac equation:
      // (H_0 + v - e)dF = -h*F
      const auto kappa_n = -1 * Fa.kappa();
      const auto hFa_dag = -1 * hw.radial_rhs(kappa_n, Fa); // <Fa|h = -h|Fa>
      const auto hFb = hw.radial_rhs(kappa_n, Fb);          // h|Fb

      const auto dFa_dag =
          DiracODE::solve_inhomog(kappa_n, Fa.en(), Vnuc, {}, alpha,
                                  -1 * hFa_dag, nullptr, nullptr, Z, m_muon);
      const auto dFb =
          DiracODE::solve_inhomog(kappa_n, Fb.en(), Vnuc, {}, alpha, -1 * hFb,
                                  nullptr, nullptr, Z, m_muon);

      // angular factors (for dipole z component)
      const auto Angular_d_a = d.rme3js(Fa.twoj(), dFb.twoj());
      const auto Angular_d_b = d.rme3js(Fb.twoj(), dFa_dag.twoj());

      // Full PNC from Mixed states: dd_a + dd_b = <δFa|d|Fb> + <Fa|d|δFb>
      const auto dd_a = Angular_d_a * d.reducedME(Fa, dFb);
      const auto dd_b = Angular_d_b * d.reducedME(dFa_dag, Fb);

      // Main term (from mixed states, by forcing orthogonality):
      const auto pnc_main_MS = Angular_d_a * d.reducedME(Fa, (dFb * F2p) * F2p);

      // // Numerical test:
      // // Main term (direct calculation):
      // const auto pnc_main_direct = Angular_d_a * d.reducedME(Fa, F2p) *
      //                              hw.radialIntegral(F2p, Fb) /
      //                              (Fb.en() - F2p.en());
      // std::cout << "Numerical test of MS method:\nmain = "
      //           << fmt::format("<{}|dz|{}><{}|hw|{}>/dE = ", Fa.shortSymbol(),
      //                          F2p.shortSymbol(), F2p.shortSymbol(),
      //                          Fb.shortSymbol())
      //           << pnc_main_direct << "\n"
      //           << fmt::format("     = <{}|dz|{}><{}|δ{}>", Fa.shortSymbol(),
      //                          F2p.shortSymbol(), F2p.shortSymbol(),
      //                          Fb.shortSymbol())
      //           << "      = " << pnc_main_MS << "\n";
      // // std::cout << "Main (direct calc): " << pnc_main_direct << "\n";
      // // std::cout << "Main (MS + orthog): " << pnc_main_MS << "\n";
      // const auto err = std::abs(pnc_main_MS / pnc_main_direct - 1.0);
      // std::cout << "Error: " << err << "\n";

      const auto pnc = dd_a + dd_b;
      const auto tail = pnc - pnc_main_MS;
      std::cout << "\n";
      fmt::print("<{}s|dz|δ{}s> = {:12.5e}\n", Fa.n(), Fb.n(), dd_a);
      fmt::print("<δ{}s|dz|{}s> = {:12.5e}\n", Fa.n(), Fb.n(), dd_b);
      fmt::print("main term   = {:12.5e}   ({:.1f}%)\n", pnc_main_MS,
                 pnc_main_MS / pnc * 100.0);
      fmt::print("rest        = {:12.5e}\n", tail);
      fmt::print("Total       = {:12.5e}  i(Qw/N)e-11\n", pnc);
    }

    // print wavefunctions:
    if (write_wf) {
      std::ofstream ofile{fmt::format("muon_fg_{}.txt", Z)};
      ofile << "R ";
      for (const auto &Fn : muon_Fs) {
        ofile << "f_{" << Fn.shortSymbol() << "} g_{" << Fn.shortSymbol()
              << "} ";
      }
      ofile << "\n";
      for (auto i = 0ul; i < radial_grid->size(); ++i) {
        ofile << radial_grid->r(i) << " ";
        for (const auto &Fn : muon_Fs) {
          ofile << Fn.f(i) << " " << Fn.g(i) << " ";
        }
        ofile << "\n";
      }
    }
  }

  //============================================================================
  // 2s-2p

  if (Z_scaling) {
    std::cout << "\n---------------------------------------------\n";
    std::cout << "Step 3: Scaling of 1s-2s muonic PNC with Z\n";
    std::cout << "Using " << nuc_type << " nuclear potential\n";
    if (Uehling) {
      std::cout << "Including Uehling potential (w/ FNS)\n";
    }
    if (elec_screening && Z > 1) {
      std::cout << "nb: screening _not_ included in this section\n";
    }
    std::cout << "(M = " << m_muon << ")\n\n";

    fmt::print(
        "{:2s} {:3s}  {:8s} {:8s} {:8s} {:8s} {:8s} {:8s} {:8s} {:8s} {:8s}\n",
        "Z", "N", "E1s", "Nrms_fm", "Arms_fm", "dE_bn", "Dz_an", "W_nb",
        "W0_nb", "sr_nb", "APVz");
    for (int t_Z = 1; t_Z <= 99; ++t_Z) {

      // Get alpha: allow non-relativistic calculations
      const auto alpha = wf.alpha();

      // Get A (for H, demand N=1)
      const int t_A = std::max(2, AtomData::defaultA(t_Z));

      // Look up nuclear rrms from Angeli tables
      double t_rrms = Nuclear::find_rrms(t_Z, t_A);
      // if t_rrms==0, means wasn't found in tables. Use approx formula
      if (t_rrms == 0.0) {
        t_rrms = Nuclear::approximate_r_rms(t_A, t_Z);
      }

      // Nuclear radius, in atomic units:
      const double t_Rn_au = std::sqrt(5.0 / 3) * t_rrms / PhysConst::aB_fm;

      // Use a different radial grid for each muonic atom
      const double tr0 = 1.0e-5 / t_Z / m_muon;
      const double trmax = 1.0e5 / t_Z / m_muon;
      const std::size_t tn_steps = n_steps;
      const double tb = trmax / 10.0;

      auto t_grid = std::make_shared<const Grid>(
          Grid{tr0, trmax, tn_steps, GridType::loglinear, tb});

      // Form nuclear potential:
      auto t_Vnuc = Nuclear::formPotential(
          Nuclear::Nucleus(t_Z, t_A, nuc_type, t_rrms, t), t_grid->r());

      // nb: Also option to include Uehling?
      if (Uehling) {
        using namespace qip::overloads;
        const auto V_ueh = [t_Z, t_Rn_au](double r) {
          return FGRP::V_Uehling(t_Z, r, t_Rn_au);
        };
        t_Vnuc += qip::apply_to(V_ueh, t_grid->r());
      }

      // Neutron number (for hW units only)
      const int t_N = t_A - t_Z;

      // Half-density radius, c, for PNC operator
      const double t_c = Nuclear::c_hdr_formula_rrms_t(t_rrms, t);

      // Weak and dipole operators
      const auto hw =
          DiracOperator::PNCnsi(t_c, t, *t_grid, t_N, "i(Qw/N)e-11");
      const auto d = DiracOperator::E1(*t_grid);
      // "Special" weak operator: constant rho(r), assumed to be larger than wavefunction
      const auto hw0 = DiracOperator::PNCnsi_const(t_Rn_au, t_N, "i(Qw/N)e-11");
      // sigma.r
      const auto sr = DiracOperator::sigma_r(*t_grid);

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
        std::cout << "# Warning: didn't converge: Z=" << t_Z << ", eps=" << teps
                  << "\n";
        continue;
      }

      // Convert RME to z-component
      const auto Angular_z = d.rme3js(Fa.twoj(), Fp.twoj());

      // Matrix elements and energy denominators:
      const auto d_ap = Angular_z * d.reducedME(Fa, Fp);
      const auto h_pb = hw.radialIntegral(Fp, Fb);
      const auto dE_bn = Fb.en() - Fp.en();

      // std::cout << Fb.en() << "\n" << Fp.en() << "\n" << dE_bn << "\n";
      // std::cin.get();

      const auto h0_pb = hw0.radialIntegral(Fp, Fb);
      const auto sr_bp = sr.radialIntegral(Fp, Fb);

      using namespace qip::overloads;

      // Atomic RMS radius (for p-state, just example):
      const auto rev_p = std::sqrt(Fp * (t_grid->r() * t_grid->r() * Fp));

      // PNC amplitude (just main term)
      const auto pnc = d_ap * h_pb / dE_bn;

      fmt::print("{:<2} {:<3} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e} "
                 "{:.2e} {:.2e}\n",
                 t_Z, t_N, Fa.en(), t_rrms, rev_p * PhysConst::aB_fm, dE_bn,
                 d_ap, h_pb, h0_pb, sr_bp, pnc);
    }
  }
}

} // namespace Module
