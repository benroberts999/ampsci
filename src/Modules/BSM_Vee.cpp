#include "Modules/BSM_Vee.hpp"
#include "Angular/Wigner369j.hpp"
#include "CI/CI_Integrals.hpp"
#include "DiracOperator/GenerateOperator.hpp"
#include "DiracOperator/Operators/Ek.hpp"
#include "DiracOperator/Operators/RadialF.hpp"
#include "DiracOperator/Operators/Vee.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "IO/InputBlock.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Modules/Modules.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Potentials/BSM_Vee_Potentials.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "ampsci/ampsci.hpp"
#include "fmt/format.hpp"
#include <cmath>
#include <gsl/gsl_sf.h>
#include <iostream>
#include <omp.h>

namespace Module {

namespace {
const Register r_BSM_Vee{
  "BSM_Vee", "BSM_Vee: New electron-electron interaction", &BSM_Vee};
} // namespace

void BSM_Vee(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check(
    {{"", "Introduces a new electron-electron "
          "interaction."},
     {"type", "'sp' (scalar-pseudoscalar), 'va' (vector-axial "
              "vector), 'ss' (scalar-scalar), "
              "'vv' (vector-vector) ['sp']"},
     {"tdhf", "Include TDHF calcs for 'sp'? [false]"},
     {"omega", "ω (for tdhf) [0]"},
     {"contact", "Consider μ->infty, i.e. a contact force [false]"},
     {"min_mu", "Minimum mediator mass to consider [1e-6]"},
     {"max_mu", "Maximum mediator mass to consider [20]"},
     {"N_mu", "Number of masses to consider [100]"},
     {"v_state", "Valence state (ket, |v>) to consider [ground]"},
     {"w_state", "Valence state (bra, <w|) to consider. [ground+1] "
                 "if type=va, otherwise [ground]."},
     {"operator",
      "Operator to get matrix element of; otherwise, all applicable."},
     {"n", "Principal quantum number for state [ground]"},
     {"kappa", "Kappa for state [ground]"},
     {"A2", "Second isotope's mass (for 'ss' or 'vv') [A+5]"},
     {"g0", "Include the gamma-0 term on both electrons? [true]"},
     {"test", "Run module testing [false]"},
     {"ci_basis", "ci_basis, to be removed!"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const bool test = input.get<bool>("test", false);

  if (test) {
    // sp_tdhf(input, wf);
    // return;
    std::cout << "\n***Using V_SP operator (test mode)***\n";
    // find_fierz(input, wf);
    test_CI(input, wf);
    return;
  }

  const auto e_handler = gsl_set_error_handler_off();

  const std::string int_type = input.get<std::string>("type", "sp");

  const double min_mu = input.get<double>("min_mu", 1.0e-4);
  const double max_mu = input.get<double>("max_mu", 1.0e4);
  const int N_mu = input.get<int>("N_mu", 100);
  const bool contact = input.get<bool>("contact", false);
  const bool tdhf = input.get<bool>("tdhf", false);
  const double omega = input.get<double>("omega", 0.0);
  const std::string type = input.get<std::string>("type", "");
  const DiracSpinor Fv = *wf.getState(
    input.get<std::string>("v_state", wf.valence()[0].shortSymbol()));

  const DiracSpinor default_Fw = type == "va" ? wf.valence()[1] : Fv;

  const DiracSpinor Fw =
    *wf.getState(input.get<std::string>("w_state", default_Fw.symbol()));

  const std::string op = input.get<std::string>("operator", "");

  if (op != "") {
    matrix_elements(wf, Fv, Fw, op, int_type, contact, tdhf, omega, min_mu,
                    max_mu, N_mu);
  } else if (int_type == "sp" || int_type == "va") {
    matrix_elements(wf, Fv, Fw, "D", int_type, contact, tdhf, omega, min_mu,
                    max_mu, N_mu);
  } else if (int_type == "ss" || int_type == "vv") {
    matrix_elements(wf, Fv, Fw, "V", int_type, contact, tdhf, omega, min_mu,
                    max_mu, N_mu);
    // V_energy_shift(input, wf);
    // matrix_elements(input, wf); // for now, testing
  } else {
    std::cout << "\nERROR: 'type = " << int_type
              << " is not valid. Ensure one of \n - 'sp' "
                 "(scalar-pseudoscalar)\n - 'ss' (scalar-scalar)\n - 'vv' "
                 "(vector-vector)\n is provided.";
  }

  // Cleanup
  gsl_set_error_handler(e_handler);
}

void test_CI(const IO::InputBlock &input, const Wavefunction &wf) {
  std::cout << "\nTesting CI;";
  std::cout << "\n\nAttempting matrix elements:\n";

  const auto basis_string = input.get<std::string>("ci_basis", "");
  const auto contact = input.get<bool>("contact", false);

  const double min_mu = input.get<double>("min_mu", 1e-4);
  const double max_mu = input.get<double>("max_mu", 1e4);
  const double N_mu = input.get<double>("N_mu", 100);

  std::vector<double> mus;

  // Populate mus
  for (double log_mu = log(min_mu); log_mu <= log(max_mu);
       log_mu += std::abs(log(max_mu) - log(min_mu)) * (1.0 / N_mu)) {
    mus.push_back(std::exp(log_mu));
  }
  const auto actual_N_mu = mus.size();

  std::vector<double> Ds(actual_N_mu);

#pragma omp parallel default(none)                                             \
  shared(actual_N_mu, wf, basis_string, contact, mus, Ds, std::cout)

  {

#pragma omp for
    for (int i = 0; i < actual_N_mu; ++i) {
      const double mu = mus[i];

      const auto D = CI_edm(wf, basis_string, mu, contact, 2, -1, 0, false);

#pragma omp critical
      {
        Ds[i] = D;
      }
    }
  }

  // print
  std::cout << "\n\n        mu          D\n";
  for (int i = 0; i < mus.size(); ++i) {
    fmt::print("{:10.7f} {:10.7f}\n", mus[i], Ds[i]);
  }
}

double CI_edm(const Wavefunction &wf, const std::string basis_string,
              const double mu, const bool contact, const int J_V,
              const int pi_V, const int i_V, const bool verbose) {

  const std::vector<DiracSpinor> ci_basis =
    CI::basis_subset(wf.basis(), basis_string, wf.coreConfiguration());

  std::string contact_str = contact ? "true" : "false";

  const auto mu_str = std::to_string(mu);
  const auto VeeOptions = "contact=" + contact_str + ";mu=" + mu_str + ";";

  const auto VeeOp = DiracOperator::Vee(wf.core(), contact, mu, "sp", false)
                       .generate(IO::InputBlock("Vee", VeeOptions), wf);

  const auto Vee_rpa = ExternalField::make_rpa("diagram", VeeOp.get(), wf.vHF(),
                                               true, wf.basis(), wf.identity());

  Vee_rpa->solve_core(0.0, 300);

  const auto d =
    DiracOperator::E1(wf.grid()).generate(IO::InputBlock("E1", ""), wf);

  const auto d_rpa = ExternalField::make_rpa("TDHF", d.get(), wf.vHF(), true,
                                             wf.basis(), wf.identity());
  d_rpa->solve_core(0.0, 300);

  if (verbose) {
    std::cout << "\nCalculating matrix element tables.." << std::flush;
  }
  auto Vee_table =
    ExternalField::me_table(ci_basis, VeeOp.get(), Vee_rpa.get(), nullptr);
  auto d_table =
    ExternalField::me_table(ci_basis, d.get(), d_rpa.get(), nullptr);

  if (verbose) {
    std::cout << "..E1 done.\n\n" << std::flush;
  }

  const auto wf_V = *wf.CIwf(J_V, pi_V);
  const auto E_V = wf_V.energy(i_V);

  const auto n_print = wf_V.num_solutions();

  if (verbose) {
    std::cout << "       V        |        N        | \n";
    std::cout << "J   π  #  conf. | J   π  #  conf. |  <V|d|N>    <N|Vee|V>\n";
  }

  double D = 0.0;

  for (auto wf_N : wf.CIwfs()) {
    for (int iN = 0; iN < n_print; ++iN) {
      if (wf_N.parity() == wf_V.parity()) {
        continue;
      }
      const auto V_rme = CI::ReducedME(wf_N, iN, wf_V, i_V, Vee_table,
                                       VeeOp->rank(), VeeOp->parity());

      const double V_me = V_rme * VeeOp->rme3js(wf_N.twoJ(), wf_V.twoJ(), 2);

      // For rme3js - is it acceptable to leave twomb=1? Should it be 2?
      // mb is projection of Jb=2, mb takes -2,-1,0,1,2 and two_mb -4,-2,0,2,4 - never 1

      const auto d_rme =
        CI::ReducedME(wf_V, 0, wf_N, iN, d_table, d->rank(), d->parity());
      const double d_me = d_rme * d->rme3js(wf_V.twoJ(), wf_N.twoJ(), 2);

      const double dE = wf_N.energy(iN) - wf_V.energy(i_V);

      if (verbose) {
        fmt::print("{} {:3}  {}  {:5s} | {} {:3}  {}  {:5s} | {:10.7f}  "
                   "{:10.7f}\n",
                   0.5 * wf_V.twoJ(), wf_V.parity(), i_V, wf_V.CSF(0).config(),
                   0.5 * wf_N.twoJ(), wf_N.parity(), iN, wf_N.CSF(iN).config(),
                   d_me, V_me);
      }

      D += 2 * d_me * V_me / dE;
    }
  }
  return D / PhysConst::alpha;
}

void find_fierz(const IO::InputBlock &input, const Wavefunction &wf) {

  // Check some random angular identities

  // const int kappa_v = -1;
  // const int kappa_n = 1;

  const auto Fv = wf.valence()[0];
  const auto Fn = wf.valence()[2];
  const auto Fa = wf.core()[0];

  if (Fn.kappa() != -Fv.kappa()) {
    throw "κn != -κv";
  }

  std::cout << "\nUsing v = " << Fv.symbol() << ", n = " << Fn.symbol()
            << ", a = " << Fa.symbol();

  std::cout << "\n\nCHECKING ANGULAR CONTRIBUTIONS\n";

  const auto A_VA_dir = std::sqrt(Fa.twojp1() * (1.0 / Fv.twojp1())) *
                        Angular::Ck_kk(0, Fn.kappa(), Fv.kappa()) *
                        Angular::Ck_kk(0, Fa.kappa(), -Fa.kappa());

  const auto A_AV_dir = std::sqrt(Fa.twojp1() * (1.0 / Fv.twojp1())) *
                        Angular::Ck_kk(0, Fn.kappa(), -Fv.kappa()) *
                        Angular::Ck_kk(0, Fa.kappa(), Fa.kappa());

  std::cout << "\nVA direct: " << A_VA_dir;
  std::cout << "\nAV direct: " << A_AV_dir;

  std::cout << "\n\nSumming for VA exchange: \n";

  double A_VA_exch = 0;
  for (int k = 0; k < 6; ++k) {
    if (k > 0) {
      std::cout << " + ";
    }
    const auto factor = (2 * k + 1) *
                        Angular::neg1pow_2(Fv.twoj() - Fa.twoj()) *
                        (1.0 / Fv.twojp1());
    const auto contribution = factor *
                              Angular::Ck_kk(k, Fn.kappa(), Fa.kappa()) *
                              Angular::Ck_kk(k, Fa.kappa(), -Fv.kappa());
    std::cout << contribution;

    A_VA_exch += contribution;
  }
  std::cout << "\nVA exchange: " << A_VA_exch;

  std::cout << "\n\nSumming for AV exchange: \n";

  double A_AV_exch = 0;
  for (int k = 0; k < 6; ++k) {
    if (k > 0) {
      std::cout << " + ";
    }
    const auto factor = (2 * k + 1) *
                        Angular::neg1pow_2(Fv.twoj() - Fa.twoj()) *
                        (1.0 / Fv.twojp1());
    const auto contribution = factor *
                              Angular::Ck_kk(k, Fn.kappa(), -Fa.kappa()) *
                              Angular::Ck_kk(k, Fa.kappa(), Fv.kappa());
    std::cout << contribution;

    A_AV_exch += contribution;
  }
  std::cout << "\nAV exchange: " << A_AV_exch;

  // std::cout << "\n\nVA + AV   direct: " << A_VA_dir + A_AV_dir;
  // std::cout << "\nVA + AV exchange: " << A_VA_exch + A_AV_exch;

  // Check the radial
  std::cout << "\n\nCHECKING RADIAL CONTRIBUTIONS\n";

  // Do it manually? Seems to work out with what I've got.

  const auto R_VA_nava = BSM_Vee::R_AV_abcd(Fa, Fn, Fa, Fv);
  const auto R_AV_nava = BSM_Vee::R_AV_abcd(Fn, Fa, Fv, Fa);
  const auto R_AV_naav = BSM_Vee::R_AV_abcd(Fn, Fa, Fa, Fv);
  const auto R_VA_naav = BSM_Vee::R_AV_abcd(Fa, Fn, Fv, Fa);

  std::cout << "\n  VA direct: " << R_VA_nava << "\n  AV direct: " << R_AV_nava
            << "\nVA exchange: " << R_VA_naav << "\nAV exchange: " << R_AV_naav
            << "\n\nVA + AV   direct: " << R_VA_nava + R_AV_nava
            << "\nVA + AV exchange: " << R_VA_naav + R_AV_naav;
}

/*
void ee_isotope_shift(const std::string int_type, const IO::InputBlock &input,
                      const Wavefunction &wf) {

  // Create second wavefunction
  int A2 = input.get<int>("A2", wf.Anuc() + 5);

  // Create new input block, warn that modules after Vee will not run for wf2.
  auto new_input =
    IO::InputBlock("ampsci", input.path(), std::fstream(input.path()));

  new_input.merge("Atom{A = " + std::to_string(A2) + ";}");

  // Currently just adds new_A without removing previous. OK because it reads the last,
  // but would be safer to remove original

  // Remove Vee module
  const auto blocks_copy = new_input.blocks();
  auto flag = false;

  for (const auto block : blocks_copy) {
    auto name = block.name();

    if (name == "Module::Vee") {
      new_input.remove_block(name);
      flag = true;
    } else if (flag == true) {
      new_input.remove_block(name);
      std::cout << "\nWARNING: removing '" << name << "' block for wf2.";
    }
  }

  // Create second wavefunction
  std::cout << "\n\nCreating wavefunction for A = " << A2 << ".\n";

  // Create second wavefunction
  const auto wf2 = ampsci(new_input);

  std::cout
    << "\nCalculating energy shifts from scalar-scalar interactions in\n";
  std::cout << wf.atom() << " and \n" << wf2.atom() << "\n";

  // Get ground state
  int i_ground = 0;
  double E_ground = wf.valence()[i_ground].en();
  const double min_mu = input.get<double>("min_mu", 1.0e-4);
  const double max_mu = input.get<double>("max_mu", 1.0e4);
  const double N_mu = input.get<double>("N_mu", 100.0);

  // Find the ground state of the given valence states
  for (int i = 1; i < wf.valence().size(); ++i) {
    const double E_test = wf.valence()[i].en();

    if (E_test < E_ground) {
      E_ground = E_test;
      i_ground = i;
    }
  }

  const auto n_ground = wf.valence()[i_ground].n();
  const auto kappa_ground = wf.valence()[i_ground].kappa();

  const int v_n = input.get<int>("n", n_ground);
  const int v_kappa = input.get<int>("kappa", kappa_ground);

  const auto Fv = *wf.getState(v_n, v_kappa);
  const auto Fv2 = *wf2.getState(v_n, v_kappa);

  //const auto Fv = wf.valence()[0];

  // Loop through mus
  std::cout
    << "\nEnergy shift due to a new scalar-scalar "
       "electron interaction (to first order).\nCalculating for the state "
    << Fv.symbol() << " with mediator mass μ."
    << "\n         μ     E1 (au)    dE1 (au)     E2 (au)    dE2 (au)     I"
       "S (au)    IS (MHz)\n";

  for (double log_mu = log(min_mu); log_mu < log(max_mu);
       log_mu += std::abs(log(max_mu) - log(min_mu)) / N_mu) {

    const auto mu = std::exp(log_mu);

    double dE1 = dE(mu, int_type, wf.core(), Fv);
    double dE2 = dE(mu, int_type, wf2.core(), Fv2);

    const auto IS = dE1 - dE2;

    fmt::print(
      "{:10.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4f}\n", mu,
      Fv.en(), dE1, Fv2.en(), dE2, IS, IS * PhysConst::Hartree_MHz);
  }
}

double dE(const double mu, const std::string int_type,
          const std::vector<DiracSpinor> &core, const DiracSpinor Fv) {
  double dE_val = 0;
  for (auto Fa : core) {

    const auto R_vava = Rk_abcd(0, mu, Fv, Fa, Fv, Fa, int_type);
    const auto A_vava = Fa.twoj();

    // Sum over |ja-jv| <= k <= ja+jv
    auto R_vaav = 0.0;
    auto A_vaav = 0.0;

    for (int twok = std::abs(Fa.twoj() - Fv.twoj());
         twok <= Fa.twoj() + Fv.twoj(); twok += 2) {
      if ((Fa.twoj() + Fv.twoj() + twok) % 4 == 0) {
        const int k = round(twok / 2);
        R_vaav = Rk_abcd(k, mu, Fv, Fa, Fa, Fv, int_type);
        A_vaav = (2.0 * k + 1) * Angular::Ck_kk(k, Fv.kappa(), Fa.kappa()) *
                 Angular::Ck_kk(k, Fa.kappa(), Fv.kappa());
      }
      const double jv = Fv.twoj() * 0.5;
      A_vaav *= std::pow(-1.0, (Fv.twoj() - Fa.twoj()) * 0.5) / (2.0 * jv + 1);

      const auto u_vava = R_vava * A_vava;
      const auto u_vaav = R_vaav * A_vaav;

      dE_val += u_vava - u_vaav;
    }
  }
  return dE_val;
}
*/
void matrix_elements(const Wavefunction &wf, const DiracSpinor &Fv,
                     const DiracSpinor &Fw, const std::string op,
                     const std::string type, const bool contact,
                     const bool tdhf, const double omega, const double min_mu,
                     const double max_mu, const int N_mu) {

  if (type == "va") {
    std::cout << "\n********WARNING: VA interaction does not yet include "
                 "spatial terms; "
                 "results likely incorrect.********\n";
  }

  const double y_sps = 1;

  // Check that limits agree

  const auto test_max_mu = 1000.0;
  const auto test_min_mu = 1e-6;

  // Test ground with all core states
  // std::cout << "\nTesting V_nv in limits for all basis states |n> with valence "
  //           << Fv.symbol()
  //           << ".\nShowing >10% discrepancies.\n\nFor the exact "
  //              "cases,\nμ->0 = "
  //           << test_min_mu << "\nμ->∞ = " << test_max_mu;

  // std::cout << "\n\nCheck OK:\n|v>  κv   |n>  κn     massless exact (μ->0) "
  //              " rel diff exact "
  //              "(μ->∞)      contact  rel diff\n";

  // for (auto Fn : wf.basis()) {
  //   if (Fn.twoj() == Fv.twoj()) {
  //     const auto V_massless = V_nv_direct(false, wf.core(), Fv, Fn, 1.0, 0.0);
  //     const auto V_exact_min =
  //         V_nv_direct(false, wf.core(), Fv, Fn, 1.0, test_min_mu);
  //     const auto V_exact_max =
  //         V_nv_direct(false, wf.core(), Fv, Fn, 1.0, test_max_mu);
  //     const auto V_contact =
  //         V_nv_direct(true, wf.core(), Fv, Fn, 1.0, test_max_mu);

  //     const auto diff_massless =
  //         std::abs((V_massless - V_exact_min) / V_massless);
  //     const auto diff_contact = std::abs((V_contact - V_exact_max) / V_contact);

  //     if ((diff_massless > 0.1) || (diff_contact > 0.1)) {
  //       fmt::print("{:3s} {:3}  {:4s} {:3} {:12.3e} {:12.3e} {:9.3f} "
  //                  "{:12.3e} {:12.3e} {:9.3f}\n",
  //                  Fv.shortSymbol(), Fv.kappa(), Fn.shortSymbol(), Fn.kappa(),
  //                  V_massless, V_exact_min, diff_massless, V_exact_max,
  //                  V_contact, diff_contact);
  //     }
  //   }
  // }

  std::vector<double> mus;

  // Populate mus
  for (double log_mu = log(min_mu); log_mu <= log(max_mu);
       log_mu += std::abs(log(max_mu) - log(min_mu)) * (1.0 / N_mu)) {
    mus.push_back(std::exp(log_mu));
  }

  const auto actual_N_mu = mus.size();
  std::vector<double> MEs_TDHF(actual_N_mu);
  std::vector<double> MEs(actual_N_mu);
  std::vector<double> MEs_contact(actual_N_mu);
  std::vector<double> Dvs_old(actual_N_mu);
  // std::vector<double> i0_maxs;
  // std::vector<double> k0_maxs;

  const auto Fn = wf.valence()[1];

  const bool debug = false;

  if (debug) {
    std::cout << "\n\nTESTING MATRIX ELEMENTS\n\n";

    double ME_manual_contact = 0.0;
    double ME_operator_contact = 0.0;

    for (auto Fn : wf.basis()) {

      for (auto Fa : wf.core()) {

        const auto Fdir = BSM_Vee::u_anav_contact(wf, Fn, Fa, Fv);
        const auto Fexch = BSM_Vee::u_anva_contact(wf, Fn, Fa, Fv);

        const auto Rdir_VA = BSM_Vee::R_AV_abcd(Fa, Fn, Fa, Fv);
        const auto Rdir_AV = BSM_Vee::R_AV_abcd(Fn, Fa, Fv, Fa);
        const auto Rexch_VA = BSM_Vee::R_AV_abcd(Fa, Fn, Fv, Fa);
        const auto Rexch_AV = BSM_Vee::R_AV_abcd(Fn, Fa, Fa, Fv);

        const auto Adir =
          (Rdir_AV + Rdir_VA) != 0 ? Fdir / (Rdir_AV + Rdir_VA) : 0.0;
        const auto Aexch =
          (Rexch_AV + Rexch_VA) != 0 ? Fexch / (Rexch_AV + Rexch_VA) : 0.0;

        if ((Fdir != 0.0) || (Fexch != 0.0)) {
          std::cout << "__________________________\n";
          std::cout << "Fn=" << Fn.shortSymbol() << "\tFv=" << Fv.shortSymbol()
                    << "\tFa=" << Fa.shortSymbol() << "\n";

          std::cout << "   Rdir_VA = " << Rdir_VA << "\n";
          std::cout << "   Rdir_AV = " << Rdir_AV << "\n";
          std::cout << "  Rexch_VA = " << Rexch_VA << "\n";
          std::cout << "  Rexch_AV = " << Rexch_AV << "\n\n";

          std::cout << "      Rdir = " << Rdir_AV + Rdir_VA << "\n";
          std::cout << "     Rexch = " << Rexch_AV + Rexch_VA << "\n";
          std::cout << "Rdir/Rexch = "
                    << (Rdir_AV + Rdir_VA) / (Rexch_AV + Rexch_VA) << "\n\n";

          std::cout << "   4π*Adir = " << 4 * M_PI * Adir << "\n";
          std::cout << "  4π*Aexch = " << 4 * M_PI * Aexch << "\n";
          std::cout << "Adir/Aexch = " << Adir / Aexch << "\n\n";

          std::cout << "      Vdir = " << Fdir << "\n";
          std::cout << "    2Vexch = " << 2 * Fexch << "\n";
          std::cout << "Vdir/Vexch = " << Fdir / Fexch << "  (-1?)\n\n";

          std::cout << "Vdir-Vexch = " << Fdir - Fexch << "\n";
          std::cout << "     2Vdir = " << 2 * Fdir << "\n";

          ME_manual_contact += Fdir - Fexch;
        }
      }

      ME_operator_contact += Fn * BSM_Vee::V_Fv(wf.core(), false, Fv, "va",
                                                Fn.kappa(), 1.0, true, 1.0);
    }
    std::cout << "\nFinal result for contact lim: " << ME_manual_contact;
    std::cout << "\n              Using operator: " << ME_operator_contact
              << "\n";

    std::cout << "\nAngular integrals with Fv = " << Fv.shortSymbol() << "\n\n";
    std::cout << " Red = diag * (-1)^k * [k] * <a||Ck||c> * <b||Ck||d>\n";
    std::cout << "Full =        (-1)^q * [k] * <a|Ck|c>  *  <b|Ck|d>\n";
    std::cout << "Dir VA: abcd -> nav(-a)\n";
    std::cout << "Dir AV: abcd -> na(-v)a\n";
    std::cout << "Exc VA: abcd -> naa(-v)\n";
    std::cout << "Exc AV: abcd -> na(-a)v\n\n";

    std::cout << "           Dir VA    Dir AV    Exc VA    Exc AV\n";
    std::cout << " Fn  Fa   Red Full  Red Full  Red Full  Red Full\n";

    // Check angular diagrams

    for (auto Fn : wf.basis()) {
      const int tmv = std::min(Fv.twoj(), Fn.twoj());
      for (auto Fa : wf.core()) {
        const auto A_direct_AV = std::sqrt(Fa.twojp1() * (1.0 / Fv.twojp1())) *
                                 Angular::Ck_kk(0, Fn.kappa(), -Fv.kappa()) *
                                 Angular::Ck_kk(0, Fa.kappa(), Fa.kappa());

        const auto A_direct_VA = std::sqrt(Fa.twojp1() * (1.0 / Fv.twojp1())) *
                                 Angular::Ck_kk(0, Fn.kappa(), Fv.kappa()) *
                                 Angular::Ck_kk(0, Fa.kappa(), -Fa.kappa());

        double A_exch_VA = 0;
        double A_exch_AV = 0;

        // Manual sum
        double A_direct_man_AV = 0;
        double A_direct_man_VA = 0;
        double A_exch_man_VA = 0;
        double A_exch_man_AV = 0;
        for (int k = 0; k <= (Fv.twoj() + Fa.twoj()); ++k) {

          const auto A_cc_VA = Angular::Ck_kk(k, Fn.kappa(), Fa.kappa()) *
                               Angular::Ck_kk(k, Fa.kappa(), -Fv.kappa());
          const auto A_cc_AV = Angular::Ck_kk(k, Fn.kappa(), -Fa.kappa()) *
                               Angular::Ck_kk(k, Fa.kappa(), Fv.kappa());

          const auto exch_factor = Angular::neg1pow_2(Fv.twoj() - Fa.twoj()) *
                                   (k * 2 + 1) * (1.0 / Fv.twojp1());

          const auto delta = Fv.twoj() == Fn.twoj() ? 1.0 : 0.0;

          A_exch_VA += A_cc_VA * exch_factor * delta;
          A_exch_AV += A_cc_AV * exch_factor * delta;

          for (int tq = -2 * k; tq <= 2 * k; tq += 2) {
            for (int tma = -Fa.twoj(); tma <= Fa.twoj(); tma += 2) {
              for (int tmn = -Fn.twoj(); tmn <= Fn.twoj(); tmn += 2) {
                const auto factor = Angular::neg1pow_2(tq) * (2 * k + 1);
                const auto dir_CC_AV =
                  Angular::Ck_kk_mmq(k, Fn.kappa(), -Fv.kappa(), tmn, tmv, tq) *
                  Angular::Ck_kk_mmq(k, Fa.kappa(), Fa.kappa(), tma, tma, -tq);

                const auto dir_CC_VA =
                  Angular::Ck_kk_mmq(k, Fn.kappa(), Fv.kappa(), tmn, tmv, tq) *
                  Angular::Ck_kk_mmq(k, Fa.kappa(), -Fa.kappa(), tma, tma, -tq);

                const auto exch_CC_AV =
                  Angular::Ck_kk_mmq(k, Fn.kappa(), -Fa.kappa(), tmn, tma, tq) *
                  Angular::Ck_kk_mmq(k, Fa.kappa(), Fv.kappa(), tma, tmv, -tq);

                const auto exch_CC_VA =
                  Angular::Ck_kk_mmq(k, Fn.kappa(), Fa.kappa(), tmn, tma, tq) *
                  Angular::Ck_kk_mmq(k, Fa.kappa(), -Fv.kappa(), tma, tmv, -tq);

                A_direct_man_AV += factor * dir_CC_AV;
                A_direct_man_VA += factor * dir_CC_VA;
                A_exch_man_VA += factor * exch_CC_VA;
                A_exch_man_AV += factor * exch_CC_AV;
              }
            }
          }
        }

        // Tabular
        if ((Fn < Fv) &&
            (abs(A_direct_AV) + abs(A_direct_man_AV) + abs(A_exch_AV) +
               abs(A_exch_man_VA) + abs(A_exch_AV) + abs(A_exch_man_AV) >
             1e-10)) {

          A_direct_man_VA = A_direct_man_VA > 1e-10 ? A_direct_man_VA : 0.0;
          std::cout << Fn.shortSymbol() << " " << Fa.shortSymbol() << "     "
                    << A_direct_VA << "    " << A_direct_man_VA << " |  "
                    << A_direct_AV << "    " << A_direct_man_AV << " |  "
                    << A_exch_VA << "    " << A_exch_man_VA << " |  "
                    << A_exch_AV << "    " << A_exch_man_AV << "\n";
        }

        // // Pretty for each state
        // if (A_exch_man_AV > 0.001) {
        //   std::cout << "_____________________________";
        //   std::cout << "\nAngular integrals with Fn=" << Fn.shortSymbol()
        //             << " Fv=" << Fv.shortSymbol() << " Fa=" << Fa.shortSymbol()
        //             << "\n";
        //   std::cout << "     Direct reduced: " << A_direct << "\n";
        //   std::cout << "        Direct full: " << A_direct_man << "\n\n";
        //   std::cout << "VA Exchange reduced: " << A_exch_VA << "\n";
        //   std::cout << "   VA Exchange full: " << A_exch_man_VA << "\n\n";
        //   std::cout << "AV Exchange reduced: " << A_exch_AV << "\n";
        //   std::cout << "   AV Exchange full: " << A_exch_man_AV << "\n";
        //   sleep(1);
        // }
      }
    }

    std::cout << "\n\nChecking individual multipoles k\n";
    std::cout << "               Dir VA    Dir AV    Exc VA    Exc AV\n";
    std::cout << " Fn  Fa  k    Red Full  Red Full  Red Full  Red Full\n";

    // Check angular diagrams

    for (auto Fn : wf.basis()) {
      const int tmv = std::min(Fv.twoj(), Fn.twoj());
      for (auto Fa : wf.core()) {
        auto A_direct_AV = std::sqrt(Fa.twojp1() * (1.0 / Fv.twojp1())) *
                           Angular::Ck_kk(0, Fn.kappa(), -Fv.kappa()) *
                           Angular::Ck_kk(0, Fa.kappa(), Fa.kappa());

        auto A_direct_VA = std::sqrt(Fa.twojp1() * (1.0 / Fv.twojp1())) *
                           Angular::Ck_kk(0, Fn.kappa(), Fv.kappa()) *
                           Angular::Ck_kk(0, Fa.kappa(), -Fa.kappa());

        for (int k = 0; k <= (Fv.twoj() + Fa.twoj()); ++k) {
          A_direct_AV = k == 0 ? A_direct_AV : 0;
          A_direct_VA = k == 0 ? A_direct_VA : 0;

          double A_exch_VA = 0;
          double A_exch_AV = 0;

          // Manual sum
          double A_direct_man_AV = 0;
          double A_direct_man_VA = 0;
          double A_exch_man_VA = 0;
          double A_exch_man_AV = 0;

          const auto A_cc_VA = Angular::Ck_kk(k, Fn.kappa(), Fa.kappa()) *
                               Angular::Ck_kk(k, Fa.kappa(), -Fv.kappa());
          const auto A_cc_AV = Angular::Ck_kk(k, Fn.kappa(), -Fa.kappa()) *
                               Angular::Ck_kk(k, Fa.kappa(), Fv.kappa());

          const auto exch_factor = Angular::neg1pow_2(Fv.twoj() - Fa.twoj()) *
                                   (k * 2 + 1) * (1.0 / Fv.twojp1());

          const auto delta = Fv.twoj() == Fn.twoj() ? 1.0 : 0.0;

          A_exch_VA += A_cc_VA * exch_factor * delta;
          A_exch_AV += A_cc_AV * exch_factor * delta;

          for (int tq = -2 * k; tq <= 2 * k; tq += 2) {
            for (int tma = -Fa.twoj(); tma <= Fa.twoj(); tma += 2) {
              for (int tmn = -Fn.twoj(); tmn <= Fn.twoj(); tmn += 2) {
                const auto factor = Angular::neg1pow_2(tq) * (2 * k + 1);
                const auto dir_CC_AV =
                  Angular::Ck_kk_mmq(k, Fn.kappa(), -Fv.kappa(), tmn, tmv, tq) *
                  Angular::Ck_kk_mmq(k, Fa.kappa(), Fa.kappa(), tma, tma, -tq);

                const auto dir_CC_VA =
                  Angular::Ck_kk_mmq(k, Fn.kappa(), Fv.kappa(), tmn, tmv, tq) *
                  Angular::Ck_kk_mmq(k, Fa.kappa(), -Fa.kappa(), tma, tma, -tq);

                const auto exch_CC_AV =
                  Angular::Ck_kk_mmq(k, Fn.kappa(), -Fa.kappa(), tmn, tma, tq) *
                  Angular::Ck_kk_mmq(k, Fa.kappa(), Fv.kappa(), tma, tmv, -tq);

                // // DIY exch_CC_AV
                // const auto exch_na_AV =
                //     Angular::neg1pow_2(Fn.twoj() - tmn) *
                //     Angular::threej_2(Fn.twoj(), 2 * k, Fa.twoj(), -tmn, tq,
                //                       tma) *
                //     Angular::Ck_kk(k, Fn.kappa(), -Fa.kappa());

                // const auto exch_av_AV =
                //     Angular::neg1pow_2(Fa.twoj() - tma) *
                //     Angular::threej_2(Fa.twoj(), 2 * k, Fv.twoj(), -tma, -tq,
                //                       tmv) *
                //     Angular::Ck_kk(k, Fa.kappa(), Fv.kappa());

                // const auto exch_CC_AV = exch_na_AV * exch_av_AV;

                const auto exch_CC_VA =
                  Angular::Ck_kk_mmq(k, Fn.kappa(), Fa.kappa(), tmn, tma, tq) *
                  Angular::Ck_kk_mmq(k, Fa.kappa(), -Fv.kappa(), tma, tmv, -tq);

                A_direct_man_AV += factor * dir_CC_AV;
                A_direct_man_VA += factor * dir_CC_VA;
                A_exch_man_VA += factor * exch_CC_VA;
                A_exch_man_AV += factor * exch_CC_AV;
              }
            }
          }

          // Tabular
          if ((Fn < Fv) &&
              (abs(A_direct_AV) + abs(A_direct_man_AV) + abs(A_exch_AV) +
                 abs(A_exch_man_VA) + abs(A_exch_AV) + abs(A_exch_man_AV) >
               1e-10)) {
            A_direct_man_VA = A_direct_man_VA > 1e-10 ? A_direct_man_VA : 0.0;
            std::cout << Fn.shortSymbol() << " " << Fa.shortSymbol() << "  "
                      << k << "     " << A_direct_VA << "    "
                      << A_direct_man_VA << " |  " << A_direct_AV << "    "
                      << A_direct_man_AV << " |  " << A_exch_VA << "    "
                      << A_exch_man_VA << " |  " << A_exch_AV << "    "
                      << A_exch_man_AV << "\n";

            if (k >= 0.5 * (Fv.twoj() + Fa.twoj())) {
              std::cout << "\n";
            }
          }
        }
      }
    }
  }

  if (tdhf) {
    std::cout << "\nRunning TDHF.\n";
  } else {
    std::cout << "\nTDHF = false, so ME_tdhf yields 'nan'.\n";
  }

  const double ME_massless =
    calc_ME(op, type, false, 0.0, false, omega, wf, Fw, Fv);

  std::vector<double> MEs_massless(actual_N_mu, ME_massless);

#pragma omp parallel default(none) shared(                                     \
    op, type, contact, tdhf, omega, wf, Fw, Fv, MEs_TDHF, MEs, actual_N_mu,    \
      ME_massless, MEs_contact, MEs_massless, Dvs_old, mus, std::cout)
  {

#pragma omp for
    for (int i = 0; i < actual_N_mu; ++i) {

      const auto mu = mus[i];
      const double ME = calc_ME(op, type, false, mu, false, omega, wf, Fw, Fv);

      //Approximations
      const double ME_contact =
        calc_ME(op, type, true, mu, false, omega, wf, Fw, Fv);

      // Old method for Dv
      double Dv_old = 0.0;
      if (op == "D") {

        // const auto i0_max = BSM_Vee::mod_sph_bessel_i(0.0, mu * wf.grid().rmax());
        // const auto k0_max = BSM_Vee::mod_sph_bessel_k(0.0, mu * wf.grid().rmax());

        // Old method
        for (auto Fn : wf.basis()) {
          if (Fn != Fv) {
            double d_wn_old = d_ab(wf.grid(), Fw, Fn);
            double V_nv_old = V_nv_direct(contact, wf.core(), Fv, Fn, 1.0, mu);

            if (Fv == Fw) {
              Dv_old += 2.0 * (d_wn_old * V_nv_old) / (Fv.en() - Fn.en());
            } else {
              Dv_old += (d_wn_old * V_nv_old) / (Fv.en() - Fn.en());
            }
          }

          if (Fv != Fw && Fn != Fw) {

            double V_wn_old = V_nv_direct(contact, wf.core(), Fn, Fw, 1.0, mu);
            double d_nv_old = d_ab(wf.grid(), Fn, Fv);

            Dv_old += (V_wn_old * d_nv_old) / (Fw.en() - Fn.en());
          }
        }

        Dv_old = Dv_old / PhysConst::alpha;
      }

      double ME_TDHF = 0;
      if (tdhf) {
        ME_TDHF = calc_ME(op, type, contact, mu, true, omega, wf, Fw, Fv);
      } else {
        ME_TDHF = NAN;
      }

// Store
#pragma omp critical
      {
        Dvs_old[i] = Dv_old;
        MEs[i] = ME;
        MEs_massless[i] = ME_massless;
        MEs_contact[i] = ME_contact;
        MEs_TDHF[i] = ME_TDHF;
      }
    }
  }

  if (op == "D") {
    std::cout
      << "\nCalculating <" << Fw.symbol() << "|D|" << Fv.symbol()
      << "> (au) with " << type
      << " interaction (mediator mass = μ).\n   μ (m_e) "
         "      D_old  D_massless   D_contact           D      D_tdhf\n";
    for (int i = 0; i < mus.size(); ++i) {
      fmt::print("{:10.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e}\n",
                 mus[i], Dvs_old[i], MEs_massless[i], MEs_contact[i], MEs[i],
                 MEs_TDHF[i]);
    }
  } else {
    std::cout << "\nCalculating <" << Fw.symbol() << "|" << op << "|"
              << Fv.symbol() << "> (au) with " << type
              << " interaction (mediator mass = μ).\n   μ (m_e)"
                 "  "
              << op << "_massless   " << op << "_contact           " << op
              << "      " << op << "_tdhf\n";
    for (int i = 0; i < mus.size(); ++i) {
      fmt::print("{:10.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e}\n", mus[i],
                 MEs_massless[i], MEs_contact[i], MEs[i], MEs_TDHF[i]);
    }
  }
}

double calc_ME(const std::string op, const std::string type, const bool contact,
               const double mu, const bool tdhf, const double omega,
               const Wavefunction &wf, const DiracSpinor &Fw,
               const DiracSpinor &Fv) {
  if (op == "D") {
    return calc_Dwv(type, contact, mu, tdhf, omega, wf, Fw, Fv);
  } else if (op == "V") {
    return calc_Vwv(type, contact, mu, tdhf, omega, wf, Fw, Fv);
  } else {
    return 0.0;
  }
}

double calc_Vwv(const std::string type, const bool contact, const double mu,
                const bool tdhf, const double omega, const Wavefunction &wf,
                const DiracSpinor &Fw, const DiracSpinor &Fv) {

  DiracOperator::Vee VeeOp(wf.core(), contact, mu, type, false);

  if (tdhf) {
    ExternalField::TDHF tdhf_Vee(&VeeOp, wf.vHF());
    tdhf_Vee.solve_core(omega, 100, false);

    // If epsilon bad, re-run and print
    if (tdhf_Vee.last_eps() > 1.0e8) {
      std::cout << "CHECK μ = " << mu << "\n";
      tdhf_Vee.solve_core(omega, 100, true);
    }

    return (VeeOp.fullME(Fw, Fv) +
            VeeOp.rme3js(Fw.twoj(), Fv.twoj()) * tdhf_Vee.dV(Fw, Fv)) /
           PhysConst::alpha;
  } else {
    return VeeOp.fullME(Fw, Fv) / PhysConst::alpha;
  }
}

// double calc_Dwv(const std::string type, const bool contact, const double mu,
//                 const bool tdhf, const Wavefunction &wf) {
//   const auto Fv = wf.valence()[0];
//   auto Fw(Fv);
//   if (type == "va") {
//     if (wf.valence().size() > 1) {
//       Fw = wf.valence()[1];
//     } else {
//       Fw = wf.core()[wf.core().size() - 1];
//     }
//   }
//   return calc_Dwv(type, contact, mu, tdhf, wf, Fv, Fw);
// }

// double calc_Dwv(const std::string type, const bool contact, const double mu,
//                 const bool tdhf, const Wavefunction &wf,
//                 const DiracSpinor &Fv) {
//   auto Fw(Fv);
//   if (type == "va") {
//     if (wf.valence().size() > 1) {
//       Fw = wf.valence()[1];
//     } else {
//       Fw = wf.core()[wf.core().size() - 1];
//     }
//   }
//   return calc_Dwv(type, contact, mu, tdhf, wf, Fv, Fw);
// }

double calc_Dwv(const std::string type, const bool contact, const double mu,
                const bool tdhf, const double omega, const Wavefunction &wf,
                const DiracSpinor &Fw, const DiracSpinor &Fv) {
  // std::cout << "In test mode";

  // For testing - normally false
  const bool eN = false;

  if (eN) {
    std::cout << "\nWARNING: Calculating 'eN' rather than 'ee' - results are "
                 "not robust.\n";
  } else if (type == "va" && contact == false) {
    std::cout << "\nWARNING: Only contact limit (via Fierz identity) supported "
                 "for VA at present. Defaulting to contact = true.\n";
  }

  double D_wv = 0.0;
  DiracOperator::Vee VeeOp(wf.core(), contact, mu, type, eN);
  DiracOperator::E1 E1(wf.grid());

  ExternalField::DiagramRPA drpa_Vee(&VeeOp, wf.basis(), wf.vHF(),
                                     wf.atomicSymbol());
  ExternalField::TDHF tdhf_d(&E1, wf.vHF());

  if (tdhf) {
    std::cout << "\nμ = " << mu << "\n";
    tdhf_d.solve_core(omega, 100, true);
    drpa_Vee.solve_core(0.0, 100, true);
  }

  for (auto Fn : wf.basis()) {
    if (Fn != Fv) {
      // For testing: y gets used as VFv = y * Fdir - (2 - y) * Fexch such that
      // y = 1.0: VFv = Fdir - Fexch
      // y = 2.0: VFv = 2.0 * Fdir

      /*
      std::cout << "\n"
                << Fn * BSM_Vee::V_Fv(wf.core(), Fv, "va", Fn.kappa(), 1.0,
                                      true, 1.0);
      std::cout << "\n"
                << Fn * BSM_Vee::V_Fv(wf.core(), Fv, "va", Fn.kappa(), 2.0,
                                      true, 1.0)
                << "\n=?\n";
                */

      // Find non-tdhf matrix elements
      double V_nv = VeeOp.fullME(Fn, Fv);
      double d_wn = E1.fullME(Fw, Fn);

      if (tdhf) {
        d_wn += E1.rme3js(Fw.twoj(), Fn.twoj()) * tdhf_d.dV(Fw, Fn);
        V_nv += VeeOp.rme3js(Fn.twoj(), Fv.twoj()) * drpa_Vee.dV(Fn, Fv);
      }

      // std::cout << "<" << Fw.kappa() << "|d|" << Fn.kappa() << "> = " << d_wn
      //           << "\t <" << Fn.kappa() << "|V|" << Fv.kappa() << "> = " << V_nv
      //           << "\n";

      if (Fv == Fw) {
        D_wv += 2.0 * (d_wn * V_nv) / (Fv.en() - Fn.en());
      } else {
        D_wv += (d_wn * V_nv) / (Fv.en() - Fn.en());
      }
    }

    if (Fv != Fw && Fn != Fw) {
      double V_wn = VeeOp.fullME(Fw, Fn);
      double d_nv = E1.fullME(Fn, Fv);

      if (tdhf) {
        V_wn += VeeOp.rme3js(Fw.twoj(), Fn.twoj()) * drpa_Vee.dV(Fw, Fn);
        d_nv += E1.rme3js(Fn.twoj(), Fv.twoj()) * tdhf_d.dV(Fn, Fv);
      }

      D_wv += (V_wn * d_nv) / (Fw.en() - Fn.en());

      // const auto d_nv = d_ab(wf.grid(), Fv, Fn);

      // <n|V|v>

      // const auto Vsps =
      // Fn * Vee::V_Fv(wf.core(), Fv, "sp", Fn.kappa(), 1.0, mu);

      // Using old implementation:
      // const auto V_old = V_nv_direct(false, wf.core(), Fv, Fn, 1, mu);

      // Using new implementation:
      // const auto Vsps = VeeOp.fullME(Fn, Fv);

      /*
      if (Vsps != V_old) {
        std::cout << "\n"
                  << Fn.symbol() << "\t" << V_old << "\t" << Vsps - V_old;
      }
      */

      // This numerically shows that <n|V|v> = <v|V|n>, i.e. V is Hermitian.
      // std::cout << "\n" << VeeOp.fullME(Fv, Fn) << "\t" << VeeOp.fullME(Fn, Fv);
    }
  }
  // Multiply by hbar c in au:
  return D_wv / PhysConst::alpha;
}

// ################ OLD FUNCTIONS ###################

double d_ab(const Grid &gr, const DiracSpinor &Fa, const DiracSpinor &Fb) {
  return DiracOperator::E1(gr).fullME(Fa, Fb);

  /* Manual - seems to match E1

  const auto direct_d_vn =
      -Angular::threej_2(Fv.twoj(), Fv.twoj(), 2, 1, -1, 0) *
      Angular::Ck_kk(1, Fv.kappa(), Fn.kappa()) *
      d.radialIntegral(Fv, Fn);

  */
}

double V_nv_direct(const bool contact, const std::vector<DiracSpinor> core,
                   const DiracSpinor &Fv, const DiracSpinor &Fn, const double y,
                   const double mu, const bool g0_both) {

  // For safety
  if (Fv.kappa() != -Fn.kappa()) {
    return 0.0;
  }

  auto u_anav = 0.0;
  auto u_naav = 0.0;
  auto u_anva = 0.0;

  for (auto Fa : core) {
    const auto R0_anav = contact == true ? R_abcd_contact(mu, Fa, Fn, Fa, Fv) :
                         mu == 0.0       ? Rk_abcd_massless(0, Fa, Fn, Fa, Fv) :
                                           Rk_abcd(0, mu, Fa, Fn, Fa, Fv);
    u_anav += R0_anav * Fa.twojp1();

    const auto phase =
      std::pow(-1.0, 0.5 * (Fv.twoj() - Fa.twoj())) * (1.0 / Fv.twojp1());

    for (int twok = std::abs(Fa.twoj() - Fv.twoj());
         twok <= Fa.twoj() + Fv.twoj(); twok += 2) {
      if ((Fa.twoj() + Fv.twoj() + twok) % 4 == 0) {
        const double k = 0.5 * twok;
        const auto A_anva = Angular::Ck_kk(k, Fa.kappa(), Fv.kappa()) *
                            Angular::Ck_kk(k, Fn.kappa(), -Fa.kappa());
        const auto R_anva = contact == true ?
                              R_abcd_contact(mu, Fa, Fn, Fv, Fa) :
                            mu == 0.0 ? Rk_abcd_massless(k, Fa, Fn, Fv, Fa) :
                                        Rk_abcd(k, mu, Fa, Fn, Fv, Fa);

        const auto A_naav = Angular::Ck_kk(k, Fn.kappa(), Fa.kappa()) *
                            Angular::Ck_kk(k, Fa.kappa(), -Fv.kappa());
        const auto R_naav = contact == true ?
                              R_abcd_contact(mu, Fn, Fa, Fa, Fv) :
                            mu == 0.0 ? Rk_abcd_massless(k, Fn, Fa, Fa, Fv) :
                                        Rk_abcd(k, mu, Fn, Fa, Fa, Fv);

        /* std::cout << "λ = " << k << "\tκa = " << Fa.kappa()
                  << "\tja = " << Fa.twoj() * 0.5
                  << "\tjv = " << Fv.twoj() * 0.5 << "\tFa = " << Fa.symbol()
                  << "\n";
        std::cout << "A_anva = " << A_anva << "\n"; */

        // std::cout << Fn * BSM_Vee::Bk_ab_v(k, contact, mu, false, "sp", Fa, Fv,
        //                                    Fa)
        //           << "\n";
        // std::cout << R_naav << "\n\n";

        u_anva += phase * (twok + 1) * A_anva * R_anva;
        u_naav += phase * (twok + 1) * A_naav * R_naav;
      }
    }
  }

  return (u_anav - u_naav - u_anva) * y;
}

double Rk_abcd(const double k, const double mu, const DiracSpinor &Fa,
               const DiracSpinor &Fb, const DiracSpinor &Fc,
               const DiracSpinor &Fd, const std::string int_type,
               const bool g0_both) {
  // Compare with yk_ab to find r> and r< functions
  // Then create radial operator?
  // Find the Rk_abcd implementation in code.

  //const auto i0 = std::max(Fa.min_pt(), Fc.min_pt());
  //const auto imax = std::min(Fa.max_pt(), Fc.max_pt());

  // const auto g_Fc = int_type == "vv" ? Fc : g0_both ? BSM_Vee::g0(Fc) : Fc;

  // const auto g_Fd = int_type == "sp" ? BSM_Vee::old_i_g0_g5(Fd) :
  //                   int_type == "ss" ? BSM_Vee::g0(Fd) :
  //                                      Fd;

  // Using new method funcs, only sp:
  const auto g_Fc = BSM_Vee::g0(Fc);
  const auto g_Fd = -1.0 * BSM_Vee::g5(Fd);

  const auto screening_function = Bk_ab(k, mu, Fb, g_Fd);

  const auto Rff =
    NumCalc::integrate(1.0, 0, Fa.grid().size(), Fa.f(), g_Fc.f(),
                       screening_function, Fa.grid().drdu());

  const auto Rgg =
    NumCalc::integrate(1.0, 0, Fa.grid().size(), Fa.g(), g_Fc.g(),
                       screening_function, Fa.grid().drdu());

  return (Rff + Rgg) * Fa.grid().du() * mu;
}

double R_abcd_contact(const double mu, const DiracSpinor &Fa,
                      const DiracSpinor &Fb, const DiracSpinor &Fc,
                      const DiracSpinor &Fd) {

  const auto &gr = Fa.grid();
  const auto &r = gr.r();
  const auto g0_Fc = BSM_Vee::g0(Fc);
  const auto ig0g5_Fd = BSM_Vee::old_i_g0_g5(Fd);

  // Delta case

  std::vector<double> integrand(gr.size());

  for (int i = 0; i < gr.size(); ++i) {
    integrand[i] = (Fa.f(i) * g0_Fc.f(i) + Fa.g(i) * g0_Fc.g(i)) *
                   (Fb.f(i) * ig0g5_Fd.f(i) + Fb.g(i) * ig0g5_Fd.g(i)) /
                   (gr.r(i) * gr.r(i));
  }

  return NumCalc::integrate(1.0, 0, gr.size(), integrand, gr.drdu()) * gr.du() /
         (mu * mu);

  /*
  // Bessel approx case 

  //const auto i0 = std::max(Fa.min_pt(), Fb.min_pt());
  //const auto imax = std::min(Fa.max_pt(), Fb.max_pt());

  // Integrate
  std::vector<double> result(gr.size());

  //Slower, more stable way
  std::vector<double> ik(gr.size());

  for (int i_mid = 0; i_mid < gr.size(); ++i_mid) {
    for (int i_gr = 0; i_gr < gr.size(); ++i_gr) {
      ik[i_gr] = std::exp(-mu * std::abs(r[i_mid] - r[i_gr])) /
                 (2 * mu * mu * r[i_gr] * r[i_mid]);
    }
    double B_ff =
        NumCalc::integrate(1.0, 0, gr.size(), ik, Fb.f(), Fd.f(), gr.drdu());

    double B_gg =
        NumCalc::integrate(1.0, 0, gr.size(), ik, Fb.g(), Fd.g(), gr.drdu());

    result[i_mid] = (B_ff + B_gg) * gr.du();
  }

  Faster, less stable way
  std::vector<double> i_k(gr.size());
  std::vector<double> k_k(gr.size());

  for (int i_gr = 0; i_gr < gr.size(); ++i_gr) {
    const auto x = mu * r[i_gr];
    //i_k[i_gr] = Vee::mod_sph_bessel_i(k, x);
    //k_k[i_gr] = Vee::mod_sph_bessel_k(k, x);

    i_k[i_gr] = std::exp(x) / (2 * x);
    k_k[i_gr] = std::exp(-x) / x;

    // For testing
    //i_k[i_gr] = 1.0;
    //k_k[i_gr] = 1.0;
  }

  for (int i_mid = 0; i_mid < gr.size(); ++i_mid) {
    double lower_ff =
        NumCalc::integrate(1.0, 0, i_mid, i_k, Fb.f(), Fd.f(), gr.drdu());

    double lower_gg =
        NumCalc::integrate(1.0, 0, i_mid, i_k, Fb.g(), Fd.g(), gr.drdu());

    // For r0 point
    if (i_mid == 0) {
      lower_ff = 0;
      lower_gg = 0;
    }

    const double upper_ff = NumCalc::integrate(1.0, i_mid, gr.size(), k_k,
                                               Fb.f(), Fd.f(), gr.drdu());

    const double upper_gg = NumCalc::integrate(1.0, i_mid, gr.size(), k_k,
                                               Fb.g(), Fd.g(), gr.drdu());

    result[i_mid] = (k_k[i_mid] * (lower_ff + lower_gg) +
                     i_k[i_mid] * (upper_ff + upper_gg)) *
                    gr.du();
  } 

  

  const auto B_bd = result;

  const auto Rff = NumCalc::integrate(1.0, 0, Fa.grid().size(), Fa.f(),
                                      ig0g5_Fc.f(), B_bd, Fa.grid().drdu());

  const auto Rgg = NumCalc::integrate(1.0, 0, Fa.grid().size(), Fa.g(),
                                      ig0g5_Fc.g(), B_bd, Fa.grid().drdu());

  return (Rff + Rgg) * Fa.grid().du() * mu;

  */
}

double Rk_abcd_massless(const double k, const DiracSpinor &Fa,
                        const DiracSpinor &Fb, const DiracSpinor &Fc,
                        const DiracSpinor &Fd) {
  const auto screening_function =
    Coulomb::yk_ab(k, Fb, BSM_Vee::old_i_g0_g5(Fd));

  const auto g0_Fc = BSM_Vee::g0(Fc);

  const auto Rff =
    NumCalc::integrate(1.0, 0, Fa.grid().size(), Fa.f(), g0_Fc.f(),
                       screening_function, Fa.grid().drdu());

  const auto Rgg =
    NumCalc::integrate(1.0, 0, Fa.grid().size(), Fa.g(), g0_Fc.g(),
                       screening_function, Fa.grid().drdu());

  return (Rff + Rgg) * Fa.grid().du() / (2 * k + 1);
}

std::vector<double> Bk_ab(const double k, const double mu,
                          const DiracSpinor &Fa, const DiracSpinor &Fb) {
  const auto &gr = Fa.grid();
  const auto &r = gr.r();
  //const auto i0 = std::max(Fa.min_pt(), Fb.min_pt());
  //const auto imax = std::min(Fa.max_pt(), Fb.max_pt());

  // Modified spherical Bessel functions
  std::vector<double> i_k(gr.size());
  std::vector<double> k_k(gr.size());

  for (int i_gr = 0; i_gr < gr.size(); ++i_gr) {
    const auto x = mu * r[i_gr];
    i_k[i_gr] = BSM_Vee::mod_sph_bessel_i(k, x);
    k_k[i_gr] = BSM_Vee::mod_sph_bessel_k(k, x);

    //i_k[i_gr] = std::exp(x) / (2 * x);
    //k_k[i_gr] = std::exp(-x) / x;

    // For testing
    //i_k[i_gr] = 1.0;
    //k_k[i_gr] = 1.0;
  }

  // Integrate
  std::vector<double> result(gr.size());

  for (int i_mid = 0; i_mid < gr.size(); ++i_mid) {
    double lower_ff =
      NumCalc::integrate(1.0, 0, i_mid, i_k, Fa.f(), Fb.f(), gr.drdu());

    double lower_gg =
      NumCalc::integrate(1.0, 0, i_mid, i_k, Fa.g(), Fb.g(), gr.drdu());

    // For r0 point
    if (i_mid == 0) {
      lower_ff = 0;
      lower_gg = 0;
    }

    const double upper_ff =
      NumCalc::integrate(1.0, i_mid, gr.size(), k_k, Fa.f(), Fb.f(), gr.drdu());

    const double upper_gg =
      NumCalc::integrate(1.0, i_mid, gr.size(), k_k, Fa.g(), Fb.g(), gr.drdu());

    result[i_mid] = (k_k[i_mid] * (lower_ff + lower_gg) +
                     i_k[i_mid] * (upper_ff + upper_gg)) *
                    gr.du();
  }

  return result;
}

void sps_testing(const Wavefunction &wf, const bool contact) {

  if (wf.valence().size() < 2) {
    std::cout << "\nERROR: Cannot run SPS testing with <2 valence states.\n";
    return;
  } else if (wf.core().size() < 2) {
    std::cout << "\nERROR: Cannot run SPS testing with <2 core states.\n";
    return;
  }

  // Pick Fv/a0 and Fv/a1 to be the lowest nκ n(-κ) pair of valence/core states

  DiracSpinor Fv0(wf.valence()[0]);
  DiracSpinor Fv1(wf.valence()[1]);
  DiracSpinor Fa0(wf.core()[0]);
  DiracSpinor Fa1(wf.core()[1]);

  if (Fv0.kappa() != -Fv1.kappa()) {
    bool flag = false;
    for (int i = 0; i < wf.valence().size(); ++i) {
      for (int j = 0; j < wf.valence().size(); ++j) {
        Fv0 = wf.valence()[i];
        Fv1 = wf.valence()[j];

        if (Fv0.kappa() == -Fv1.kappa()) {
          flag = true;
          break;
        }
      }
      if (flag == true) {
        break;
      }
    }
    if (flag == false) {
      std::cout << "\n ERROR: No valence states exist with nonzero V_nv. "
                   "To run "
                   "testing, ensure that two valence states exist such "
                   "that κ1 = "
                   "-κ2.\n";
      return;
    }
  }

  // Create spinor F1 satisfying (ff + gg) = 1 (i.e, f = g = 1/sqrt(2))
  DiracSpinor F1(Fv0);
  std::vector<double> ones(wf.grid().size(), 1.0 / std::sqrt(2.0));
  F1.f() = ones;
  F1.g() = ones;

  // Create spinor Fr satisfying (ff + gg) = r^2 (i.e., f = g = r/sqrt(2))
  DiracSpinor Fr(F1);

  for (int i = 0; i < wf.grid().size(); ++i) {
    Fr.f(i) = ones[i] * wf.grid().r(i);
    Fr.g(i) = ones[i] * wf.grid().r(i);
  }

  std::cout << "\nRunning checks on the SPS module with the following "
               "states\nname  state  κ   j\n"
            << " Fv0 " << Fv0.symbol(true) << " " << Fv0.kappa() << " "
            << Fv0.twoj() * 0.5 << "\n Fv1 " << Fv1.symbol(true) << " "
            << Fv1.kappa() << " " << Fv1.twoj() * 0.5 << "\n Fa0 "
            << Fa0.symbol() << " " << Fa0.kappa() << " " << Fa0.twoj() * 0.5
            << "\n Fa1 " << Fa1.symbol() << " " << Fa1.kappa() << " "
            << Fa1.twoj() * 0.5;

  std::cout << "\nAnd the special states 'F1' with f = g = 1/sqrt(2) such that "
               "F1*F1 = ff+gg "
               "= 1,\n\tand 'γF1' where f = -g = 1/sqrt(2) such that γF1*F1 = "
               "-fg+gf = 1,\n\tand 'Fr' where f = g = r/sqrt(2) such that "
               "Fr*Fr = ff+gg = r^2\n";

  // <v|d|v> - should be 0
  const auto vdv = d_ab(wf.grid(), Fv0, Fv0);

  // <v0|d|v1> and <v1|d|v0> - should be symmetric
  const auto v0dv1 =
    d_ab(wf.grid(), Fv0, Fv1); // should be 5.144 for Fr (<7s_1/2|d|7p_1/2>)
  const auto v1dv0 = d_ab(wf.grid(), Fv1, Fv0);

  // <v0||d||v1> - RME should be 5.144 for Fr (<7s_1/2||d||7p_1/2), https://arxiv.org/pdf/2212.11490 Table II
  const auto z_ab =
    v0dv1 / (std::pow(-1.0, (Fv0.twoj() - 1) * 0.5) *
             Angular::threej_2(Fv0.twoj(), 2, Fv1.twoj(), -1, 0, 1));

  std::cout << "\nDipole operator checks for v0 = " << Fv0.symbol()

            << " and v1 = " << Fv1.symbol() << "\n<v0|d|v0> = " << vdv
            << "\t\t=0?\n<v0|d|v1> = " << v0dv1
            << "\t=<v1|d|v0>?\n<v1|d|v0> = " << v1dv0
            << "\t=<v0|d|v1>?\n<v0||d||v1> = " << z_ab << "\tif Fr, =5.144?\n";

  // Check that the radial contact integral Rδ_abcd is returning expected orthonormality

  // R_aaaa = R_abab = R_abac = 0 due to γ5
  const auto R_aaaa = R_abcd_contact(1.0, Fv0, Fv0, Fv0, Fv0);
  const auto R_abab = R_abcd_contact(1.0, Fv0, Fv1, Fv0, Fv1);
  const auto R_abac = R_abcd_contact(1.0, Fv0, Fv1, Fv0, Fa0);
  const auto R_1a1b = R_abcd_contact(1.0, Fr, Fv0, Fr, Fv1);

  // Need another special function

  // R_(γ5*1)a1a = R_(γ5*a)1a1 = 1
  const auto R_g1a1a =
    R_abcd_contact(1.0, BSM_Vee::old_i_g0_g5(Fr), Fv0, Fr, Fv0);
  const auto R_ga1a1 =
    R_abcd_contact(1.0, BSM_Vee::old_i_g0_g5(Fv0), Fr, Fv0, Fr);

  std::cout << "\nRadial contact approximation checks\nR_aaaa = " << R_aaaa
            << "\t=0?\nR_abab = " << R_abab << "\t=0?\nR_abac = " << R_abac
            << "\t=0?\nR_1a1b = " << R_1a1b << "\t=0?\n"
            << "R_(γ1)a1a = " << R_g1a1a << "\t=1?\nR_(γa)1a1 = " << R_ga1a1
            << "\t=1?\n\n";

  // Check that λ = 0 Bessel function integral Bλ_ab produces expected results
  // Set (fafb + gagb) = 1 for all r
  // Then,
  // B0_11(r) = k0(r) int_0^r dr' i0(r') + i0(r) int_r^inf k0(r') dr'
  //          =       k0(r) * Shi[r]     -      i0(r) * Ei[-r]

  const auto B0_11 = Bk_ab(0.0, 1.0, F1, F1);

  // Check 10 values
  std::cout << "Bessel function integration (Bλ_ab) checks with λ=0 and "
               "(fafb + gagb) = 1 for all r\n         r    f    g numeric "
               "  exact \n";
  for (int i = 0; i < wf.grid().size();
       i += std::round(wf.grid().size() / 10)) {
    const auto ri = wf.grid().r(i);
    const auto i0 = BSM_Vee::mod_sph_bessel_i(0.0, ri);
    const auto k0 = BSM_Vee::mod_sph_bessel_k(0.0, ri);

    // B0_11(r) = k0(r)\int_0^r dr' i0(r) + i0(r)\int_r^\infty dr' k0(r)
    // B0_11(r) = k0(r)Shi[r] - i0(r)Ei[-r]
    const auto B0_11_exact =
      k0 * (gsl_sf_Shi(ri) - gsl_sf_Shi(0.000001)) -
      i0 * (gsl_sf_expint_Ei(-ri) - gsl_sf_expint_Ei(-150));

    fmt::print("{:10.6f} {:4.2f} {:4.2f} {:7.4f} {:7.4f}\n", ri, F1.f(i),
               F1.g(i), B0_11_exact, B0_11[i]);
  }

  // Check that radial integral is integrating Bk_ab as expected.
  // With (-fagc + fcga) = (fbfd + gbgd) = 1, should have
  // R0_(γ1)111 = \int_0^inf B0_bd

  const auto R0_g1111 = Rk_abcd(0.0, 1.0, BSM_Vee::old_i_g0_g5(F1), F1, F1, F1);
  const auto R0_g1111_manual =
    NumCalc::integrate(1.0, 0.0, wf.grid().size(), B0_11, wf.grid().drdu()) *
    wf.grid().du();

  std::cout << "\nRadial function integration with λ=0 and (-fagc + "
               "fcga)=(fbfd "
               "+ gbgd)=1 for all r\nR0_(γ1)1111\t = "
            << R0_g1111 << "\nExact \t = " << R0_g1111_manual
            << "\nAre these equal?\n";

  // Check radial integral's orthonormality is as expected. Same as contact case
  // Consider λ = 1 and μ = 1

  // R1_aaaa = R1_abab = R1_abac = 0 due to γ5
  const auto R1_aaaa = Rk_abcd(1.0, 1.0, Fv0, Fv0, Fv0, Fv0);
  const auto R1_abab = Rk_abcd(1.0, 1.0, Fv0, Fv1, Fv0, Fv1);
  const auto R1_abac = Rk_abcd(1.0, 1.0, Fv0, Fv1, Fv0, Fa0);
  const auto R1_1a1b = Rk_abcd(1.0, 1.0, F1, Fv0, F1, Fv1);

  std::cout << "\nRadial integration orthonormality checks with λ = μ = "
               "1\nR1_aaaa = "
            << R1_aaaa << "\t=0?\nR1_abab = " << R1_abab
            << "\t=0?\nR1_abac = " << R1_abac << "\t=0?\nR1_1a1b = " << R1_1a1b
            << "\t=0?\n";

  // V_nv_direct checks

  // Compare contact approximation with large μ

  const auto max_mu = 1000.0;
  const auto min_mu = 1e-6;

  // Test ground with all core states
  std::cout << "\nTesting V_nv in limits for all valence |v> and basis |n> "
               "states.\nShowing >10% discrepancies.\n\nFor the exact "
               "cases,\nμ->0 = "
            << min_mu << "\nμ->∞ = " << max_mu;

  std::cout << "\n\nCheck OK:\n|v>  κv   |n>  κn     massless exact (μ->0) "
               " rel diff exact "
               "(μ->∞)      contact  rel diff\n";

  for (auto Fv : wf.valence()) {
    for (auto Fn : wf.basis()) {
      if (Fn.twoj() == Fv.twoj()) {
        const auto V_massless = V_nv_direct(false, wf.core(), Fv, Fn, 1.0, 0.0);
        const auto V_exact_min =
          V_nv_direct(false, wf.core(), Fv, Fn, 1.0, min_mu);
        const auto V_exact_max =
          V_nv_direct(false, wf.core(), Fv, Fn, 1.0, max_mu);
        const auto V_contact =
          V_nv_direct(true, wf.core(), Fv, Fn, 1.0, max_mu);

        const auto diff_massless =
          std::abs((V_massless - V_exact_min) / V_massless);
        const auto diff_contact =
          std::abs((V_contact - V_exact_max) / V_contact);

        if ((diff_massless > 0.1) || (diff_contact > 0.1)) {
          fmt::print("{:3s} {:3}  {:4s} {:3} {:12.3e} {:12.3e} {:9.3f} "
                     "{:12.3e} {:12.3e} {:9.3f}\n",
                     Fv.shortSymbol(), Fv.kappa(), Fn.shortSymbol(), Fn.kappa(),
                     V_massless, V_exact_min, diff_massless, V_exact_max,
                     V_contact, diff_contact);
        }
      }
    }
    std::cout << "\n";
  }

  std::cout << "\nCalculating V_nv matrix elements for range of μ for <"
            << Fv1.symbol() << "|V|" << Fv0.symbol()
            << ">\n         mu"
               "    massless        full     contact  i0(150μ)  k0(150μ)\n";

  const auto V_v0v1_massless =
    V_nv_direct(false, wf.core(), Fv0, Fv1, 1.0, 0.0);

  for (double log_mu = log(min_mu); log_mu < log(max_mu);
       log_mu += std::abs(log(max_mu) - log(min_mu)) / 100) {
    const auto mu = exp(log_mu);
    const auto V_v0v1_contact = V_nv_direct(true, wf.core(), Fv0, Fv1, 1.0, mu);

    auto V_v0v1_full = V_nv_direct(false, wf.core(), Fv0, Fv1, 1.0, mu);

    const auto i0_max = BSM_Vee::mod_sph_bessel_i(0.0, mu * 150.0);
    const auto k0_max = BSM_Vee::mod_sph_bessel_k(0.0, mu * 150.0);
    /*if (i0_max == 0 & k0_max == 0) {
      V_v0v1_full = 0.0;
    } */

    fmt::print("{:11.6f}  {:10.3e}  {:10.3e}  {:10.3e}  {:8.1e}  {:8.1e}\n", mu,
               V_v0v1_massless, V_v0v1_full, V_v0v1_contact, i0_max, k0_max);
  }
}

// Dv checks

} // namespace Module