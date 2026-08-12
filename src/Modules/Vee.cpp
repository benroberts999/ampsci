#include "DiracOperator/Operators/Vee.hpp"
#include "Angular/Wigner369j.hpp"
#include "CI/CI_Integrals.hpp"
#include "DiracOperator/GenerateOperator.hpp"
#include "DiracOperator/Operators/Ek.hpp"
#include "DiracOperator/Operators/RadialF.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "IO/InputBlock.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Modules/Modules.hpp"
#include "Modules/Vee.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Potentials/Vee_Potentials.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "ampsci/ampsci.hpp"
#include "fmt/format.hpp"
#include <cmath>
#include <gsl/gsl_sf.h>
#include <iostream>
#include <omp.h>

namespace Module {

namespace {
const Register r_Vee{
  "Vee",
  "Vee: Calculates matrix elements of two-body electron-electron interactions.",
  &Vee};
} // namespace

void Vee(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check(
    {{"",
      "Calculates matrix elements of two-body electron-electron interactions."},
     {"type", "'sp' (scalar-pseudoscalar), 'va' (vector-axial "
              "vector), 'ss' (scalar-scalar), "
              "'vv' (vector-vector) ['sp']"},
     {"tdhf", "Include TDHF calcs for 'sp'? [false]"},
     {"omega", "ω (for tdhf) [0]"},
     {"approx",
      "Use 'contact' (μ->infty) or 'massless' (μ->0) approximation? [None]"},
     {"min_mu", "Minimum mediator mass [1e-6]"},
     {"max_mu", "Maximum mediator mass [20]"},
     {"N_mu", "Number of masses [100] (in massless limit, [1])"},
     {"v", "Initial state (ket, |v>) [ground]. Uses short-symbol, e.g. 6s+"},
     {"w", "Final state (bra, <w|) [v] (if type=va "
           "AND >1 valence states then [first-excited]). Uses short-symbol, "
           "e.g. 6p-"},
     {"operator",
      "Operator for matrix elements with Vee mixing [E1]."}, // TODO ?
     {"ci", "Use CI method with configuration-state wavefunctions? Requires CI "
            "block [false]"},
     {"n", "Principal quantum number for state [ground]"},
     {"kappa", "Kappa for state [ground]"},
     {"A2", "Second isotope's mass (for 'ss' or 'vv') [A+5]"},
     {"g0", "Include the gamma-0 term on both electrons? [true]"},
     {"ci_basis", "ci_basis, to be removed!"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }
  const auto e_handler = gsl_set_error_handler_off();

  // TODO: Use enum
  const std::string int_type = input.get<std::string>("type", "sp");
  const std::string approx = input.get<std::string>("approx", "None");
  const bool contact = approx == "contact";
  const bool massless = approx == "massless";

  if (!contact && !massless && (approx != "None")) {
    std::cout << "\n\nWARNING: unidentified approximation '" << approx
              << "', using 'None'";
  }

  const double min_mu = massless ? 1e-10 : input.get<double>("min_mu", 1.0e-4);
  const int N_mu = massless ? 1 : input.get<int>("N_mu", 100);
  const double max_mu = input.get<double>("max_mu", 1.0e4);
  const bool tdhf = input.get<bool>("tdhf", false);
  const double omega = input.get<double>("omega", 0.0);

  std::string v_sym = input.get<std::string>("v", "");
  std::string w_sym = input.get<std::string>("w", "");

  if (v_sym == "") {
    if (wf.valence().size() == 0) {
      std::cout << "\n\nERROR: No valence states available!";
      return;
    }
    v_sym = wf.valence()[0].shortSymbol();
  }

  const DiracSpinor Fv = *wf.getState(v_sym);

  if (w_sym == "") {
    if (wf.valence().size() == 0) {
      std::cout << "\n\nERROR: No valence states available!";
      return;
    } else if (int_type != "va") {
      w_sym = v_sym;
    } else if (wf.valence().size() == 1) {
      w_sym = v_sym;
      std::cout << "\n\nWARNING: 'va' interaction chosen but only one valence "
                   "state exists.";
    } else {
      w_sym = wf.valence()[1].shortSymbol();
    }
  }

  const DiracSpinor Fw = *wf.getState(w_sym);
  const std::string op = input.get<std::string>("operator", "D");

  if (op != "D") {
    std::cout << "\n\nERROR: Only mixing with E1 currently supported.";
    return;
  }

  if ((int_type == "va") && (!contact)) {
    std::cout << "\n\nWARNING: Only contact limit currently supported for 'va' "
                 "interaction.";
  }

  if ((int_type != "sp") && (int_type != "va")) {
    std::cout << "\n\nERROR: Unsupported interaction" << int_type;
    std::cout << "\nOnly 'va' and 'sp' interactions currently supported.";
    return;
  }

  const bool ci = input.get<bool>("ci", false);

  if (ci) {
    test_CI(input, wf);
    return;
  }

  std::cout << "\nCalculating matrix elements <" << Fw.shortSymbol() << "|"
            << op << "|" << Fv.shortSymbol() << "> "
            << "perturbed by Vee\n";

  if (contact) {
    std::cout << "mu = \u221E (contact limit)";
  } else if (massless) {
    std::cout << "mu = 0 (massless limit)";
  } else {
    std::cout << "mu = " << min_mu;
    if (N_mu > 1) {
      std::cout << "..." << max_mu << " (N = " << N_mu << ")";
    }
  }

  std::cout << "\n\n";

  std::vector<double> mus(N_mu);
  std::vector<double> MEs_TDHF(N_mu);
  std::vector<double> MEs(N_mu);
  std::vector<double> MEs_contact(N_mu);

  if (tdhf) {
    std::cout << "\nRunning TDHF for each mu.\n";
  }

  const double ME_massless =
    contact ? NAN : calc_ME(op, int_type, false, 0.0, false, omega, wf, Fw, Fv);

  std::vector<double> MEs_massless(N_mu, ME_massless);

#pragma omp parallel default(none)                                             \
  shared(op, int_type, contact, massless, tdhf, N_mu, min_mu, max_mu, omega,   \
           wf, Fw, Fv, MEs_TDHF, MEs, ME_massless, MEs_contact, MEs_massless,  \
           mus, std::cout)
  {

#pragma omp for
    for (int i = 0; i < N_mu; ++i) {

      const auto mu = N_mu > 1 ?
                        std::exp(log(min_mu) + i * (log(max_mu) - log(min_mu)) *
                                                 (1.0 / (N_mu - 1))) :
                        min_mu;

      // w/o TDHF
      const double ME = massless || contact ? NAN :
                                              calc_ME(op, int_type, false, mu,
                                                      false, omega, wf, Fw, Fv);

      // w/o TDHF
      const double ME_contact =
        massless ? NAN :
                   calc_ME(op, int_type, true, mu, false, omega, wf, Fw, Fv);

      const double ME_TDHF =
        tdhf ? calc_ME(op, int_type, contact, mu, true, omega, wf, Fw, Fv) :
               NAN;
// Store
#pragma omp critical
      {
        mus[i] = mu;
        MEs[i] = ME;
        MEs_massless[i] = ME_massless;
        MEs_contact[i] = ME_contact;
        MEs_TDHF[i] = ME_TDHF;
      }
    }
  }

  std::cout << "\nCalculating <" << Fw.symbol() << "|" << op << "|"
            << Fv.symbol() << "> (au) with " << int_type
            << " interaction (mediator mass = μ).\n";

  std::cout << "   μ (m_e)";

  if (!contact) {
    std::cout << "    Massless";
  }
  if (!massless) {
    std::cout << "     Contact";
  }
  if (!contact && !massless) {
    std::cout << "     Full ME";
  }
  if (tdhf) {
    std::cout << "       +TDHF";
  }
  std::cout << "\n";

  for (int i = 0; i < mus.size(); ++i) {
    fmt::print("{:10.4e} ", mus[i]);

    if (!contact) {
      fmt::print("{:11.4e} ", MEs_massless[i]);
    }

    if (!massless) {
      fmt::print("{:11.4e} ", MEs_contact[i]);
    }

    if (!contact && !massless) {
      fmt::print("{:11.4e} ", MEs[i]);
    }

    if (tdhf) {
      fmt::print("{:11.4e} ", MEs_TDHF[i]);
    }
    std::cout << "\n";
  }

  // Cleanup
  gsl_set_error_handler(e_handler);
}

double calc_ME(const std::string op, const std::string type, const bool contact,
               const double mu, const bool tdhf, const double omega,
               const Wavefunction &wf, const DiracSpinor &Fw,
               const DiracSpinor &Fv) {
  if (op == "V") {
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
  } else {
    // std::cout << "In test mode";

    // For testing - normally false
    const bool eN = false;

    if (eN) {
      std::cout << "\nWARNING: Calculating 'eN' rather than 'ee' - results are "
                   "not robust.\n";
    } else if (type == "va" && contact == false) {
      std::cout
        << "\nWARNING: Only contact limit (via Fierz identity) supported "
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

    // Sum over states
    for (auto Fn : wf.basis()) {
      if (Fn != Fv) {

        // Find non-tdhf matrix elements
        double V_nv = VeeOp.fullME(Fn, Fv);
        double d_wn = E1.fullME(Fw, Fn);

        if (tdhf) {
          d_wn += E1.rme3js(Fw.twoj(), Fn.twoj()) * tdhf_d.dV(Fw, Fn);
          V_nv += VeeOp.rme3js(Fn.twoj(), Fv.twoj()) * drpa_Vee.dV(Fn, Fv);
        }

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
      }
    }
    // Multiply by hbar c in au:
    return D_wv / PhysConst::alpha;
  }
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

} // namespace Module