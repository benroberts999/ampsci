#include "Amplitudes/DCP.hpp"
#include "Amplitudes/MatrixElements.hpp"
#include "Coulomb/include.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/include.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/Modules.hpp"
#include "Wavefunction/Wavefunction.hpp"

//==============================================================================
namespace Module {

// Declare, register, then define below.
void dcp(const IO::InputBlock &input, const Wavefunction &wf);
namespace {
const Register r_dcp{"dcp", "Module for testing dcp calculations", &dcp};
} // namespace

void dcp(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check({{"initial", ""},
               {"final", ""},
               {"rpa", ""},
               {"K", ""},
               {"s", ""},
               {"t", ""},
               {"omega", ""}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  using namespace std::string_literals;

  // Parse which valence states for transition
  const auto [ni, ki] = AtomData::parse_symbol(input.get("initial", ""s));
  const auto [nf, kf] = AtomData::parse_symbol(input.get("final", ""s));
  // Lookup states in wavefunction, check if exist
  const auto i_ptr = wf.getState(ni, ki);
  const auto f_ptr = wf.getState(nf, kf);
  if (!i_ptr || !f_ptr) {
    std::cerr << "\nFAIL 22\n"
              << "Couldn't find requested state. Is it in valence list?\n";
    return;
  }
  // Transition v->w => <w|...|v>
  const auto v = *i_ptr;
  const auto w = *f_ptr;

  std::cout << "\n";
  std::cout << v << " -> " << w << " : <" << w << "|D|" << v << ">\n";

  // Frequencies of t and s: energy conservation, e_w = e_v + omega_t + omega_s
  const auto omega_default = w.en() - v.en();
  const auto omega_t = input.get("omega", omega_default);
  const auto omega_s = omega_default - omega_t;

  std::cout << "E_w - E_v = " << omega_default << "\n";
  std::cout << "omega_t   = " << omega_t << "\n";
  std::cout << "omega_s   = " << omega_s << "\n";

  // Get required operators:
  const auto t_oper = input.get("t", ""s);
  const auto s_oper = input.get("s", ""s);
  const auto t = DiracOperator::generate(t_oper, {}, wf);
  const auto s = DiracOperator::generate(s_oper, {}, wf);

  const auto ks = s->rank();
  const auto kt = t->rank();

  const int K = input.get("K", -1);

  const auto [kmin, kmax] = std::array{std::abs(ks - kt), std::abs(ks + kt)};
  if (K < kmin || K > kmax) {
    std::cout << "\nFail 299.\n"
                 "K = "
              << K << " outside possible range: [" << kmin << ", " << kmax
              << "]\n";
    return;
  }

  // Solve RPA; if required
  auto dVt = ExternalField::TDHF(t.get(), wf.vHF());
  auto dVs = ExternalField::TDHF(s.get(), wf.vHF());
  const auto rpaQ = input.get("rpa", true);
  if (rpaQ) {
    dVt.solve_core(omega_t);
    dVs.solve_core(omega_s);
  }

  int tm = std::min(w.twoj(), v.twoj());

  const auto C_K = [](int KK, int k1, int k2, int tjw, int tjn, int tjv) {
    return std::sqrt(2.0 * KK + 1.0) * Angular::neg1pow_2(tjw + tjv) *
           Angular::sixj_2(2 * KK, 2 * k2, 2 * k1, tjn, tjw, tjv);
  };

  double AK = 0.0;
  double A_zz = 0.0;
  for (const auto &n : wf.spectrum()) {
    // Each dV is the component at the frequency of its operator
    const auto Twn = t->reducedME(w, n) + dVt.dV(w, n, omega_t < 0.0);
    const auto Snv = s->reducedME(n, v) + dVs.dV(n, v, omega_s < 0.0);
    const auto de_ts = v.en() + omega_s - n.en();

    const auto C_ts =
      C_K(K, kt, ks, w.twoj(), n.twoj(), v.twoj()) * Angular::neg1pow(K);

    const auto Swn = s->reducedME(w, n) + dVs.dV(w, n, omega_s < 0.0);
    const auto Tnv = t->reducedME(n, v) + dVt.dV(n, v, omega_t < 0.0);
    const auto de_st = v.en() + omega_t - n.en();

    const auto C_st =
      C_K(K, ks, kt, w.twoj(), n.twoj(), v.twoj()) * Angular::neg1pow(kt + ks);

    const auto AK_n = C_ts * Twn * Snv / de_ts + C_st * Swn * Tnv / de_st;
    AK += AK_n;

    const auto cz_ts = t->rme3js(w, n, tm) * s->rme3js(n, v, tm);
    const auto cz_st = s->rme3js(w, n, tm) * t->rme3js(n, v, tm);

    const auto Az_n = cz_ts * Twn * Snv / de_ts + cz_st * Swn * Tnv / de_st;
    A_zz += Az_n;
  }

  const auto z_comp = Angular::neg1pow_2(w.twoj() - tm) *
                      Angular::threej_2(w.twoj(), 2 * K, v.twoj(), -tm, 0, tm);

  std::cout << "\n";
  std::cout << "A^K  : " << AK << "\n";
  std::cout << "A^K_0: " << AK * z_comp << "\n";
  std::cout << "A_zz : " << A_zz << "\n";

  // Specific output:
  std::cout << "\n";
  if (K == 0 && s->name() == "E1") {
    std::cout << "alpha: " << AK / std::sqrt(3.0 * w.twojp1()) << "\n";
  }

  if (K == 1 && s->name() == "E1") {
    const auto ss = 2.0 * std::sqrt(2.0) * Angular::S_kk(w.kappa(), v.kappa());
    std::cout << "beta: " << AK / ss << "\n";
  }

  if (K == 1 && s->name() == "pnc") {
    std::cout << "Epnc: " << AK * z_comp << "\n";
  }

  const auto [core, excited] =
    DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel());

  std::cout << "\nCalculating DCP:\n";
  std::cout << "Filling tables.." << std::flush;
  // Each table carries the dV component at the frequency of its operator
  const auto t_ce = Amplitudes::me_table(
    wf.basis(), t.get(), rpaQ ? &dVt : nullptr, nullptr, omega_t);
  const auto s_ce = Amplitudes::me_table(
    wf.basis(), s.get(), rpaQ ? &dVs : nullptr, nullptr, omega_s);
  std::cout << "..done\n" << std::flush;

  double cc1 = Amplitudes::dcp_C1(w, v, omega_t, omega_s, K, t->rank(),
                                  s->rank(), t_ce, s_ce, core, excited);

  double cc2 = Amplitudes::dcp_C2(w, v, omega_t, omega_s, K, t->rank(),
                                  s->rank(), t_ce, s_ce, core, excited);

  double cr12 = Amplitudes::dcp_R12(w, v, omega_t, omega_s, K, t->rank(),
                                    s->rank(), t_ce, s_ce, core, excited);

  double cl12 = Amplitudes::dcp_L12(w, v, omega_t, omega_s, K, t->rank(),
                                    s->rank(), t_ce, s_ce, core, excited);

  const auto delta_AK = cc1 + cc2 + cr12 + cl12;

  std::cout << "\ndelta(AK):\n";
  fmt::print(" C1    = {:.4e}\n C2    = {:.4e}\n "
             "R     = {:.4e}\n L     = {:.4e}\n "
             "Total = {:.4e}\n\n",
             cc1, cc2, cr12, cl12, delta_AK);

  std::cout << "A^K  : " << AK << " + " << delta_AK << "\n";

  // Specific output:
  std::cout << "\n";
  if (K == 0 && s->name() == "E1") {
    std::cout << "alpha: " << AK / std::sqrt(3.0 * w.twojp1()) << " + "
              << delta_AK / std::sqrt(3.0 * w.twojp1()) << "\n";
  }

  if (K == 1 && s->name() == "E1") {
    const auto ss = 2.0 * std::sqrt(2.0) * Angular::S_kk(w.kappa(), v.kappa());
    std::cout << "beta: " << AK / ss << " + " << delta_AK / ss << "\n";
  }

  if (K == 1 && s->name() == "pnc") {
    std::cout << "Epnc: " << AK * z_comp << " + " << delta_AK * z_comp << "\n";
  }
}

} // namespace Module
