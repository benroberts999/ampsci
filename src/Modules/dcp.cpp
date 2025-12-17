#include "Modules/dcp.hpp"
#include "Coulomb/include.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/include.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"

//==============================================================================
double C1(const DiracSpinor &w, const DiracSpinor &v, double omega, int K,
          int kt, int ks, const Coulomb::meTable<double> &T,
          const Coulomb::meTable<double> &S,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited) {

  double Ak = 0.0;
#pragma omp parallel for reduction(+ : Ak) collapse(2)
  for (const auto &m : excited) {
    for (const auto &n : excited) {

      const auto s_wm = Angular::neg1pow_2(w.twoj() + m.twoj());
      const auto Wmwvn = Coulomb::Wk_abcd(K, w, m, v, n);
      if (Wmwvn == 0.0)
        continue;

      for (const auto &a : core) {

        // C1_ts:
        const auto sj_ts = Angular::SixJ(kt, ks, K, m, n, a);
        if (sj_ts != 0.0) {
          const auto Tna = T.getv(n, a);
          const auto Sam = S.getv(a, m);
          const auto de_ts = (a.en() - n.en() + omega) * (a.en() - m.en());
          const auto s_ts = s_wm * Angular::neg1pow(kt + ks);
          Ak += Tna * Sam * Wmwvn / de_ts * s_ts * sj_ts;
        }

        // C1_st:
        const auto sj_st = ks == kt ? sj_ts : Angular::SixJ(ks, kt, K, m, n, a);
        if (sj_st != 0.0) {
          const auto Sna = S.getv(n, a);
          const auto Tam = T.getv(a, m);
          const auto de_st = (a.en() - n.en()) * (a.en() - m.en() - omega);
          const auto s_st = s_wm * Angular::neg1pow(K);
          Ak += Sna * Tam * Wmwvn / de_st * s_st * sj_st;
        }
      }
    }
  }
  return Ak / std::sqrt(2 * K + 1);
}

//==============================================================================
double C2(const DiracSpinor &w, const DiracSpinor &v, double omega, int K,
          int kt, int ks, const Coulomb::meTable<double> &T,
          const Coulomb::meTable<double> &S,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited) {

  const auto sK = Angular::neg1pow(K + 1);
  const auto sts = Angular::neg1pow(kt + ks + 1);

  double Ak = 0.0;
#pragma omp parallel for reduction(+ : Ak) collapse(2)
  for (const auto &a : core) {
    for (const auto &b : core) {

      const auto s_wa = Angular::neg1pow_2(w.twoj() + a.twoj());
      const auto Wwavb = Coulomb::Wk_abcd(K, w, a, v, b);
      if (Wwavb == 0.0)
        continue;

      for (const auto &n : excited) {

        // C2_ts
        const auto sj_ts = Angular::SixJ(kt, ks, K, b, a, n);
        if (sj_ts != 0.0) {
          const auto Tna = T.getv(n, a);
          const auto Sbn = S.getv(b, n);
          const auto de_ts = (a.en() - n.en() + omega) * (b.en() - n.en());
          Ak += s_wa * sK * sj_ts * Tna * Sbn * Wwavb / de_ts;
        }

        // C2_st
        const auto sj_st = ks == kt ? sj_ts : Angular::SixJ(ks, kt, K, b, a, n);
        if (sj_st != 0.0) {
          const auto Sna = S.getv(n, a);
          const auto Tbn = T.getv(b, n);
          const auto de_st = (a.en() - n.en()) * (b.en() - n.en() - omega);
          Ak += s_wa * sts * sj_st * Sna * Tbn * Wwavb / de_st;
        }
      }
    }
  }
  return Ak / std::sqrt(2 * K + 1);
}

//==============================================================================
double R12(const DiracSpinor &w, const DiracSpinor &v, double omega, int K,
           int kt, int ks, const Coulomb::meTable<double> &T,
           const Coulomb::meTable<double> &S,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited) {

  // Both R1 and R2 share a Wk integral; faster to combine

  const auto sK = Angular::neg1pow(K);
  const auto sts = Angular::neg1pow(kt + ks);

  double Ak = 0.0;
#pragma omp parallel for reduction(+ : Ak) collapse(2)
  for (const auto &a : core) {
    for (const auto &n : excited) {

      const auto Wwavn = Coulomb::Wk_abcd(K, w, a, v, n);
      if (Wwavn == 0.0)
        continue;

      // R1:
      const auto s_wa = Angular::neg1pow_2(w.twoj() + a.twoj());
      for (const auto &m : excited) {

        //R1_ts
        const auto sj_ts = Angular::SixJ(kt, ks, K, n, a, m);
        if (sj_ts != 0.0) {
          const auto Tma = T.getv(m, a);
          const auto Snm = S.getv(n, m);
          const auto de_ts =
              (a.en() - m.en() + omega) * (a.en() - n.en() + omega);
          Ak += s_wa * sK * sj_ts * Tma * Snm * Wwavn / de_ts;
        }

        //R1_st
        const auto sj_st = ks == kt ? sj_ts : Angular::SixJ(ks, kt, K, n, a, m);
        if (sj_st != 0.0) {
          const auto Sma = S.getv(m, a);
          const auto Tnm = T.getv(n, m);
          const auto de_st = (a.en() - m.en()) * (a.en() - n.en() + omega);
          Ak += s_wa * sts * sj_st * Sma * Tnm * Wwavn / de_st;
        }
      }

      // R2 (note -ve sign cf R1):
      for (const auto &b : core) {

        //R2_ts
        const auto sj_ts = Angular::SixJ(kt, ks, K, a, n, b);
        if (sj_ts != 0.0) {
          const auto Tnb = T.getv(n, b);
          const auto Sba = S.getv(b, a);
          const auto de_ts =
              (b.en() - n.en() + omega) * (a.en() - n.en() + omega);
          Ak -= s_wa * sts * sj_ts * Tnb * Sba * Wwavn / de_ts;
        }

        //R2_st
        const auto sj_st = ks == kt ? sj_ts : Angular::SixJ(ks, kt, K, a, n, b);
        if (sj_st != 0.0) {
          const auto Snb = S.getv(n, b);
          const auto Tba = T.getv(b, a);
          const auto de_st = (b.en() - n.en()) * (a.en() - n.en() + omega);
          Ak -= s_wa * sK * sj_st * Snb * Tba * Wwavn / de_st;
        }
      }
    }
  }
  return Ak / std::sqrt(2 * K + 1);
}

//==============================================================================
double L12(const DiracSpinor &w, const DiracSpinor &v, double omega, int K,
           int kt, int ks, const Coulomb::meTable<double> &T,
           const Coulomb::meTable<double> &S,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited) {

  // Both L1 and L2 share a Wk integral; faster to combine

  const auto sK = Angular::neg1pow(K);
  const auto sts = Angular::neg1pow(kt + ks);

  double Ak = 0.0;
#pragma omp parallel for reduction(+ : Ak) collapse(2)
  for (const auto &a : core) {
    for (const auto &n : excited) {

      const auto s_wn = Angular::neg1pow_2(w.twoj() + n.twoj());
      const auto Wwnva = Coulomb::Wk_abcd(K, w, n, v, a);
      if (Wwnva == 0.0)
        continue;

      // L1:
      for (const auto &m : excited) {

        //L1_ts
        const auto sj_ts = Angular::SixJ(kt, ks, K, a, n, m);
        if (sj_ts != 0.0) {
          const auto Tmn = T.getv(m, n);
          const auto Sam = S.getv(a, m);
          const auto de_ts = (a.en() - n.en() - omega) * (a.en() - m.en());
          Ak += s_wn * sK * sj_ts * Tmn * Sam * Wwnva / de_ts;
        }

        //L1_st
        const auto sj_st = ks == kt ? sj_ts : Angular::SixJ(ks, kt, K, a, n, m);
        if (sj_st != 0.0) {
          const auto Smn = S.getv(m, n);
          const auto Tam = T.getv(a, m);
          const auto de_st =
              (a.en() - n.en() - omega) * (a.en() - m.en() - omega);
          Ak += s_wn * sts * sj_st * Smn * Tam * Wwnva / de_st;
        }
      }

      // L2 (note -ve sign cf R1):
      for (const auto &b : core) {

        //L2_ts
        const auto sj_ts = Angular::SixJ(kt, ks, K, n, a, b);
        if (sj_ts != 0.0) {
          const auto Tab = T.getv(a, b);
          const auto Sbn = S.getv(b, n);
          const auto de_ts = (a.en() - n.en() - omega) * (b.en() - n.en());
          Ak -= s_wn * sts * sj_ts * Tab * Sbn * Wwnva / de_ts;
        }

        //L2_st
        const auto sj_st = ks == kt ? sj_ts : Angular::SixJ(ks, kt, K, n, a, b);
        if (sj_st != 0.0) {
          const auto Sab = S.getv(a, b);
          const auto Tbn = T.getv(b, n);
          const auto de_st =
              (a.en() - n.en() - omega) * (b.en() - n.en() - omega);
          Ak -= s_wn * sK * sj_st * Sab * Tbn * Wwnva / de_st;
        }
      }
    }
  }
  return Ak / std::sqrt(2 * K + 1);
}

//==============================================================================
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

  // Transition frequency (omega)
  const auto omega_default = w.en() - v.en();
  const auto omega = input.get("omega", omega_default);

  std::cout << "E_w - E_v = " << omega_default << "\n";
  std::cout << "omega     = " << omega << "\n";

  // Get required operators:
  // t is time-dependent, s is static
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
    dVt.solve_core(omega);
    dVs.solve_core(0.0);
  }

  int tm = std::min(w.twoj(), v.twoj());

  const auto C_K = [](int KK, int k1, int k2, int tjw, int tjn, int tjv) {
    return std::sqrt(2.0 * KK + 1.0) * Angular::neg1pow_2(tjw + tjv) *
           Angular::sixj_2(2 * KK, 2 * k2, 2 * k1, tjn, tjw, tjv);
  };

  double AK = 0.0;
  double A_zz = 0.0;
  for (const auto &n : wf.spectrum()) {
    const auto Twn = t->reducedME(w, n) + dVt.dV(w, n);
    const auto Snv = s->reducedME(n, v) + dVs.dV(n, v);
    const auto de_ts = v.en() - n.en();

    const auto C_ts =
        C_K(K, kt, ks, w.twoj(), n.twoj(), v.twoj()) * Angular::neg1pow(K);

    const auto Swn = s->reducedME(w, n) + dVs.dV(w, n);
    const auto Tnv = t->reducedME(n, v) + dVt.dV(n, v);
    const auto de_st = w.en() - n.en();

    const auto C_st = C_K(K, ks, kt, w.twoj(), n.twoj(), v.twoj()) *
                      Angular::neg1pow(kt + ks);

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
    std::cout << "alpha: " << -AK / std::sqrt(3.0 * w.twojp1()) << "\n";
  }

  if (K == 1 && s->name() == "E1") {
    const auto ss = 2.0 * std::sqrt(2.0) * Angular::S_kk(w.kappa(), v.kappa());
    std::cout << "beta: " << AK / ss << "\n";
  }

  if (K == 1 && s->name() == "pnc-nsi") {
    std::cout << "Epnc: " << AK * z_comp << "\n";
  }

  const auto [core, excited] =
      DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel());

  std::cout << "\nCalculating DCP:\n";
  std::cout << "Filling tables.." << std::flush;
  const auto t_ce =
      ExternalField::me_table(wf.basis(), t.get(), rpaQ ? &dVt : nullptr);
  const auto s_ce =
      ExternalField::me_table(wf.basis(), s.get(), rpaQ ? &dVs : nullptr);
  std::cout << "..done\n" << std::flush;

  double cc1 =
      C1(w, v, omega, K, t->rank(), s->rank(), t_ce, s_ce, core, excited);

  double cc2 =
      C2(w, v, omega, K, t->rank(), s->rank(), t_ce, s_ce, core, excited);

  double cr12 =
      R12(w, v, omega, K, t->rank(), s->rank(), t_ce, s_ce, core, excited);

  double cl12 =
      L12(w, v, omega, K, t->rank(), s->rank(), t_ce, s_ce, core, excited);

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
    std::cout << "alpha: " << -AK / std::sqrt(3.0 * w.twojp1()) << " + "
              << -delta_AK / std::sqrt(3.0 * w.twojp1()) << "\n";
  }

  if (K == 1 && s->name() == "E1") {
    const auto ss = 2.0 * std::sqrt(2.0) * Angular::S_kk(w.kappa(), v.kappa());
    std::cout << "beta: " << AK / ss << " + " << delta_AK / ss << "\n";
  }

  if (K == 1 && s->name() == "pnc-nsi") {
    std::cout << "Epnc: " << AK * z_comp << " + " << delta_AK * z_comp << "\n";
  }
}