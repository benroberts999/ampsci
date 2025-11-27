#include "Modules/dcp.hpp"
#include "Coulomb/include.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/include.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"

//==============================================================================
double C1(const DiracSpinor &w, const DiracSpinor &v, double omega, int K,
          int kt, int ks, const Coulomb::meTable<double> &t_me,
          const Coulomb::meTable<double> &s_me,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited) {

  //

  // const auto ks = s->rank();
  // const auto kt = t->rank();

  std::vector<std::pair<std::size_t, std::size_t>> mn_index;
  for (std::size_t im = 0; im < excited.size(); im++) {
    for (std::size_t in = 0; in < excited.size(); in++) {
      const auto &m = excited[im];
      const auto &n = excited[in];

      const auto sj_ts_OK =
          Angular::sixjTriads(2 * kt, 2 * ks, 2 * K, m.twoj(), n.twoj(), {});
      const auto sj_st_OK =
          Angular::sixjTriads(2 * ks, 2 * kt, 2 * K, m.twoj(), n.twoj(), {});

      const auto Wk_ok = Coulomb::Qk_abcd_SR(K, w, m, v, n) ||
                         Coulomb::Pk_abcd_SR(K, w, m, v, n);

      if ((sj_ts_OK || sj_st_OK) && Wk_ok)
        mn_index.push_back({im, in});
    }
  }

  double Ak = 0.0;
#pragma omp parallel for reduction(+ : Ak)
  for (const auto &[im, in] : mn_index) {
    const auto &m = excited[im];
    const auto &n = excited[in];

    const auto ss = Angular::neg1pow_2(w.twoj() + m.twoj());

    const auto Wmwvn = Coulomb::Wk_abcd(K, w, m, v, n);

    for (const auto &a : core) {

      double AK_ts_anm{0.0}, AK_st_anm{0.0};
      const auto sj_ts =
          Angular::sixj_2(2 * kt, 2 * ks, 2 * K, m.twoj(), n.twoj(), a.twoj());

      // C1_ts:
      if (sj_ts != 0.0) {
        const auto Tna = t_me.getv(n, a);
        const auto Sam = s_me.getv(a, m);
        const auto de_ts =
            (a.en() - n.en() + omega) * (a.en() - m.en() - omega);
        const auto s_ts = ss * Angular::neg1pow(kt + ks);
        AK_ts_anm = Tna * Sam * Wmwvn / de_ts * s_ts * sj_ts;
      }

      const auto sj_st = kt == ks ?
                             sj_ts :
                             Angular::sixj_2(2 * ks, 2 * kt, 2 * K, m.twoj(),
                                             n.twoj(), a.twoj());

      // C1_st:
      if (sj_st != 0.0) {
        const auto Sna = s_me.getv(n, a);
        const auto Tam = t_me.getv(a, m);
        const auto de_st =
            (a.en() - n.en() + omega) * (a.en() - m.en() - omega);
        const auto s_st = ss * Angular::neg1pow(K);
        AK_st_anm = Sna * Tam * Wmwvn / de_st * s_st * sj_st;
      }

      Ak += AK_ts_anm + AK_st_anm;
    }
  }
  return Ak / std::sqrt(2 * K + 1);
}

//==============================================================================
double C2(const DiracSpinor &w, const DiracSpinor &v, double omega, int K,
          int kt, int ks, const Coulomb::meTable<double> &t_me,
          const Coulomb::meTable<double> &s_me,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited) {

  std::vector<std::pair<std::size_t, std::size_t>> ab_index;
  for (std::size_t ia = 0; ia < core.size(); ia++) {
    for (std::size_t ib = 0; ib < core.size(); ib++) {
      const auto &b = excited[ib];
      const auto &a = excited[ia];

      const auto sj_ts_OK =
          Angular::sixjTriads(2 * kt, 2 * ks, 2 * K, b.twoj(), a.twoj(), {});
      const auto sj_st_OK =
          Angular::sixjTriads(2 * ks, 2 * kt, 2 * K, b.twoj(), a.twoj(), {});

      const auto Wk_ok = Coulomb::Qk_abcd_SR(K, w, a, v, b) ||
                         Coulomb::Pk_abcd_SR(K, w, a, v, b);

      if ((sj_ts_OK || sj_st_OK) && Wk_ok)
        ab_index.push_back({ia, ib});
    }
  }

  double Ak = 0.0;
#pragma omp parallel for reduction(+ : Ak)
  for (const auto &[ia, ib] : ab_index) {
    const auto &a = excited[ia];
    const auto &b = excited[ib];

    // const auto ss = Angular::neg1pow_2(w.twoj() + m.twoj());

    // const auto Wmwvn = Coulomb::Wk_abcd(K, w, m, v, n);

    // for (const auto &a : core) {

    //   double AK_ts_anm{0.0}, AK_st_anm{0.0};
    //   const auto sj_ts =
    //       Angular::sixj_2(2 * kt, 2 * ks, 2 * K, m.twoj(), n.twoj(), a.twoj());

    //   // C1_ts:
    //   if (sj_ts != 0.0) {
    //     const auto Tna = t_me.getv(n, a);
    //     const auto Sam = s_me.getv(a, m);
    //     const auto de_ts =
    //         (a.en() - n.en() + omega) * (a.en() - m.en() - omega);
    //     const auto s_ts = ss * Angular::neg1pow(kt + ks);
    //     AK_ts_anm = Tna * Sam * Wmwvn / de_ts * s_ts * sj_ts;
    //   }

    //   const auto sj_st = kt == ks ?
    //                          sj_ts :
    //                          Angular::sixj_2(2 * ks, 2 * kt, 2 * K, m.twoj(),
    //                                          n.twoj(), a.twoj());

    //   // C1_st:
    //   if (sj_st != 0.0) {
    //     const auto Sna = s_me.getv(n, a);
    //     const auto Tam = t_me.getv(a, m);
    //     const auto de_st =
    //         (a.en() - n.en() + omega) * (a.en() - m.en() - omega);
    //     const auto s_st = ss * Angular::neg1pow(K);
    //     AK_st_anm = Sna * Tam * Wmwvn / de_st * s_st * sj_st;
    //   }

    //   Ak += AK_ts_anm + AK_st_anm;
    // }
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

  // Transition frequency (omega)
  const auto omega_default = w.en() - v.en();
  const auto omega = input.get("omega", omega_default);

  // Get required operators:
  // t is time-dependent, s is static
  const auto t_oper = input.get("t", ""s);
  const auto s_oper = input.get("s", ""s);
  const auto t = DiracOperator::generate(t_oper, {}, wf);
  const auto s = DiracOperator::generate(s_oper, {}, wf);

  const auto ks = s->rank();
  const auto kt = t->rank();

  const int K = input.get("K", -1);

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

    const auto cz_ts =
        t->rme3js(w.twoj(), n.twoj(), tm) * s->rme3js(n.twoj(), v.twoj(), tm);
    const auto cz_st =
        s->rme3js(w.twoj(), n.twoj(), tm) * t->rme3js(n.twoj(), v.twoj(), tm);

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

  if (K == 1 && s->name() == "pnc") {
    std::cout << "Epnc: " << AK * z_comp << "\n";
  }

  const auto [core, excited] =
      DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel());

  const auto t_ce = ExternalField::me_table(core, excited, t.get(), &dVt);
  const auto s_ce = ExternalField::me_table(core, excited, s.get(), &dVs);

  double cc1 =
      C1(w, v, omega, K, t->rank(), s->rank(), t_ce, s_ce, core, excited);
  std::cout << cc1 << "\n";
}