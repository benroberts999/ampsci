#include "Modules/VQE.hpp"
#include "Angular/Angular.hpp"
#include "CI/CI.hpp"
#include "Coulomb/Coulomb.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "IO/InputBlock.hpp"
#include "LinAlg/Matrix.hpp"
#include "MBPT/Sigma2.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include "fmt/ostream.hpp"
#include "qip/Vector.hpp"
#include <array>
#include <fstream>
#include <vector>

namespace Module {

using nkIndex = DiracSpinor::Index;

//==============================================================================
//==============================================================================
//==============================================================================
// This is the actual module that runs:
void VQE(const IO::InputBlock &input, const Wavefunction &wf) {

  // Check input options:
  input.check(
      {{"ci_basis",
        "Basis used for CI expansion; must be a sub-set of full ampsci basis "
        "[default: 10spdf]"},
       {"J", "List of total angular momentum J for CI solutions (comma "
             "separated). Must be integers (two-electron only). [0]"},
       {"J+", "As above, but for EVEN CSFs only (takes precedence over J)."},
       {"J-", "As above, but for ODD CSFs (takes precedence over J)."},
       {"num_solutions", "Number of CI solutions to find (for each J/pi) [5]"},
       {"sigma1", "Include one-body MBPT correlations? [false]"},
       {"sigma2", "Include two-body MBPT correlations? [false]"},
       {"cis2_basis",
        "The subset of ci_basis for which the two-body MBPT corrections are "
        "calculated. Must be a subset of ci_basis. If existing sk file has "
        "more integrals, they will be used. [default: Nspdf, where N is "
        "maximum n for core + 3]"},
       {"s1_basis",
        "Basis used for the one-body MBPT diagrams (Sigma^1). These are the "
        "most important, so in general the default (all basis states) should "
        "be used. Must be a subset of full ampsci basis. [default: full "
        "basis]"},
       {"s2_basis",
        "Basis used for internal lines of the two-body MBPT diagrams "
        "(Sigma^2). Must be a subset of s1_basis. [default: s1_basis]"},
       {"n_min_core", "Minimum n for core to be included in MBPT [1]"},
       {"max_k",
        "Maximum k (multipolarity) to include when calculating new "
        "Coulomb integrals. Higher k often contribute negligably. Note: if qk "
        "file already has higher-k terms, they will be included. Set negative "
        "(or very large) to include all k. [6]"},
       {"write_integrals",
        "Writes orbitals, CSFs, CI matrix, and 1 and 2 particle "
        "integrals to plain text file [true]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  // Decide if should write single-particle integrals to file:
  const auto write_integrals = input.get("write_integrals", true);

  //----------------------------------------------------------------------------
  // Single-particle basis:
  std::cout << "\nConstruct single-particle basis:\n";

  // Determine the sub-set of basis to use in CI:
  const auto basis_string = input.get("ci_basis", std::string{"10spdf"});

  // Select from wf.basis() [MBPT basis], those which match input 'basis_string'
  const std::vector<DiracSpinor> ci_sp_basis =
      CI::basis_subset(wf.basis(), basis_string, wf.coreConfiguration());

  // Print info re: basis to screen:
  std::cout << "\nUsing " << DiracSpinor::state_config(ci_sp_basis) << " = "
            << ci_sp_basis.size() << " relativistic single-particle orbitals\n";

  // Write orbital list
  std::map<nkm, int> orbital_map;
  if (write_integrals) {
    std::string fname = wf.atomicSymbol() + "_orbitals_" +
                        DiracSpinor::state_config(ci_sp_basis) + ".txt";
    std::ofstream of(fname);
    int index = 0;
    fmt::print(of, "# index Symbol n kappa m\n");
    for (const auto &Fn : ci_sp_basis) {
      for (int twom = -Fn.twoj(); twom <= Fn.twoj(); twom += 2) {

        fmt::print(of, "{} {} {} {} {}/2\n", index, Fn.shortSymbol(), Fn.n(),
                   Fn.kappa(), twom);
        orbital_map.insert({{Fn.n(), Fn.kappa(), twom}, index});
        ++index;
      }
    }
    std::cout << "( = " << index
              << " single-particle orbitals, including m projections)\n";
  }

  //----------------------------------------------------------------------------

  // check if including MBPT corrections
  const auto include_Sigma1 = input.get("sigma1", false);
  const auto include_Sigma2 = input.get("sigma2", false);
  const auto max_k_Coulomb = input.get("max_k", 6);
  const auto include_MBPT = include_Sigma1 || include_Sigma2;

  // s1 and s2 MBPT basis
  const auto s1_basis_string = input.get("s1_basis");
  const auto &s1_basis = s1_basis_string ?
                             CI::basis_subset(wf.basis(), *s1_basis_string) :
                             wf.basis();
  const auto s2_basis_string = input.get("s2_basis");
  const auto &s2_basis = s2_basis_string ?
                             CI::basis_subset(wf.basis(), *s2_basis_string) :
                             s1_basis;

  // Ensure s2_basis is subset of s1_basis
  assert(s2_basis.size() <= s1_basis.size() &&
         "s2_basis must be a subset of s1_basis");

  // Split basis' into core/excited (for MBPT evaluations)
  const auto n_min_core = input.get("n_min_core", 1);
  const auto [core_s1, excited_s1] =
      MBPT::split_basis(s1_basis, wf.FermiLevel(), n_min_core);
  const auto [core_s2, excited_s2] =
      MBPT::split_basis(s2_basis, wf.FermiLevel(), n_min_core);

  // S2 corrections are included only for this subset of the CI basis:
  const auto Ncore = DiracSpinor::max_n(wf.core()) + 3;
  const auto cis2_basis_string =
      input.get("cis2_basis", std::to_string(Ncore) + "spdf");
  const auto &cis2_basis = CI::basis_subset(ci_sp_basis, cis2_basis_string);

  //----------------------------------------------------------------------------

  if (include_MBPT) {
    std::cout << "\nIncluding MBPT: "
              << (include_Sigma1 && include_Sigma2 ? "Σ_1 + Σ_2" :
                  include_Sigma1                   ? "Σ_1" :
                                                     "Σ_2")
              << "\n";
    std::cout << "Including core excitations from n ≥ " << n_min_core << "\n";
    if (max_k_Coulomb >= 0 && max_k_Coulomb < 50) {
      std::cout << "Including k ≤ " << max_k_Coulomb
                << " in Coulomb integrals (unless already calculated)\n";
    }
    if (include_Sigma1) {
      // if (wf.Sigma()) {
      //   std::cout << "With existing Correlation Potential for Σ_1:\n";
      //   wf.Sigma()->print_info();
      // } else
      {
        std::cout << "With basis for Σ_1: "
                  << DiracSpinor::state_config(s1_basis) << "\n";
      }
    }
    if (include_Sigma2) {
      std::cout << "With basis for Σ_2: " << DiracSpinor::state_config(s2_basis)
                << "\n";

      std::cout << "Including Σ_2 correction to Coulomb integrals up to: "
                << DiracSpinor::state_config(cis2_basis) << "\n";
    }
    std::cout << "\n";
  }

  //----------------------------------------------------------------------------

  // Lookup table; stores all qk's
  Coulomb::QkTable qk;
  {
    std::cout << "Calculate two-body Coulomb integrals: Q^k_abcd\n";

    const auto qk_filename = input.get("qk_file", wf.atomicSymbol() + ".qk");

    // Try to read from disk (may already have calculated Qk)
    qk.read(qk_filename);
    const auto existing = qk.count();
    {

      // Try to limit number of Coulomb integrals we calculate
      // use whole basis (these are used inside Sigma_2)
      // If not including MBPT, only need to caculate smaller set of integrals

      // First, calculate the integrals between ci basis states:
      {
        std::cout << "For: " << DiracSpinor::state_config(ci_sp_basis) << "\n"
                  << std::flush;
        const auto yk = Coulomb::YkTable(ci_sp_basis);
        qk.fill(ci_sp_basis, yk, max_k_Coulomb, false);
      }

      // Selection function for which Qk's to calculate.
      // For Sigma, we only need those with 1 or 2 core electrons
      // i.e., no Q_vwxy, Q_vabc, or Q_abcd
      // Note: we *do* need Q_vwxy for the CI part (but with smaller basis)
      const auto select_Q_sigma =
          [eF = wf.FermiLevel()](int, const DiracSpinor &s,
                                 const DiracSpinor &t, const DiracSpinor &u,
                                 const DiracSpinor &v) {
            // Only calculate Coulomb integrals with 1 or 2 electrons in the core
            auto num = MBPT::number_below_Fermi(s, t, u, v, eF);
            return num == 1 || num == 2;
          };

      // Then, add those required for Sigma_1 (unless we have matrix!)
      if (include_Sigma1 /*&& !wf.Sigma()*/) {
        const auto temp_basis = qip::merge(core_s1, excited_s1);
        std::cout << "and: " << DiracSpinor::state_config(temp_basis) << "\n"
                  << std::flush;
        const auto yk = Coulomb::YkTable(temp_basis);
        qk.fill_if(temp_basis, yk, select_Q_sigma, max_k_Coulomb, false);
      }

      // Then, add those required for Sigma_2 (unless we already did Sigma_1)
      if (include_Sigma2 && !(include_Sigma1 /*&& !wf.Sigma()*/)) {
        const auto temp_basis = qip::merge(core_s2, excited_s2);
        std::cout << "and: " << DiracSpinor::state_config(temp_basis) << "\n"
                  << std::flush;
        const auto yk = Coulomb::YkTable(temp_basis);
        qk.fill_if(temp_basis, yk, select_Q_sigma, max_k_Coulomb, false);
      }

      // print summary
      qk.summary();

      // If we calculated new integrals, write to disk
      const auto total = qk.count();
      assert(total >= existing);
      const auto new_integrals = total - existing;
      std::cout << "Calculated " << new_integrals << " new Coulomb integrals\n";
      if (new_integrals > 0) {
        qk.write(qk_filename);
      }
    }
    std::cout << "\n" << std::flush;
  }

  //----------------------------------------------------------------------------

  // Create lookup table for one-particle matrix elements, h1
  const auto h1 =
      /* 
          wf.Sigma() ?
          CI::calculate_h1_table(ci_sp_basis, *wf.Sigma(), include_Sigma1) :
     */
      CI::calculate_h1_table(ci_sp_basis, core_s1, excited_s1, qk,
                             include_Sigma1);

  //----------------------------------------------------------------------------

  // print all single-body integrals to file:
  std::string one_file = wf.atomicSymbol() + "_h1_" +
                         DiracSpinor::state_config(ci_sp_basis) + ".txt";
  if (write_integrals) {
    std::cout << "Writing one-particle (h1) integrals to file: " << one_file
              << "\n";
    std::ofstream h1_file(one_file);
    h1_file << "# a  b  h1_ab   ## (h1_ab = h1_ba)\n";
    for (const auto &a : ci_sp_basis) {
      for (const auto &b : ci_sp_basis) {
        // h1 is scalar operator, so kappa's must be equal!
        if (b.kappa() != a.kappa())
          continue;
        for (int twom = -a.twoj(); twom <= a.twoj(); twom += 2) {
          // m_a = m_b (scalar operator)

          const auto index_a = orbital_map[nkm{a.n(), a.kappa(), twom}];
          const auto index_b = orbital_map[nkm{a.n(), a.kappa(), twom}];
          // symmetric; only store 'smallest' index set:
          if (index_b < index_a)
            continue;

          const auto value = h1.getv(a, b);
          //should never happen:
          if (value == 0.0)
            continue;
          fmt::print(h1_file, "{} {} {:.8e}\n", index_a, index_b, value);
        }
      }
    }
  }

  //----------------------------------------------------------------------------

  // Writes g integrals to text file
  // Modify this to include Sigma_2!
  if (write_integrals) {
    std::string g_file = wf.atomicSymbol() + "_h2_" +
                         DiracSpinor::state_config(ci_sp_basis) + ".txt";

    std::cout << "Writing two-particle (h2 = g_vwxy) integrals to file: "
              << g_file << "\n";
    write_CoulombIntegrals(g_file, ci_sp_basis, orbital_map, qk);
  }

  //----------------------------------------------------------------------------
  // Calculate MBPT corrections to two-body Coulomb integrals
  // Fix filename: account for n_min_core!
  const auto Sk_filename = wf.atomicSymbol() + "_" +
                           std::to_string(n_min_core) + "_" +
                           DiracSpinor::state_config(excited_s2) +
                           (max_k_Coulomb >= 0 && max_k_Coulomb < 50 ?
                                "_" + std::to_string(max_k_Coulomb) :
                                "") +
                           ".sk";

  Coulomb::LkTable Sk;
  if (include_Sigma2) {
    std::cout << "Calculate two-body MBPT integrals: Σ^k_abcd\n";

    std::cout << "For: " << DiracSpinor::state_config(cis2_basis) << ", using "
              << DiracSpinor::state_config(excited_s2) << "\n";

    Sk = CI::calculate_Sk(Sk_filename, cis2_basis, core_s2, excited_s2, qk,
                          max_k_Coulomb, false);
    std::cout << "\n" << std::flush;
  }

  if (include_Sigma2 && write_integrals) {
    std::string s_file = wf.atomicSymbol() + "_s2_" +
                         DiracSpinor::state_config(cis2_basis) + ".txt";

    std::cout << "Writing two-particle correlation corrections (s2 = S_vwxy) "
                 "integrals to file: "
              << s_file << "\n";
    write_CoulombIntegrals(s_file, cis2_basis, orbital_map, Sk);
  }

  //----------------------------------------------------------------------------
  const auto J_list = input.get("J", std::vector<int>{0});
  const auto J_even_list = input.get("J+", J_list);
  const auto J_odd_list = input.get("J-", J_list);
  const auto num_solutions = input.get("num_solutions", 5);

  //----------------------------------------------------------------------------
  if (write_integrals) {
    // Print CSFs:
    std::cout << "\n";
    for (const auto J : J_even_list) {
      CI::PsiJPi psi{2 * J, +1, ci_sp_basis};
      std::string csf_file = wf.atomicSymbol() + "_" + std::to_string(J) + "+";
      fmt::print("{} CSFs for J={}, even parity: {}\n", psi.CSFs().size(), J,
                 csf_file);
      write_CSFs(psi.CSFs(), 2 * J, orbital_map, csf_file);

      // Construct the CI matrix:
      const auto Hci = include_Sigma2 ? CI::construct_Hci(psi, h1, qk, &Sk) :
                                        CI::construct_Hci(psi, h1, qk);
      write_H(Hci, csf_file);
      std::cout << "\n";
    }
    for (const auto J : J_odd_list) {
      CI::PsiJPi psi{2 * J, -1, ci_sp_basis};

      std::string csf_file = wf.atomicSymbol() + "_" + std::to_string(J) + "-";
      fmt::print("{} CSFs for J={}, odd parity: {}\n", psi.CSFs().size(), J,
                 csf_file);
      write_CSFs(psi.CSFs(), 2 * J, orbital_map, csf_file);

      // Construct the CI matrix:
      const auto Hci = include_Sigma2 ? CI::construct_Hci(psi, h1, qk, &Sk) :
                                        CI::construct_Hci(psi, h1, qk);
      write_H(Hci, csf_file);
      std::cout << "\n";
    }
  }

  //----------------------------------------------------------------------------
  std::cout << "\nRun actual CI:\n";

  // even parity:
  for (const auto J : J_even_list) {
    CI::run_CI(ci_sp_basis, int(std::round(2 * J)), +1, num_solutions, h1, qk,
               Sk, include_Sigma2);
  }

  // odd parity:
  for (const auto J : J_odd_list) {
    CI::run_CI(ci_sp_basis, int(std::round(2 * J)), -1, num_solutions, h1, qk,
               Sk, include_Sigma2);
  }
}

//==============================================================================
//==============================================================================

//==============================================================================
void write_H(const LinAlg::Matrix<double> &Hci, const std::string &csf_fname) {
  std::string ci_fname = csf_fname + "_H.txt";
  std::cout << "Writing CI matrix to file: " << ci_fname << "\n";
  std::ofstream ci_file(ci_fname);
  ci_file << "# in matrix/table form: \n";
  for (std::size_t iA = 0; iA < Hci.rows(); ++iA) {
    for (std::size_t iX = 0; iX < Hci.cols(); ++iX) {
      fmt::print(ci_file, "{:+.6e} ", Hci(iA, iX));
    }
    ci_file << "\n";
  }
}

//==============================================================================
void write_CSFs(const std::vector<CI::CSF2> &CSFs, int twoJ,
                const std::map<nkm, int> &orbital_map,
                const std::string &csf_fname) {

  std::cout << "Writing CSFs and projections to files: {csf/proj}-" << csf_fname
            << "\n";
  std::ofstream csf_file(csf_fname + "_csf.txt");
  std::ofstream proj_file(csf_fname + "_proj.txt");

  csf_file << "# csf_index  a  b\n";
  proj_file << "# proj_index  a  b  CGC\n";
  int csf_count = 0;
  int proj_count = 0;
  for (const auto &csf : CSFs) {

    const auto v = csf.state(0);
    const auto w = csf.state(1);

    const auto [nv, kv] = Angular::index_to_nk(v);
    const auto [nw, kw] = Angular::index_to_nk(w);

    const auto tjv = Angular::twoj_k(kv);
    const auto tjw = Angular::twoj_k(kw);

    fmt::print(csf_file, "{} {} {}\n", csf_count,
               DiracSpinor::shortSymbol(nv, kv),
               DiracSpinor::shortSymbol(nw, kw));
    ++csf_count;

    // Each individual m projection:
    for (int two_m_v = -tjv; two_m_v <= tjv; two_m_v += 2) {
      const auto two_m_w = twoJ - two_m_v;
      if (std::abs(two_m_w) > tjw)
        continue;
      const auto cgc = Angular::cg_2(tjv, two_m_v, tjw, two_m_w, twoJ, twoJ);
      if (cgc == 0.0)
        continue;

      const auto iv = orbital_map.at(nkm{nv, kv, two_m_v});
      const auto iw = orbital_map.at(nkm{nw, kw, two_m_w});
      if (iw > iv)
        continue;

      const auto eta = v == w ? 1 / std::sqrt(2.0) : 1.0;
      const auto d_proj = eta * cgc; //?

      fmt::print(proj_file, "{} {} {} {:.8f}\n", proj_count, iv, iw, d_proj);
      ++proj_count;
    }
  }
  std::cout << "Writing " << csf_count
            << " CSFs to file: " << (csf_fname + "_csf.txt") << "\n";
  std::cout << "Writing " << proj_count
            << " projections to file: " << (csf_fname + "_csf.txt") << "\n";
}

//==============================================================================

} // namespace Module
