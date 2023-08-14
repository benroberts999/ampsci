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
      {{"frozen_core", "Core states that are not included in CI expansion. By "
                       "default, this is the same as the HartreeFock{} core"},
       {"basis", "Basis used for CI expansion; must be a sub-set of MBPT basis "
                 "[default: 5spd]"},
       {"J", "List of J angular symmetry for CSFs (comma separated). For "
             "half-integer, enter as floats: '0.5' not '1/2' [default: 0]"},
       {"J+", "As above, but for EVEN CSFs only (takes precedence over J)."},
       {"J-", "As above, but for ODD CSFs (takes precedence over J)."},
       {"num_solutions", "Number of CI solutions to find (for each J/pi) [5]"},
       {"sigma2", "Include two-body correlations? [false]"},
       {"e0", "Optional: ground-state energy (in 1/cm) for relative energies. "
              "If not given, will assume lowest J+"},
       {"write_integrals",
        "Writes orbitals, CSFs, CI matrix, and 1 and 2 particle "
        "integrals to plain text file [false]"},
       {"ci_input", "Input text file contiaining list of CI expansion "
                    "coeficients for each CSF (each on "
                    "new line); uses these to calculate CI energy"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  // Decide if should write single-particle integrals to file:
  const auto write_integrals = input.get("write_integrals", false);

  //----------------------------------------------------------------------------
  // Single-particle basis:
  std::cout << "\nConstruct single-particle basis:\n";

  // Determine the sub-set of basis to use in CI:
  const auto basis_string = input.get("basis", std::string{"5spd"});

  const auto frozen_core_string =
      input.get("frozen_core", wf.coreConfiguration());

  // Select from wf.basis() [MBPT basis], those which match input 'basis_string'
  const std::vector<DiracSpinor> ci_sp_basis =
      CI::basis_subset(wf.basis(), basis_string, frozen_core_string);

  // Print info re: basis to screen:
  std::cout << "\nUsing " << DiracSpinor::state_config(ci_sp_basis) << " = "
            << ci_sp_basis.size() << " relativistic single-particle orbitals\n";

  // Write orbital list
  std::map<nkm, int> orbital_map;
  if (write_integrals) {
    std::string fname =
        "orbitals_" + DiracSpinor::state_config(ci_sp_basis) + ".txt";
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

  std::cout << "\nSingle-particle matrix elements:\n";

  // Correlation potential (one-particle Sigma matrix):
  const auto Sigma = wf.Sigma();

  // Create lookup table for one-particle matrix elements, h1
  Coulomb::meTable<double> h1;

  // Calculate + store all 1-body integrals
  for (const auto &v : ci_sp_basis) {
    for (const auto &w : ci_sp_basis) {
      if (w > v)
        continue;
      if (w.kappa() != v.kappa())
        continue;
      const auto h0_vw = v == w ? v.en() : 0.0;
      // const auto h0_vw = wf.Hab(v, w);
      // nb: This makes a difference: energy denominators!
      // Take lowest n (from valence) of that state?
      const auto Sigma_vw = Sigma ? v * Sigma->SigmaFv(w) : 0.0;
      // if (w == v && Sigma_vw != 0.0) {
      //   std::cout << v << " " << w << " " << Sigma_vw << "\n";
      // }
      // const auto Sigma_vw = Sigma ? Sigma->Sigma_vw(v, w) : 0.0;
      h1.add(v, w, h0_vw + Sigma_vw);
      if (v != w)
        h1.add(w, v, h0_vw + Sigma_vw);
    }
  }

  // print all single-body integrals to file:
  std::string one_file =
      "h1_" + DiracSpinor::state_config(ci_sp_basis) + ".txt";
  if (write_integrals) {
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

  std::cout << "\nCalculate two-body Coulomb integrals: Q^k_abcd\n";

  // Lookup table; stores all qk's
  Coulomb::QkTable qk;
  const auto qk_filename =
      wf.identity() + DiracSpinor::state_config(wf.basis()) + ".qk";
  // Try to read from disk (may already have calculated Qk)
  const auto read_from_file_ok = qk.read(qk_filename);
  if (!read_from_file_ok) {
    // if didn't find Qk file to read in, calculate from scratch:

    // use whole basis (these are used iside Sigma_2)
    const Coulomb::YkTable yk(wf.basis());
    qk.fill(wf.basis(), yk);
    qk.write(qk_filename);
  }
  std::cout << std::flush;

  // Writes Rk integrals to text file
  // Modify this to include Sigma_2!
  if (write_integrals) {
    write_CoulombIntegrals(ci_sp_basis, orbital_map, qk);
  }

  //----------------------------------------------------------------------------
  const auto J_list = input.get("J", std::vector<double>{});
  const auto J_even_list = input.get("J+", J_list);
  const auto J_odd_list = input.get("J-", J_list);
  const auto num_solutions = input.get("num_solutions", 5);
  const auto include_Sigma2 = input.get("sigma2", false);
  const auto ci_input = input.get("ci_input", std::string{""});
  // double e0 = input.get("e0", 0.0) / PhysConst::Hartree_invcm;
  // even parity:
  std::vector<CIlevel> levels;
  for (auto &J : J_even_list) {
    const auto t_levels = run_CI(wf.atomicSymbol(), ci_sp_basis, orbital_map,
                                 int(std::round(2 * J)), +1, num_solutions, h1,
                                 qk, write_integrals, include_Sigma2,
                                 wf.basis(), wf.FermiLevel(), 1, ci_input);
    levels.insert(levels.begin(), t_levels.begin(), t_levels.end());
    // if (e0 == 0.0)
    // e0 = e1;
  }
  // odd parity:
  for (auto &J : J_odd_list) {
    const auto t_levels = run_CI(wf.atomicSymbol(), ci_sp_basis, orbital_map,
                                 int(std::round(2 * J)), -1, num_solutions, h1,
                                 qk, write_integrals, include_Sigma2,
                                 wf.basis(), wf.FermiLevel(), 1, ci_input);
    levels.insert(levels.begin(), t_levels.begin(), t_levels.end());
  }

  std::sort(levels.begin(), levels.end(),
            [](const auto &a, const auto &b) { return a.e < b.e; });

  std::cout << "\nLevel Summry:\n\n";
  const auto e0 = levels.at(0).e;

  std::cout
      << "config.   Jπ   Energy(au)    Energy(/cm)  Level(/cm)  Term    gJ\n";
  for (const auto &[config, twoj, parity, e, gJ, L, tSp1] : levels) {

    char pm = parity == 1 ? '+' : '-';

    fmt::print("{:<8s} {:>2}{}  {:+11.8f}  {:+11.2f}  "
               "{:11.2f}   {}^{}_{}{}",
               config, twoj / 2, pm, e, e * PhysConst::Hartree_invcm,
               (e - e0) * PhysConst::Hartree_invcm, uint(std::round(tSp1)),
               AtomData::L_symbol((int)std::round(L)), twoj / 2, pm);
    if (gJ != 0.0) {
      fmt::print("  {:.4f}", gJ);
    }
    std::cout << "\n";
  }
}

//==============================================================================
//==============================================================================

//==============================================================================
void write_CSFs(const std::vector<CI::CSF2> &CSFs, int twoJ,
                const std::map<nkm, int> &orbital_map,
                const std::string &csf_fname) {

  std::cout << "Writing CSFs and projections to files: {csf/proj}-" << csf_fname
            << "\n";
  std::ofstream csf_file("csf-" + csf_fname);
  std::ofstream proj_file("proj-" + csf_fname);

  csf_file << "# csf_index  a  b\n";
  proj_file << "# proj_index  a  b  CGC\n";
  int csf_count = 0;
  int proj_count = 0;
  for (const auto &csf : CSFs) {

    const auto &v = *csf.state(0);
    const auto &w = *csf.state(1);

    fmt::print(csf_file, "{} {} {}\n", csf_count, v.shortSymbol(),
               w.shortSymbol());
    ++csf_count;

    // Each individual m projection:
    for (int two_m_v = -v.twoj(); two_m_v <= v.twoj(); two_m_v += 2) {
      const auto two_m_w = twoJ - two_m_v;
      if (std::abs(two_m_w) > w.twoj())
        continue;
      const auto cgc =
          Angular::cg_2(v.twoj(), two_m_v, w.twoj(), two_m_w, twoJ, twoJ);
      if (cgc == 0.0)
        continue;

      const auto iv = orbital_map.at(nkm{v.n(), v.kappa(), two_m_v});
      const auto iw = orbital_map.at(nkm{w.n(), w.kappa(), two_m_w});
      if (iw > iv)
        continue;

      // if (v == w && two_m_v == two_m_w)
      //   continue;

      const auto eta = v == w ? 1 / std::sqrt(2.0) : 1.0;
      const auto d_proj = eta * cgc; //?

      fmt::print(proj_file, "{} {} {} {:.8f}\n", proj_count, iv, iw, d_proj);
      ++proj_count;
    }
  }
}

//==============================================================================
void write_CoulombIntegrals(const std::vector<DiracSpinor> &ci_sp_basis,
                            const std::map<nkm, int> &orbital_map,
                            const Coulomb::QkTable &qk) {

  // Modify this to include Sigma_2!

  std::string two_file_R =
      "h2_" + DiracSpinor::state_config(ci_sp_basis) + ".txt";

  std::cout << "\nWriting radial Coulomb intgrals 'g' to file: " << two_file_R
            << "\n";
  // print all two-particle integrals Rk to file:
  std::ofstream g_file(two_file_R);
  g_file << "# a  b  c  d  g_abcd    ## (nb: abcd = badc = cdab = dcba; only "
            "'smallest' is written)\n";

  for (const auto &a : ci_sp_basis) {
    for (const auto &b : ci_sp_basis) {
      for (const auto &c : ci_sp_basis) {
        for (const auto &d : ci_sp_basis) {

          // Parity and triangle j selection rules:
          const auto [k0, k1] = Coulomb::k_minmax_Q(a, b, c, d);
          if (k1 < k0)
            continue;

          for (int tma = -a.twoj(); tma <= a.twoj(); tma += 2) {
            for (int tmb = -b.twoj(); tmb <= b.twoj(); tmb += 2) {
              for (int tmc = -c.twoj(); tmc <= c.twoj(); tmc += 2) {
                for (int tmd = -d.twoj(); tmd <= d.twoj(); tmd += 2) {

                  // m = j_z selection rules:
                  if (tmc - tma != tmb - tmd)
                    continue;

                  const auto ia =
                      (uint16_t)orbital_map.at(nkm{a.n(), a.kappa(), tma});
                  const auto ib =
                      (uint16_t)orbital_map.at(nkm{b.n(), b.kappa(), tmb});
                  const auto ic =
                      (uint16_t)orbital_map.at(nkm{c.n(), c.kappa(), tmc});
                  const auto id =
                      (uint16_t)orbital_map.at(nkm{d.n(), d.kappa(), tmd});

                  // Equivilant integrals:
                  // abcd = badc = cdab = dcba
                  // nb: this only works if largest of (ia,ib,ic,id)
                  // is smaller than 2^16, which is always true
                  const auto indexify = [](uint16_t w, uint16_t x, uint16_t y,
                                           uint16_t z) {
                    return ((uint64_t)w << 48) + ((uint64_t)x << 32) +
                           ((uint64_t)y << 16) + (uint64_t)z;
                  };
                  uint64_t i1 = indexify(ia, ib, ic, id);
                  uint64_t i2 = indexify(ib, ia, id, ic);
                  uint64_t i3 = indexify(ic, id, ia, ib);
                  uint64_t i4 = indexify(id, ic, ib, ia);
                  // Only include the unique ones:
                  if (i1 != std::min({i1, i2, i3, i4}))
                    continue;
                  const auto g = qk.g(a, b, c, d, tma, tmb, tmc, tmc);
                  // Note sure why zero values slip through?
                  // Missing SR? or mistake?
                  if (g == 0.0)
                    continue;

                  fmt::print(g_file, "{} {} {} {} {:.8e}\n", ia, ib, ic, id, g);
                  // print extra detail:
                  // fmt::print(
                  //     g_file, "{} {} {} {} {:.8e}  {} {} {} {}  {} {} {} {}\n",
                  //     ia, ib, ic, id, g, a.shortSymbol(), b.shortSymbol(),
                  //     c.shortSymbol(), d.shortSymbol(), tma, tmb, tmc, tmd);
                }
              }
            }
          }
        }
      }
    }
  }
}

//==============================================================================
std::vector<CIlevel>
run_CI(const std::string &atom_name,
       const std::vector<DiracSpinor> &ci_sp_basis,
       const std::map<nkm, int> &orbital_map, int twoJ, int parity,
       int num_solutions, const Coulomb::meTable<double> &h1,
       const Coulomb::QkTable &qk, bool write_integrals, bool include_Sigma2,
       const std::vector<DiracSpinor> &mbpt_basis, double E_Fermi, int min_n,
       const std::string &ci_input) {
  //----------------------------------------------------------------------------

  std::vector<CIlevel> out;

  auto printJ = [](int twoj) {
    return twoj % 2 == 0 ? std::to_string(twoj / 2) :
                           std::to_string(twoj) + "/2";
  };
  auto printPi = [](int pi) { return pi > 0 ? "even" : "odd"; };

  fmt::print("\nForm CSFs for J={}, {} parity\n", printJ(twoJ),
             printPi(parity));
  std::cout << std::flush;

  if (twoJ < 0) {
    std::cout << "twoJ must >=0\n";
    return out;
  }
  if (twoJ % 2 != 0) {
    std::cout << "twoJ must be even for two-electron CSF\n";
  }

  std::vector<CI::CSF2> CSFs = CI::form_CSFs(twoJ, parity, ci_sp_basis);
  std::cout << "Total CSFs: " << CSFs.size() << "\n";
  std::cout << std::flush;

  //----------------------------------------------------------------------------

  std::string output_prefix =
      atom_name + "_" + printJ(twoJ) + "_" + printPi(parity);

  // Write CSFs (just labels) to file:
  if (write_integrals)
    write_CSFs(CSFs, twoJ, orbital_map, output_prefix + ".txt");

  //----------------------------------------------------------------------------
  fmt::print("Construct CI matrix for J={}, {} parity:\n", printJ(twoJ),
             printPi(parity));
  std::cout << std::flush;

  LinAlg::Matrix Hci(CSFs.size(), CSFs.size());

  {
    IO::ChronoTimer t("Fill matrix");
#pragma omp parallel for collapse(2)
    for (std::size_t iA = 0; iA < CSFs.size(); ++iA) {
      // go to iB <= iA only: symmetric matrix
      // for (std::size_t iB = 0; iB <= iA; ++iB) { // work with collapse? how?
      for (std::size_t iB = 0; iB < CSFs.size(); ++iB) {
        if (iB > iA)
          continue;
        const auto &A = CSFs.at(iA);
        const auto &B = CSFs.at(iB);

        const auto E_AB = Hab(A, B, twoJ, h1, qk);
        Hci(iA, iB) = E_AB;
        // fill other half of symmetric matrix:
        if (iB != iA) {
          Hci(iB, iA) = E_AB;
        }
      }
    }
  }
  std::cout << std::flush;

  if (include_Sigma2 && !mbpt_basis.empty()) {
    LinAlg::Matrix H_sigma(CSFs.size(), CSFs.size());

    const auto [core, excited] = MBPT::split_basis(mbpt_basis, E_Fermi, 1);

    Angular::SixJTable sjt(DiracSpinor::max_tj(mbpt_basis));

    IO::ChronoTimer t("Add Sigma Sigma matrix");
#pragma omp parallel for
    for (std::size_t iA = 0; iA < CSFs.size(); ++iA) {
      for (std::size_t iB = 0; iB <= iA; ++iB) {
        // go to iB <= iA only: symmetric matrix

        const auto &A = CSFs.at(iA);
        const auto &B = CSFs.at(iB);

        const auto dE_AB = CI::Sigma2_AB(A, B, twoJ, qk, core, excited, sjt);
        H_sigma(iA, iB) = dE_AB;
        // Add to other half of symmetric matrix:
        if (iB != iA) {
          H_sigma(iB, iA) = dE_AB;
        }
      }
    }
    Hci += H_sigma;
  }
  std::cout << std::flush;

  // Write CI matrix (H matrix) to file
  if (write_integrals) {
    std::string ci_fname = "ci-" + output_prefix + ".txt";
    std::cout << "Writing CI matrix to file: " << ci_fname << "\n";
    std::ofstream ci_file(ci_fname);
    ci_file << "# in matrix/table form: \n";
    for (std::size_t iA = 0; iA < CSFs.size(); ++iA) {
      for (std::size_t iX = 0; iX < CSFs.size(); ++iX) {
        fmt::print(ci_file, "{:+.6e} ", Hci(iA, iX));
      }
      ci_file << "\n";
    }
  }

  //----------------------------------------------------------------------------
  std::cout << std::flush;

  IO::ChronoTimer t2("Diagonalise:");
  // const auto [val, vec] = LinAlg::symmhEigensystem(Hci);
  const auto [val, vec] = LinAlg::symmhEigensystem(Hci, num_solutions);
  std::cout << "T=" << t2.lap_reading_str() << "\n";
  const auto E0 = val(0);

  fmt::print("Full CI for J={}, pi={} : E0 = {:.1f} cm^-1\n\n", printJ(twoJ),
             printPi(parity), E0 * PhysConst::Hartree_invcm);
  std::cout << std::flush;

  // For calculating g-factors
  // (Use non-rel formula? Or relativistic M1?)
  DiracOperator::M1nr m1{};
  // DiracOperator::M1 m1{wf.grid(), wf.alpha(), 0.0};

  int l1, l2;

  for (std::size_t i = 0; i < val.size() && int(i) < num_solutions; ++i) {

    fmt::print(
        "{:<2} {} {:+1}  {:+11.8f} au  {:+11.2f} cm^-1  {:11.2f} cm^-1\n", i,
        0.5 * twoJ, parity, val(i), val(i) * PhysConst::Hartree_invcm,
        (val(i) - E0) * PhysConst::Hartree_invcm);

    std::size_t max_j = 0;
    double max_cj = 0.0;
    for (std::size_t j = 0ul; j < vec.cols(); ++j) {
      const auto cj = 100.0 * std::pow(vec(i, j), 2);
      if (cj > max_cj) {
        max_cj = cj;
        max_j = j;
        l1 = CSFs.at(j).state(0)->l();
        l2 = CSFs.at(j).state(1)->l();
      }
      if (cj > 1.0) {
        fmt::print("  {:>4s},{:<4s} {:5.3f}%\n",
                   CSFs.at(j).state(0)->shortSymbol(),
                   CSFs.at(j).state(1)->shortSymbol(), cj);
      }
    }

    // Calculate g-factors, for line identification. Only defined for J!=0
    double gJ = 0.0;

    // g_J <JJz|J|JJz> = <JJz|L + 2*S|JJz>
    // take J=Jz, <JJz|J|JJz> = J
    // then: g_J = <JJ|L + 2*S|JJ> / J
    // And: <JJ|L + 2*S|JJ> = 3js * <A||L+2S||A> (W.E. Theorem)
    const auto m1AA = CI::ReducedME(vec.row(i), CSFs, twoJ, &m1);
    const auto tjs = Angular::threej_2(twoJ, twoJ, 2, twoJ, -twoJ, 0);

    if (twoJ != 0) {
      gJ = tjs * m1AA / (0.5 * twoJ);
      std::cout << "gJ = " << gJ << "\n";
    }

    // Determine Term Symbol, from g-factor
    const auto min_L = std::abs(l1 - l2);
    const auto max_L = std::abs(l1 + l2);
    const auto min_S = 0;
    const auto max_S = 1;
    int L = min_L;
    int S = min_L;
    double best_del = 2.0;
    if (twoJ != 0) {
      for (int tL = min_L; tL <= max_L; ++tL) {
        for (int tS = min_S; tS <= max_S; ++tS) {
          auto gJNR = 1.5 + (tS * (tS + 1.0) - tL * (tL + 1.0)) /
                                (twoJ * (0.5 * twoJ + 1.0));
          if (std::abs(gJ - gJNR) < best_del) {
            best_del = std::abs(gJ - gJNR);
            L = tL;
            S = tS;
          }
        }
      }
    }
    const auto tSp1 = 2.0 * S + 1.0;

    const auto config =
        fmt::format("{:s},{:s}", CSFs.at(max_j).state(0)->shortSymbol(),
                    CSFs.at(max_j).state(1)->shortSymbol());
    out.emplace_back(
        CIlevel{config, twoJ, parity, val(i), gJ, double(L), tSp1});

    std::cout << "\n";
  }

  //----------------------------------------------------------------------------
  // std::cout
  //     << "\n`Direct' energy calculation: E = Σ_{IJ} c_I * c_J * <I|H|J>:\n";
  // std::cout
  //     << "(Using the CI expansion coefficients from full CI, just a test)\n";
  // // Energy calculation for ground state:
  // double E_direct1 = 0.0, E_direct2 = 0.0;
  // const auto Nci = vec.rows(); // number of CSFs
  // // Energy:  E = Sum_ij c_i * c_j * <i|H|j>
  // for (std::size_t i = 0ul; i < Nci; ++i) {
  //   const auto &csf_i = CSFs.at(i); // the ith CSF
  //   const auto ci = vec.at(0, i);   // the ith CI coefficient (for 0th e.val)
  //   for (std::size_t j = 0ul; j < Nci; ++j) {
  //     const auto &csf_j = CSFs.at(j); // jth CSF
  //     const auto cj = vec.at(0, j);   // the jth CI coefficient (for 0th e.val)
  //     // use pre-calculated CI matrix:
  //     E_direct1 += ci * cj * Hci.at(i, j);
  //     // Calculate MEs on-the-fly
  //     E_direct2 += ci * cj * Hab(csf_i, csf_j, twoJ, h1, qk);
  //   }
  // }

  // std::cout << "E0 = " << val.at(0) * PhysConst::Hartree_invcm
  //           << " cm^-1  (from diagonalisation)\n";
  // std::cout << "E0 = " << E_direct1 * PhysConst::Hartree_invcm
  //           << " cm^-1  (uses pre-calculated CI matrix)\n";
  // std::cout << "E0 = " << E_direct2 * PhysConst::Hartree_invcm
  //           << " cm^-1  (calculates H matrix elements from scratch)\n";

  // // std::ifstream ci("")
  // std::ifstream is(ci_input);
  // if (is) {
  //   std::istream_iterator<double> start(is), end;
  //   std::vector<double> in_CI(start, end);
  //   std::cout << "\nCalculate energy from input CI coeficients:\n";
  //   std::cout << "Read " << in_CI.size() << " CI coeficients from: " << ci_input
  //             << "\n";
  //   double E_input = 0.0;
  //   for (std::size_t i = 0ul; i < std::min(Nci, in_CI.size()); ++i) {
  //     const auto &csf_i = CSFs.at(i); // the ith CSF
  //     const auto ci = in_CI.at(i);    // the ith CI coefficient (for 0th e.val)
  //     for (std::size_t j = 0ul; j < std::min(Nci, in_CI.size()); ++j) {
  //       const auto &csf_j = CSFs.at(j); // jth CSF
  //       const auto cj = in_CI.at(j); // the jth CI coefficient (for 0th e.val)
  //       E_input += ci * cj * Hab(csf_i, csf_j, twoJ, h1, qk);
  //     }
  //   }
  //   std::cout << "E = " << E_input * PhysConst::Hartree_invcm << "\n";
  // }

  return out;
}

} // namespace Module
