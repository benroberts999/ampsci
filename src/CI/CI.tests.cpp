#include "Amplitudes/MatrixElements.hpp"
#include "Angular/Wigner369j.hpp"
#include "CI_Integrals.hpp"
#include "ConfigurationInteraction.hpp"
#include "ConfigurationSpace.hpp"
#include "Coulomb/include.hpp"
#include "DiracOperator/include.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/InputBlock.hpp"
#include "LinAlg/include.hpp"
#include "MBPT/Sigma2.hpp"
#include "MixedStates.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Random.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <tuple>
#include <variant>
#include <vector>

namespace UnitTest {
// CI tables for Li (core [He]): orbitals 1s and the excited 5spd basis
// states; Coulomb Q^k; one-body (orbital energies, arbitrary off-diagonal
// elements); Breit and Sigma_2 (arbitrary multiples of Q^k; Sigma_2 over
// n <= 3, with extrapolation ratios and a dS^k/dE0 table); a dSigma/dE
// correction. As in "CI: three-electron Hamiltonian"
struct CITables {
  std::vector<DiracSpinor> orbitals{};
  Coulomb::QkTable qk{};
  Coulomb::meTable<double> h1{};
  Coulomb::WkTable Bk{};
  Coulomb::WkTable Sk{};
  Coulomb::WkTable dSk{};
  std::vector<double> hk{};
  CI::Sigma1Correction s1c{};
  // Shift of E0 from where Sk was tabulated
  double dE0{0.0};
};
CITables Li_CI_tables();
// CI input options for Be (core [He]), Sigma_2 in a small basis, two levels
// of 0+ and 1-, no files written
std::string Be_sigma2_input();

// Slater determinants, for the hole CSF tests. The spin-orbitals (nk_index,
// 2m) of a list of orbitals are numbered in order; a determinant is a sorted
// list of spin-orbital numbers, |D> = a+_D[0] a+_D[1] ... |0>; a state is a
// map from determinants to amplitudes
using SpinOrbital = std::pair<DiracSpinor::Index, int>;
using Determinant = std::vector<std::size_t>;
using Expansion = std::map<Determinant, double>;

// A small atom with hole CSFs: the Hartree-Fock core, the core orbitals
// that may hold a hole, the excited orbitals (the CI basis), all orbitals
// (core first), and Coulomb Q^k and Breit B^k tables over all of them
struct HoleAtom {
  Wavefunction wf;
  std::vector<DiracSpinor> core{}, holes{}, excited{}, all{};
  Coulomb::QkTable qk{};
  Coulomb::WkTable Bk{};
};
HoleAtom hole_atom(const std::string &atom, const std::string &core,
                   const std::vector<std::string> &hole_shells,
                   const std::string &basis);
// Li (1s hole), Be (1s and 2s holes), Na (2s and 2p holes), with their
// (2J, parity) blocks
struct HoleAtomInput {
  std::string atom, core;
  std::vector<std::string> hole_shells;
  std::string basis;
  std::vector<std::pair<int, int>> blocks;
};
std::vector<HoleAtomInput> hole_atoms();

// Spin-orbitals of the orbitals, in order; the number of the one (nk, 2m)
std::vector<SpinOrbital>
spin_orbitals(const std::vector<DiracSpinor> &orbitals);
std::size_t spin_orbital(const std::vector<SpinOrbital> &so,
                         DiracSpinor::Index nk, int twom);
// Orbital energy of each spin-orbital
std::vector<double> energies(const std::vector<SpinOrbital> &so,
                             const std::vector<DiracSpinor> &orbitals);
// a+_i |D> and a_i |D>: the sign (zero if not possible) and the result
std::pair<double, Determinant> create(std::size_t i, Determinant D);
std::pair<double, Determinant> annihilate(std::size_t i, Determinant D);
// <ij|g|kl> between spin-orbitals, as g(i n + j, k n + l): Coulomb from the
// Q^k table, plus Breit from the B^k table if given
LinAlg::Matrix<double> g_table(const std::vector<SpinOrbital> &so,
                               const std::vector<DiracSpinor> &orbitals,
                               const Coulomb::QkTable &qk,
                               const Coulomb::WkTable *Bk = nullptr);
// <D'|H|D> for H = sum_i h0(i) + sum_{i<j} g(i,j), by the Slater-Condon
// rules, with h0 = e - V_HF and V_HF from the core determinant: h0 + V_HF =
// e exactly, as the CI assumes (V^{N-M})
double H_det(const Determinant &Dp, const Determinant &D,
             const Determinant &core, const std::vector<double> &e,
             const LinAlg::Matrix<double> &g);
// A CSF (one electron, or two particles and a hole) at M = J
Expansion expand(const CI::CSF &csf, int twoJ, const Determinant &core,
                 const std::vector<SpinOrbital> &so);
// A particle-hole pair of the core, |(w a~) K, M = K> = sum C(j_w m_w, j_a m
// | K K) a+_{w m_w} h+_{a m} |0>, with h+_{a m} = (-1)^(j_a - m) a_{a,-m}
Expansion particle_hole(DiracSpinor::Index w, DiracSpinor::Index a, int twoK,
                        const Determinant &core,
                        const std::vector<SpinOrbital> &so);
// Matrix <A_i|op|B_j> of a determinant operator between two sets of states
template <class Op>
LinAlg::Matrix<double> matrix(const std::vector<Expansion> &A,
                              const std::vector<Expansion> &B, const Op &op) {
  LinAlg::Matrix<double> M(A.size(), B.size());
  for (std::size_t i = 0; i < A.size(); ++i) {
    for (std::size_t j = 0; j < B.size(); ++j) {
      double sum = 0.0;
      for (const auto &[Dp, amp_i] : A[i]) {
        for (const auto &[D, amp_j] : B[j]) {
          sum += amp_i * amp_j * op(Dp, D);
        }
      }
      M(i, j) = sum;
    }
  }
  return M;
}
// The determinants with the given number of electrons (particles minus
// holes), M and parity, outside the CI orbitals: one to three holes, at least
// one of them in a core shell not in `holes` (every excited orbital is a CI
// orbital here)
std::vector<Determinant> Q_space(const std::vector<SpinOrbital> &so,
                                 const Determinant &core,
                                 const std::vector<DiracSpinor> &holes,
                                 int twoM, int parity, int num_electrons);
// Second order between the states, through the determinants Q:
// (1/2)[1/(E_A - E_Q) + 1/(E_B - E_Q)] <A|H|Q><Q|H|B>, with orbital-sum
// energies (the symmetrised Rayleigh-Schrodinger form of the S^k)
LinAlg::Matrix<double> second_order(const std::vector<Expansion> &states,
                                    const std::vector<Determinant> &Q,
                                    const Determinant &core,
                                    const std::vector<double> &e,
                                    const LinAlg::Matrix<double> &g);
} // namespace UnitTest

//==============================================================================
TEST_CASE("CI: Configuration Interaction unit tests", "[CI][unit]") {

  std::cout << "CI, unit tests (not meant to be accurate)\n";

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"He", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[]", std::nullopt, 1.0e-5);
  wf.formBasis(SplineBasis::Parameters("6spd", 20, 6, {1.0e-2}, 20.0));

  std::string qk_filename = "deleteme_" + qip::random_string(3) + ".qk.abf";
  std::string ci_filename = "deleteme_" + qip::random_string(3) + ".ci.abf";

  const IO::InputBlock input{"CI", "ci_basis = 6spd;"
                                   "n_min_core = 1;"
                                   "J + = 0, 1;"
                                   "J - = 0, 1;"
                                   "qk_file = " +
                                     qk_filename +
                                     ";"
                                     "num_solutions = 2;"
                                     "ci_file = " +
                                     ci_filename + ";"};

  // Initial CI calculation (writes to file):
  const auto CIWFs = CI::configuration_interaction(input, wf).levels;

  // Expected values: (not accurate, just simple unit test)
  std::vector J{0, 1, 0, 1};
  std::vector parity{1, 1, -1, -1};
  std::vector energy{std::vector{-2.84834207, -2.13740814},
                     {-2.17164001, -2.0657259},
                     {-2.13100411, -2.05532175},
                     {-2.13100096, -2.12092003}};
  std::vector gj{
    std::vector{0.0, 0.0}, {1.9999, 1.9999}, {0.0, 0.0}, {1.5, 1.0}};
  // Terms: {1S0, 1S0}, {3S1, 3S1}, {3P0, 3P0}, {3P1, 1P1}
  std::vector L2{std::vector{0.0, 0.0}, {0.0, 0.0}, {2.0, 2.0}, {2.0, 2.0}};
  std::vector S2{std::vector{0.0, 0.0}, {2.0, 2.0}, {2.0, 2.0}, {2.0, 0.0}};

  for (std::size_t i = 0; i < CIWFs.size(); ++i) {
    const auto &ci_wf = CIWFs.at(i);
    REQUIRE(ci_wf.twoJ() == 2 * J[i]);
    REQUIRE(ci_wf.parity() == parity[i]);
    for (std::size_t j = 0ul; j < ci_wf.num_solutions(); ++j) {
      REQUIRE(ci_wf.energy(j) == Approx(energy[i][j]).epsilon(1.0e-2));
      REQUIRE(ci_wf.info(j).gJ == Approx(gj[i][j]).epsilon(1.0e-2));
      REQUIRE(ci_wf.info(j).L2 == Approx(L2[i][j]).margin(0.05));
      REQUIRE(ci_wf.info(j).S2 == Approx(S2[i][j]).margin(0.05));
    }
  }

  //-----------------------------------------------------------------------
  std::cout << "\nRe-run CI, reading in from file:\n";

  // Re-"run" identical CI calculation
  // This time: should read in prev solution from disk
  // Read from ci file; compare energies, gJ, and coefficients
  auto input2 = input;
  input2.add("read_only = true;");
  const auto solution_2 = CI::configuration_interaction(input2, wf);
  const auto CIWFs_2 = solution_2.levels;

  // read_only must still read in the integrals: the modules need them to
  // re-construct the CI Hamiltonian of any other J/parity
  REQUIRE(solution_2.integrals.availableQ());
  REQUIRE(solution_2.integrals.qk.count() ==
          CI::configuration_interaction(input, wf).integrals.qk.count());

  REQUIRE(CIWFs.size() == CIWFs_2.size());

  for (std::size_t i = 0; i < CIWFs.size(); ++i) {
    const auto &w = CIWFs.at(i);
    const auto &r = CIWFs_2.at(i);
    REQUIRE(w.twoJ() == r.twoJ());
    REQUIRE(w.parity() == r.parity());
    REQUIRE(w.num_solutions() == r.num_solutions());
    REQUIRE(w.CSFs().size() == r.CSFs().size());
    for (std::size_t s = 0; s < w.num_solutions(); ++s) {
      REQUIRE(w.energy(s) == Approx(r.energy(s)));
      REQUIRE(w.info(s).gJ == Approx(r.info(s).gJ));
      const auto wc = w.coefs(s);
      const auto rc = r.coefs(s);
      REQUIRE(wc.size() == rc.size());
      REQUIRE(wc.size() == w.CSFs().size());
      for (std::size_t k = 0; k < wc.size(); ++k) {
        REQUIRE(wc[k] == Approx(rc[k]));
      }
    }
  }

  //-----------------------------------------------------------------------
  std::cout << "\nRe-run CI, reading in from file: but ask for more states\n";
  // Now, request 3 solutions:
  // ci file has 2, so must re-solve and overwrite
  const IO::InputBlock input3{"CI", "ci_basis = 6spd;"
                                    "n_min_core = 1;"
                                    "J + = 0, 1;"
                                    "J - = 0, 1;"
                                    "qk_file = " +
                                      qk_filename +
                                      ";"
                                      "num_solutions = 3;"
                                      "print_details = false;"
                                      "ci_file = " +
                                      ci_filename + ";"};
  // don't read in
  const IO::InputBlock input4{"CI", "ci_basis = 6spd;"
                                    "n_min_core = 1;"
                                    "J + = 0, 1;"
                                    "J - = 0, 1;"
                                    "qk_file = false;"
                                    "num_solutions = 3;"
                                    "print_details = false;"
                                    "ci_file = " +
                                      ci_filename + ";"};

  // Reads from file (2 stored, 3 required, re-calculates)
  const auto CIWFs_3a = CI::configuration_interaction(input3, wf).levels;
  // calculate 3 from scratch
  const auto CIWFs_3b = CI::configuration_interaction(input3, wf).levels;

  REQUIRE(CIWFs_3a.size() == CIWFs.size());
  REQUIRE(CIWFs_3b.size() == CIWFs.size());
  for (std::size_t i = 0; i < CIWFs_3a.size(); ++i) {
    const auto &a = CIWFs_3a.at(i);
    const auto &b = CIWFs_3b.at(i);
    REQUIRE(a.num_solutions() == 3);
    REQUIRE(b.num_solutions() == 3);
    // First 2 solutions match CIWFs
    for (std::size_t s = 0; s < 2; ++s) {
      REQUIRE(a.energy(s) == Approx(CIWFs.at(i).energy(s)));
      REQUIRE(a.info(s).gJ == Approx(CIWFs.at(i).info(s).gJ));
    }
    // All 3 identical between solve and read
    for (std::size_t s = 0; s < 3; ++s) {
      REQUIRE(a.energy(s) == Approx(b.energy(s)));
      REQUIRE(a.info(s).gJ == Approx(b.info(s).gJ));
      const auto ac = a.coefs(s);
      const auto bc = b.coefs(s);
      REQUIRE(ac.size() == bc.size());
      for (std::size_t k = 0; k < ac.size(); ++k) {
        REQUIRE(ac[k] == Approx(bc[k]));
      }
    }
  }

  //-----------------------------------------------------------------------
  // Sigma_2 extrapolation

  {
    const auto &ints = solution_2.integrals;

    // Table with S^k = Q^k for a small subset: ratio S/Q = 1 exactly
    Coulomb::WkTable Sk_test;
    const auto S_eq_Q = [&ints](int k, const DiracSpinor &a,
                                const DiracSpinor &b, const DiracSpinor &c,
                                const DiracSpinor &d) {
      return ints.qk.Q(k, a, b, c, d);
    };
    const auto Qk_SR = [](int k, const DiracSpinor &a, const DiracSpinor &b,
                          const DiracSpinor &c, const DiracSpinor &d) {
      return Coulomb::Qk_abcd_SR(k, a, b, c, d);
    };
    const auto sub_basis = CI::basis_subset(ints.ci_basis, "3sp");
    Sk_test.fill(sub_basis, S_eq_Q, Qk_SR, 8, false);

    // Average correction ratios: exactly 1 for every k with data
    const auto hk = MBPT::average_hk(Sk_test, ints.qk, sub_basis, 8);
    bool any_data = false;
    for (const auto &h : hk) {
      if (h != 0.0) {
        any_data = true;
        REQUIRE(h == Approx(1.0));
      }
    }
    REQUIRE(any_data);

    // Extrapolate with hk = 1 for all k: then S^k = Q^k everywhere, and
    // Sigma2_AB must equal the Coulomb part for every CSF pair.
    // nb: Sk_test only covers sub_basis, so all other pairs are extrapolated
    const std::vector<double> hk_1(10, 1.0);
    for (const auto &psi0 : CIWFs) {
      const auto twoJ = psi0.twoJ();
      const auto n_csf = std::min(psi0.CSFs().size(), std::size_t(6));
      for (std::size_t i = 0; i < n_csf; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
          const auto &A = psi0.csf(i);
          const auto &B = psi0.csf(j);
          const auto s2_x = CI::Sigma2_AB(A, B, twoJ, Sk_test, &ints.qk, hk_1);
          const auto coulomb = CI::CSF2_Coulomb(ints.qk, A.state(0), A.state(1),
                                                B.state(0), B.state(1), twoJ);
          REQUIRE(s2_x == Approx(coulomb).margin(1.0e-12));
        }
      }
    }
  }

  //-----------------------------------------------------------------------
  // Multiple (J, pi) blocks in one ci solutions file

  {
    const auto fname = "deleteme_" + qip::random_string(3) + ".ci.abf";
    const auto t_basis = CI::basis_subset(wf.basis(), "4sp");
    CI::PsiJPi psi_a{0, +1, t_basis};
    psi_a.set_solution(-1.0, LinAlg::Vector<double>(psi_a.CSFs().size()));
    REQUIRE(psi_a.read_write(fname, IO::FRW::write, std::cout));
    // second block is appended; the first must survive
    CI::PsiJPi psi_b{2, +1, t_basis};
    psi_b.set_solution(-0.5, LinAlg::Vector<double>(psi_b.CSFs().size()));
    REQUIRE(psi_b.read_write(fname, IO::FRW::write, std::cout));
    CI::PsiJPi read_a{0, +1, t_basis};
    REQUIRE(read_a.read_write(fname, IO::FRW::read, std::cout));
    REQUIRE(read_a.energy(0) == Approx(-1.0));
    CI::PsiJPi read_b{2, +1, t_basis};
    REQUIRE(read_b.read_write(fname, IO::FRW::read, std::cout));
    REQUIRE(read_b.energy(0) == Approx(-0.5));
    // a block that was never written:
    CI::PsiJPi read_c{0, -1, t_basis};
    REQUIRE(!read_c.read_write(fname, IO::FRW::read, std::cout));
  }

  //-----------------------------------------------------------------------
  // basic/Misc tests

  // Term(int two_J, int L, int two_S, int parity)

  REQUIRE(CI::Term_Symbol(2, 3, 2, +1) == "3^F_1");
  REQUIRE(CI::Term_Symbol(2, 3, 2, -1) == "3^F°_1");
  REQUIRE(CI::Term_Symbol(1, 3, 2, -1) == "3^F°_1/2");
  REQUIRE(CI::Term_Symbol(1, 2, 1, -1) == "2^D°_1/2");
  REQUIRE(CI::Term_Symbol(6, 0, 0, 1) == "1^S_3");
}

//==============================================================================
TEST_CASE("CI: derivative correction", "[CI][unit]") {

  std::cout << "CI, derivative (dSigma/dE) correction for Sigma_1\n";

  //-----------------------------------------------------------------------
  // corrected_Sigma: formula and guard

  // Sigma = 0: no correction
  REQUIRE(CI::corrected_Sigma(0.0, 0.1, 0.5) == 0.0);

  // dE = 0: Sigma unchanged (Sigma^2/Sigma = Sigma)
  REQUIRE(CI::corrected_Sigma(-0.1, 0.05, 0.0) == Approx(-0.1).margin(1.0e-14));

  // Exact resummed value: Sigma*Sigma/(Sigma - dE*dSigma)
  REQUIRE(CI::corrected_Sigma(-0.1, -0.02, -0.5) ==
          Approx(0.01 / -0.11).margin(1.0e-14));

  // Small dE: matches linear expansion, Sigma + dE*dSigma
  {
    const auto dE = 1.0e-4;
    const auto lin = -0.1 + dE * 0.05;
    REQUIRE(CI::corrected_Sigma(-0.1, 0.05, dE) == Approx(lin).margin(1.0e-8));
  }

  // Guard: corrected value exceeds |Sigma| (near divergence): uncorrected
  // -0.1 - (-0.5)*(0.02) = -0.09, so |Sig^2/denom| = 0.111 > 0.1
  REQUIRE(CI::corrected_Sigma(-0.1, 0.02, -0.5) ==
          Approx(-0.1).margin(1.0e-14));

  //-----------------------------------------------------------------------
  // corrected_Sk: formula and pole guard

  // Sk = 0: no shift
  REQUIRE(CI::corrected_Sk(0.0, 0.1, 0.5) == 0.0);

  // dE0 = 0: Sk unchanged (Sk^2/Sk = Sk)
  REQUIRE(CI::corrected_Sk(-0.1, 0.05, 0.0) == Approx(-0.1).margin(1.0e-14));

  // Exact resummed value: Sk*Sk/(Sk - dE0*dSk)
  REQUIRE(CI::corrected_Sk(-0.1, -0.02, -0.5) ==
          Approx(0.01 / -0.11).margin(1.0e-14));

  // Unlike corrected_Sigma, enhancement is allowed:
  // 0.1 - 0.5*0.05 = 0.075 (same sign as Sk): |corrected| = 0.133 > 0.1
  REQUIRE(CI::corrected_Sk(0.1, 0.05, 0.5) ==
          Approx(0.01 / 0.075).margin(1.0e-14));

  // Small dE0: matches linear expansion, Sk + dE0*dSk
  {
    const auto dE0 = 1.0e-4;
    const auto lin = -0.1 + dE0 * 0.05;
    REQUIRE(CI::corrected_Sk(-0.1, 0.05, dE0) == Approx(lin).margin(1.0e-8));
  }

  // Exact for a single denominator: Sk = N/D, dSk = -N/D^2
  // => corrected = N/(D + dE0)
  {
    const auto N = 0.35, D = -1.2, dE0 = 0.4;
    const auto Sk = N / D;
    const auto dSk = -N / (D * D);
    REQUIRE(CI::corrected_Sk(Sk, dSk, dE0) ==
            Approx(N / (D + dE0)).margin(1.0e-14));
  }

  // Pole guard: resummed denominator changes sign (dE0 past the pole at
  // Sk/dSk = 2.0): return unshifted
  REQUIRE(CI::corrected_Sk(0.1, 0.05, 3.0) == Approx(0.1).margin(1.0e-14));
  // Exactly at the pole:
  REQUIRE(CI::corrected_Sk(0.1, 0.05, 2.0) == Approx(0.1).margin(1.0e-14));

  //-----------------------------------------------------------------------
  // Sigma1Correction::delta

  const auto ia = static_cast<DiracSpinor::Index>(Angular::nk_to_index(1, -1));
  const auto ib = static_cast<DiracSpinor::Index>(Angular::nk_to_index(2, -1));

  CI::Sigma1Correction corr;
  corr.E0 = -2.0;
  corr.e_sigma[-1] = -1.0;
  corr.en[ia] = -1.0;
  corr.en[ib] = -0.5;
  corr.S1.add(ia, ia, -0.1);
  corr.S1.add(ia, ib, -0.05);
  corr.S1.add(ib, ia, -0.05);
  corr.S1.add(ib, ib, -0.08);
  corr.dS1.add(ia, ia, -0.02);
  corr.dS1.add(ia, ib, -0.02);
  corr.dS1.add(ib, ia, -0.02);
  corr.dS1.add(ib, ib, -0.02);

  REQUIRE(!corr.empty());

  // delta(a,a; spectator b): dE = E0 - e_sigma - en(b) = -2 + 1 + 0.5 = -0.5
  // corrected = 0.01/(-0.1 - 0.01) = -0.0909..., delta = +0.00909...
  REQUIRE(corr.delta_h1(ia, ia, ib) ==
          Approx(0.01 / -0.11 + 0.1).margin(1.0e-14));

  // Symmetric in (a, b):
  REQUIRE(corr.delta_h1(ia, ib, ia) == Approx(corr.delta_h1(ib, ia, ia)));

  // Missing orbital/spectator: no correction
  const auto ix = static_cast<DiracSpinor::Index>(Angular::nk_to_index(3, 1));
  REQUIRE(corr.delta_h1(ix, ix, ia) == 0.0);
  REQUIRE(corr.delta_h1(ia, ia, ix) == 0.0);

  // Zero derivative: no correction (any dE)
  CI::Sigma1Correction corr0 = corr;
  corr0.dS1 = {};
  corr0.dS1.add(ia, ia, 0.0);
  REQUIRE(corr0.delta_h1(ia, ia, ib) == Approx(0.0).margin(1.0e-14));

  // Default (empty) correction:
  REQUIRE(CI::Sigma1Correction{}.empty());

  //-----------------------------------------------------------------------
  // read_write round-trip

  const auto ds1_fname = "deleteme_" + qip::random_string(3) + ".ds1.abf";
  REQUIRE(corr.read_write(ds1_fname, IO::FRW::write));

  CI::Sigma1Correction corr2;
  corr2.E0 = corr.E0; // E0 is not stored (per J/pi; set at solve time)
  REQUIRE(corr2.read_write(ds1_fname, IO::FRW::read));

  REQUIRE(corr2.S1.getv(ia, ib) == corr.S1.getv(ia, ib));
  REQUIRE(corr2.dS1.getv(ib, ib) == corr.dS1.getv(ib, ib));
  REQUIRE(corr2.e_sigma == corr.e_sigma);
  REQUIRE(corr2.en == corr.en);
  REQUIRE(corr2.delta_h1(ia, ia, ib) == corr.delta_h1(ia, ia, ib));

  // Missing file: read fails
  REQUIRE(!CI::Sigma1Correction{}.read_write("deleteme_no_such_file.abf",
                                             IO::FRW::read));
}

//==============================================================================
// Tests the Sigma_1 correction tables built from a correlation potential
// holding dSigma/dE, against the second-order (qk-route) tables
TEST_CASE("CI: dSdE from correlation potential", "[CI][integration]") {

  std::cout << "CI, dSigma/dE correction tables from correlation potential\n";

  Wavefunction wf({800, 1.0e-5, 40.0, 10.0, "loglinear"},
                  {"Na", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-7);
  wf.solve_valence("3sp");
  // wf.formBasis(SplineBasis::Parameters("12spd", 25, 7, {1.1e-3, 3.2e-2}, 30.0));
  wf.formBasis(SplineBasis::Parameters("12spd", 25, 7, {1.0e-3}, 30.0));

  // Goldstone correlation potential, with dSigma/dE (derivative = true);
  // written to file for the round-trip check below
  const auto sig_fname = "deleteme_" + qip::random_string(3) + ".s2.abf";
  wf.formSigma(1, 1.0e-4, 30.0, 2, false, false, false, 0, {}, {}, {}, true,
               sig_fname, false, false, false, 6, -0.2, 0.01, 1.5, false,
               std::nullopt, "", true);
  REQUIRE(wf.Sigma() != nullptr);
  REQUIRE(wf.Sigma()->has_derivative());

  // Tables from the correlation potential:
  const auto ci_basis =
    CI::basis_subset(wf.basis(), "5sp", wf.coreConfiguration());
  const auto corr_CP = CI::calculate_dSdE_correction(ci_basis, *wf.Sigma());
  REQUIRE(!corr_CP.empty());

  // Same tables via the second-order qk route:
  Coulomb::QkTable qk;
  const auto yk = Coulomb::YkTable(wf.basis());
  qk.fill(wf.basis(), yk, 8, false);
  const auto [core_s1, excited_s1] =
    DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel(), 1);
  const auto corr_G =
    CI::calculate_dSdE_correction(ci_basis, core_s1, excited_s1, qk);

  // Same physics (both 2nd-order Goldstone), different machinery (sub-grid
  // matrix application vs direct basis sums): agree to sub-grid accuracy
  std::cout << "     S1(CP)      S1(qk)      dS1(CP)     dS1(qk)\n";
  for (const auto &v : ci_basis) {
    const auto S_cp = corr_CP.S1.getv(v, v);
    const auto S_g = corr_G.S1.getv(v, v);
    const auto dS_cp = corr_CP.dS1.getv(v, v);
    const auto dS_g = corr_G.dS1.getv(v, v);
    fmt::print("{:4s} {:+.8f} {:+.8f} {:+.7f} {:+.7f}\n", v.shortSymbol(), S_cp,
               S_g, dS_cp, dS_g);
    REQUIRE(S_cp == Approx(S_g).epsilon(0.05).margin(1.0e-6));
    REQUIRE(dS_cp == Approx(dS_g).epsilon(0.05).margin(1.0e-5));
  }

  // e_sigma: energy Sigma was formed at (HF valence) vs first ci-basis state
  // (basis approximates HF): close, not identical
  for (const auto &[kappa, e_cp] : corr_CP.e_sigma) {
    REQUIRE(e_cp == Approx(corr_G.e_sigma.at(kappa)).epsilon(1.0e-3));
  }

  //-----------------------------------------------------------------------
  // Round-trip: dSigma/dE read back from the sigma file

  const MBPT::CorrelationPotential Sig2(sig_fname, wf.vHF(), wf.basis(), 1.0e-4,
                                        30.0, 2, 1,
                                        MBPT::SigmaMethod::Goldstone);
  REQUIRE(Sig2.has_derivative());
  for (const auto &v : ci_basis) {
    REQUIRE(v * Sig2.dSigmaFv(v) ==
            Approx(v * wf.Sigma()->dSigmaFv(v)).margin(1.0e-12));
  }

  //-----------------------------------------------------------------------
  // fk stored alongside Sigma: storage, average_fk, and file round-trip

  {
    const auto fk_test = std::vector{0.7, 0.8, 0.9};
    const auto fk_fname = "deleteme_" + qip::random_string(3) + ".s2.abf";
    MBPT::CorrelationPotential Sig3(fk_fname, wf.vHF(), wf.basis(), 1.0e-4,
                                    30.0, 4, 1, MBPT::SigmaMethod::Goldstone,
                                    false, false, 0, {}, false, fk_test);
    Sig3.formSigma(-1, -0.19, 3);
    Sig3.formSigma(1, -0.11, 3);
    REQUIRE(Sig3.get(-1) != nullptr);
    REQUIRE(Sig3.get(-1)->fk == fk_test);
    // both stored sets identical: average is exact
    REQUIRE(Sig3.average_fk() == fk_test);
    // l_max = 0: only the s set
    REQUIRE(Sig3.average_fk(0) == fk_test);
    Sig3.write(fk_fname);

    const MBPT::CorrelationPotential Sig4(fk_fname, wf.vHF(), wf.basis(),
                                          1.0e-4, 30.0, 4, 1,
                                          MBPT::SigmaMethod::Goldstone);
    REQUIRE(Sig4.get(-1) != nullptr);
    REQUIRE(Sig4.get(-1)->fk == fk_test);
    REQUIRE(Sig4.average_fk() == fk_test);
  }

  // No stored fk: average is empty
  REQUIRE(wf.Sigma()->average_fk().empty());
}

//==============================================================================
// Three-electron CI Hamiltonian (CSF3), against independent results:
//  1. Closed s1/2 pair plus one electron: the block is the energy of the pair,
//     plus h1, plus the Hartree-Fock potential of the pair (direct and
//     exchange), which is computed independently by the HF routines.
//  2. A third electron that does not interact with the other two: the
//     spectrum is the two-electron CI spectrum plus its orbital energy. Covers
//     Coulomb, Breit, Sigma_2 (incl. extrapolation and the Brillouin-Wigner
//     pair-energy shift), and the dSigma/dE correction.
//  3. Centre-of-gravity sum rule: sum_J [J] Tr(H_J) over all CSFs of one
//     configuration equals the trace over its determinants, which is formed
//     from the two-electron (pair) traces.
TEST_CASE("CI: three-electron Hamiltonian", "[CI][3VCI][unit]") {

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"Li", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[He]", std::nullopt, 1.0e-5, false);
  auto b_params = SplineBasis::Parameters("5spd", 20, 6, {1.0e-2}, 20.0);
  b_params.verbose = false;
  wf.formBasis(b_params);

  // Closed 1s pair, and excited basis states
  const auto &a = wf.core().front();
  REQUIRE(a.n() == 1);
  std::vector<DiracSpinor> orbitals{a};
  for (const auto &Fn : wf.basis()) {
    if (Fn.n() > 1 && Fn.l() <= 2) {
      orbitals.push_back(Fn);
    }
  }

  Coulomb::QkTable qk;
  qk.fill(orbitals, Coulomb::YkTable(orbitals), -1, false);

  // One-body table: orbital energies, plus (arbitrary) off-diagonal elements
  Coulomb::meTable<double> h1;
  for (const auto &v : orbitals) {
    for (const auto &w : orbitals) {
      if (v.kappa() == w.kappa()) {
        h1.add(v, w, v == w ? v.en() : 0.01 / (v.n() + w.n()));
      }
    }
  }

  //----------------------------------------------------------------------------
  // 1. Closed pair a^2 plus one electron z: J = j_z, one CSF per z
  double worst_closed = 0.0;
  {
    const auto ia = a.nk_index();
    const auto E_pair =
      2.0 * h1.getv(ia, ia) + CI::CSF2_Coulomb(qk, ia, ia, ia, ia, 0);

    // Direct potential of the closed pair
    auto v_direct = Coulomb::yk_ab(0, a, a);
    for (auto &v : v_direct) {
      v *= a.twojp1();
    }

    for (const auto &z : orbitals) {
      if (z == a)
        continue;
      const auto csf_z = CI::form_CSF3s({ia, ia, z.nk_index()}, z.twoj());
      REQUIRE(csf_z.size() == 1);
      for (const auto &zp : orbitals) {
        if (zp == a || zp.kappa() != z.kappa())
          continue;
        const auto csf_zp = CI::form_CSF3s({ia, ia, zp.nk_index()}, zp.twoj());
        const auto H = CI::Hab(csf_z[0], csf_zp[0], z.twoj(), h1, qk);
        const auto V_HF =
          z * (v_direct * zp) + z * HF::vexFa(zp, std::vector{a});
        const auto expected = (z == zp ? E_pair : 0.0) + h1.getv(z, zp) + V_HF;
        worst_closed = std::max(worst_closed, std::abs(H - expected));
      }
    }
  }

  //----------------------------------------------------------------------------
  // 2. Third electron z (3d3/2) does not interact with the pair space S:
  // no integrals involving z, and no one-body coupling to S
  double worst_spectator = 0.0;
  {
    std::vector<DiracSpinor> S;
    std::vector<DiracSpinor> S_small; // Sk tabulated here only
    for (const auto &Fn : wf.basis()) {
      if (Fn.n() > 1 && Fn.l() <= 1) {
        S.push_back(Fn);
        if (Fn.n() <= 3) {
          S_small.push_back(Fn);
        }
      }
    }
    const auto &z =
      *std::find_if(wf.basis().cbegin(), wf.basis().cend(), [](const auto &Fn) {
        return Fn.n() == 3 && Fn.kappa() == 2;
      });

    Coulomb::QkTable qk_S;
    qk_S.fill(S, Coulomb::YkTable(S), -1, false);

    // Breit and Sigma_2 tables: arbitrary multiples of Q^k
    const auto all_k = [](int, const DiracSpinor &, const DiracSpinor &,
                          const DiracSpinor &,
                          const DiracSpinor &) { return true; };
    Coulomb::WkTable Bk;
    Bk.fill(
      S,
      [&qk_S](int k, const DiracSpinor &v, const DiracSpinor &w,
              const DiracSpinor &x,
              const DiracSpinor &y) { return 0.3 * qk_S.Q(k, v, w, x, y); },
      all_k, -1, false);
    Coulomb::WkTable Sk, dSk;
    Sk.fill(
      S_small,
      [&qk_S](int k, const DiracSpinor &v, const DiracSpinor &w,
              const DiracSpinor &x,
              const DiracSpinor &y) { return 0.2 * qk_S.Q(k, v, w, x, y); },
      all_k, -1, false);
    dSk.fill(
      S_small,
      [&qk_S](int k, const DiracSpinor &v, const DiracSpinor &w,
              const DiracSpinor &x,
              const DiracSpinor &y) { return 0.05 * qk_S.Q(k, v, w, x, y); },
      all_k, -1, false);
    // Sigma_2 extrapolation ratios, beyond S_small
    const std::vector<double> hk{0.1, 0.12, 0.14, 0.16};

    // dSigma/dE correction: no entries for z
    CI::Sigma1Correction s1c_2e;
    for (const auto &v : S) {
      s1c_2e.en[v.nk_index()] = v.en();
      s1c_2e.e_sigma.insert({v.kappa(), v.en()});
      for (const auto &w : S) {
        if (v.kappa() == w.kappa()) {
          s1c_2e.S1.add(v, w, -0.02 / (v.n() * w.n()));
          s1c_2e.dS1.add(v, w, -0.005 / (v.n() + w.n()));
        }
      }
    }
    s1c_2e.en[z.nk_index()] = z.en();
    s1c_2e.E0 = -0.6;
    const auto dE0_2e = 0.03;

    // Three-electron: the total energy includes the third electron
    auto s1c_3e = s1c_2e;
    s1c_3e.E0 = s1c_2e.E0 + z.en();
    const auto dE0_3e = dE0_2e + z.en();

    for (const int parity : {1, -1}) {
      for (const int twoJ : {1, 3, 5}) {

        std::vector<CI::CSF> csfs_3e;
        for (std::size_t i = 0; i < S.size(); ++i) {
          for (std::size_t j = i; j < S.size(); ++j) {
            for (auto &csf : CI::form_CSF3s(
                   {S[i].nk_index(), S[j].nk_index(), z.nk_index()}, twoJ)) {
              if (csf.parity() == parity) {
                csfs_3e.push_back(std::move(csf));
              }
            }
          }
        }
        const auto H_3e = CI::construct_Hci(csfs_3e, twoJ, h1, qk_S, &Bk, &Sk,
                                            &s1c_3e, hk, &dSk, dE0_3e);
        const auto E_3e = LinAlg::symmhEigensystem(H_3e).first;

        // Two-electron spectra of each J2 that couples with j_z to J
        std::vector<double> expected;
        for (int twoJ2 = std::abs(twoJ - z.twoj()); twoJ2 <= twoJ + z.twoj();
             twoJ2 += 2) {
          const CI::PsiJPi psi_2e(twoJ2, parity * z.parity(), S);
          if (psi_2e.CSFs().empty())
            continue;
          const auto H_2e = CI::construct_Hci(psi_2e, h1, qk_S, &Bk, &Sk,
                                              &s1c_2e, hk, &dSk, dE0_2e);
          for (const auto e : LinAlg::symmhEigensystem(H_2e).first) {
            expected.push_back(e + z.en());
          }
        }
        std::sort(expected.begin(), expected.end());

        REQUIRE(E_3e.size() == expected.size());
        for (std::size_t i = 0; i < expected.size(); ++i) {
          worst_spectator =
            std::max(worst_spectator, std::abs(E_3e[i] - expected[i]));
        }
      }
    }
  }

  //----------------------------------------------------------------------------
  // 3. Centre-of-gravity sum rule, for each type of configuration
  double worst_trace = 0.0;
  {
    const auto index = [&wf](int n, int kappa) {
      return std::find_if(wf.basis().cbegin(), wf.basis().cend(),
                          [n, kappa](const auto &Fn) {
                            return Fn.n() == n && Fn.kappa() == kappa;
                          })
        ->nk_index();
    };

    // Trace of g over the determinants of orbitals {x, y}
    const auto pair_trace = [&qk](DiracSpinor::Index x, DiracSpinor::Index y) {
      const auto tjx = Angular::nkindex_to_twoj(x);
      const auto tjy = Angular::nkindex_to_twoj(y);
      double trace = 0.0;
      for (int twoJ2 = std::abs(tjx - tjy); twoJ2 <= tjx + tjy; twoJ2 += 2) {
        if (x == y && twoJ2 % 4 != 0)
          continue;
        trace += (twoJ2 + 1) * CI::CSF2_Coulomb(qk, x, y, x, y, twoJ2);
      }
      return trace;
    };

    const std::vector<std::array<DiracSpinor::Index, 3>> configs{
      {index(2, -1), index(2, 1), index(3, 2)},
      {index(2, 1), index(3, 1), index(4, 1)},
      {index(2, -2), index(2, -2), index(3, -1)},
      {index(2, -1), index(2, -2), index(2, -2)},
      {index(2, -2), index(2, -2), index(2, -2)},
      {index(3, -3), index(3, -3), index(3, -3)},
      {index(3, -3), index(3, -3), index(4, -3)}};

    for (auto config : configs) {
      std::sort(config.begin(), config.end());
      const auto [x, y, z] = config;
      const auto nx = Angular::nkindex_to_twoj(x) + 1;
      const auto ny = Angular::nkindex_to_twoj(y) + 1;
      const auto nz = Angular::nkindex_to_twoj(z) + 1;

      // Number of determinants, and two-body trace over them
      double num_det = 0.0;
      double trace_g = 0.0;
      if (x != y && y != z) {
        num_det = nx * ny * nz;
        trace_g =
          nz * pair_trace(x, y) + ny * pair_trace(x, z) + nx * pair_trace(y, z);
      } else if (x == y && y != z) {
        num_det = nx * (nx - 1) / 2.0 * nz;
        trace_g = nz * pair_trace(x, x) + (nx - 1) * pair_trace(x, z);
      } else if (x != y && y == z) {
        num_det = nx * ny * (ny - 1) / 2.0;
        trace_g = nx * pair_trace(y, y) + (ny - 1) * pair_trace(x, y);
      } else {
        num_det = nx * (nx - 1) * (nx - 2) / 6.0;
        trace_g = (nx - 2) * pair_trace(x, x);
      }
      const auto trace_h1 =
        num_det * (h1.getv(x, x) + h1.getv(y, y) + h1.getv(z, z));
      const auto expected = trace_h1 + trace_g;

      double trace = 0.0;
      for (int twoJ = 1; twoJ <= nx + ny + nz - 3; twoJ += 2) {
        for (const auto &csf : CI::form_CSF3s(config, twoJ)) {
          trace += (twoJ + 1) * CI::Hab(csf, csf, twoJ, h1, qk);
        }
      }
      worst_trace =
        std::max(worst_trace, std::abs(trace - expected) / std::abs(expected));
    }
  }

  REQUIRE(worst_closed < 1.0e-10);
  REQUIRE(worst_spectator < 1.0e-10);
  REQUIRE(worst_trace < 1.0e-12);
}

//==============================================================================
// Sparse CI matrix against the dense one, with all corrections (Breit,
// Sigma_2 with extrapolation and BW shift, dSigma/dE), for two and three
// electrons: every element equal, and every non-zero element stored. The
// number stored is bounded by the pairs of CSFs that can interact, which
// equals a direct count over all pairs
TEST_CASE("CI: sparse Hamiltonian", "[CI][3VCI][unit]") {

  const auto t = UnitTest::Li_CI_tables();

  double worst_element = 0.0;
  // {2J, parity, number of electrons}
  for (const auto &[twoJ, parity, num_electrons] :
       std::vector<std::array<int, 3>>{{1, 1, 3}, {3, -1, 3}, {2, 1, 2}}) {
    const CI::PsiJPi psi{twoJ, parity, t.orbitals, num_electrons};
    const auto H_dense = CI::construct_Hci(psi, t.h1, t.qk, &t.Bk, &t.Sk,
                                           &t.s1c, t.hk, &t.dSk, t.dE0);
    const auto H_sparse =
      CI::construct_Hci<LinAlg::SparseSymmetricMatrix<double>>(
        psi, t.h1, t.qk, &t.Bk, &t.Sk, &t.s1c, t.hk, &t.dSk, t.dE0);
    const auto H_sparse_dense = H_sparse.dense();

    const auto N = psi.num_CSFs();
    std::size_t num_nonzero = 0;
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < N; ++j) {
        worst_element = std::max(
          worst_element, std::abs(H_sparse_dense(i, j) - H_dense(i, j)));
        if (j <= i && H_dense(i, j) != 0.0) {
          ++num_nonzero;
        }
      }
    }

    // Pairs that can interact, counted directly
    std::size_t num_pairs = 0;
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        if (CI::num_different(psi.csf(i), psi.csf(j)) <= 2) {
          ++num_pairs;
        }
      }
    }
    const auto num_interacting = CI::num_interacting_pairs(psi.CSFs());

    REQUIRE(H_sparse.num_stored() == num_nonzero);
    REQUIRE(num_interacting == num_pairs);
    REQUIRE(H_sparse.num_stored() <= num_interacting);
    if (num_electrons == 3) {
      // Three electrons: not every pair interacts
      REQUIRE(num_interacting < N * (N + 1) / 2);
    } else {
      REQUIRE(num_interacting == N * (N + 1) / 2);
    }
  }
  REQUIRE(worst_element == 0.0);
}

//==============================================================================
// Lowest five solutions of a three-electron CI matrix (with all
// corrections): Davidson, on the sparse matrix, vs direct diagonalisation
TEST_CASE("CI: Davidson vs dense", "[CI][3VCI]") {

  const auto t = UnitTest::Li_CI_tables();

  double worst_E = 0.0, worst_overlap = 0.0;
  for (const auto &[twoJ, parity] :
       std::vector<std::pair<int, int>>{{1, 1}, {3, -1}}) {
    CI::PsiJPi psi_dense{twoJ, parity, t.orbitals, 3};
    auto psi_iterative = psi_dense;
    psi_dense.solve(CI::construct_Hci(psi_dense, t.h1, t.qk, &t.Bk, &t.Sk,
                                      &t.s1c, t.hk, &t.dSk, t.dE0),
                    5);
    const auto convergence = psi_iterative.solve(
      CI::construct_Hci<LinAlg::SparseSymmetricMatrix<double>>(
        psi_iterative, t.h1, t.qk, &t.Bk, &t.Sk, &t.s1c, t.hk, &t.dSk, t.dE0),
      5);
    REQUIRE(convergence.converged);
    REQUIRE(convergence.iterations > 0);
    REQUIRE(psi_iterative.num_solutions() == 5);
    for (std::size_t i = 0; i < 5; ++i) {
      worst_E = std::max(
        worst_E, std::abs(psi_iterative.energy(i) - psi_dense.energy(i)));
      double overlap = 0.0;
      for (std::size_t j = 0; j < psi_dense.num_CSFs(); ++j) {
        overlap += psi_iterative.coef(i, j) * psi_dense.coef(i, j);
      }
      worst_overlap = std::max(worst_overlap, 1.0 - std::abs(overlap));
    }
  }
  fmt::print("E vs dense:          max diff = {:.1e}\n", worst_E);
  fmt::print("|<c|c_dense>| vs 1:  max diff = {:.1e}\n", worst_overlap);
  REQUIRE(worst_E < 1.0e-12);
  REQUIRE(worst_overlap < 1.0e-12);

  // All solutions below an energy: the number found is not known in
  // advance (here 20: found 8, 16, then 32, of which 20 are below)
  {
    CI::PsiJPi psi_dense{1, 1, t.orbitals, 3};
    auto psi_iterative = psi_dense;
    const auto H_dense = CI::construct_Hci(psi_dense, t.h1, t.qk, &t.Bk, &t.Sk,
                                           &t.s1c, t.hk, &t.dSk, t.dE0);
    const auto H_sparse =
      CI::construct_Hci<LinAlg::SparseSymmetricMatrix<double>>(
        psi_iterative, t.h1, t.qk, &t.Bk, &t.Sk, &t.s1c, t.hk, &t.dSk, t.dE0);
    const auto E_all = LinAlg::symmhEigensystem(H_dense).first;
    const auto all_below_cm =
      0.5 * (E_all[19] + E_all[20]) * PhysConst::Hartree_invcm;
    psi_dense.solve(H_dense, 0, all_below_cm);
    const auto convergence = psi_iterative.solve(H_sparse, 0, all_below_cm);
    REQUIRE(convergence.converged);
    // Found iteratively
    REQUIRE(convergence.iterations > 0);
    REQUIRE(psi_dense.num_solutions() == 20);
    REQUIRE(psi_iterative.num_solutions() == 20);
    double worst_E_below = 0.0;
    for (std::size_t i = 0; i < 20; ++i) {
      worst_E_below = std::max(
        worst_E_below, std::abs(psi_iterative.energy(i) - psi_dense.energy(i)));
    }
    fmt::print("All below, E vs dense: max diff = {:.1e}\n", worst_E_below);
    REQUIRE(worst_E_below < 1.0e-12);

    // Threshold above every solution: as many as Davidson finds, N - 1
    const auto above_all =
      (E_all[E_all.size() - 1] + 1.0) * PhysConst::Hartree_invcm;
    std::cout << "Expect two warnings below (deliberate): many solutions\n";
    REQUIRE(psi_iterative.solve(H_sparse, 0, above_all).converged);
    REQUIRE(psi_iterative.num_solutions() == psi_iterative.num_CSFs() - 1);
    REQUIRE(psi_iterative.energy(5) == Approx(E_all[5]).epsilon(1.0e-12));

    // More than N/20 solutions requested (slower than dense, but works)
    const auto many = int(psi_iterative.num_CSFs() / 20) + 1;
    const auto convergence_many = psi_iterative.solve(H_sparse, many);
    REQUIRE(convergence_many.converged);
    REQUIRE(psi_iterative.num_solutions() == std::size_t(many));
    REQUIRE(psi_iterative.energy(std::size_t(many) - 1) ==
            Approx(E_all[std::size_t(many) - 1]).epsilon(1.0e-12));
  }

  // A block of a single CSF (1s^2 2p-, J=1/2-, from the orbitals 1s and
  // 2p-), as the first pass of iterate_E0 requests: the one solution, without
  // an eigensolve
  {
    const auto p_minus =
      std::find_if(t.orbitals.cbegin(), t.orbitals.cend(), [](const auto &Fn) {
        return Fn.n() == 2 && Fn.kappa() == 1;
      });
    REQUIRE(p_minus != t.orbitals.cend());
    const std::vector<DiracSpinor> two_orbitals{t.orbitals.front(), *p_minus};
    CI::PsiJPi psi_single{1, -1, two_orbitals, 3};
    REQUIRE(psi_single.num_CSFs() == 1);
    const auto H_single =
      CI::construct_Hci<LinAlg::SparseSymmetricMatrix<double>>(psi_single, t.h1,
                                                               t.qk);
    const auto convergence_single = psi_single.solve(H_single, 1);
    REQUIRE(convergence_single.iterations == 0);
    REQUIRE(psi_single.num_solutions() == 1);
    REQUIRE(psi_single.energy(0) == H_single.diagonal()[0]);
    REQUIRE(std::abs(psi_single.coef(0, 0)) == 1.0);
  }
}

//==============================================================================
// Self-consistent dSigma/dE reference energy E0, and the CI matrix built
// with it: sparse (Davidson) vs dense
TEST_CASE("CI: iterate_E0 sparse", "[CI][3VCI]") {

  const auto t = UnitTest::Li_CI_tables();

  std::ostringstream discard;
  auto s1c_unset = t.s1c;
  s1c_unset.E0 = 0.0;
  CI::PsiJPi psi_dense{1, 1, t.orbitals, 3};
  auto psi_iterative = psi_dense;
  const auto H_dense = CI::iterate_E0(&psi_dense, s1c_unset, t.h1, t.qk, &t.Bk,
                                      &t.Sk, t.hk, &t.dSk, -7.4, discard);
  const auto H_sparse = CI::iterate_E0<LinAlg::SparseSymmetricMatrix<double>>(
                          &psi_iterative, s1c_unset, t.h1, t.qk, &t.Bk, &t.Sk,
                          t.hk, &t.dSk, -7.4, discard)
                          .dense();
  double worst_H = 0.0;
  for (std::size_t i = 0; i < psi_dense.num_CSFs(); ++i) {
    for (std::size_t j = 0; j < psi_dense.num_CSFs(); ++j) {
      worst_H = std::max(worst_H, std::abs(H_sparse(i, j) - H_dense(i, j)));
    }
  }
  const auto dE0 = std::abs(psi_iterative.energy(0) - psi_dense.energy(0));
  fmt::print("E0 vs dense:  diff     = {:.1e}\n", dE0);
  fmt::print("H vs dense:   max diff = {:.1e}\n", worst_H);
  REQUIRE(dE0 < 1.0e-12);
  REQUIRE(worst_H < 1.0e-12);
}

//==============================================================================
// E0 for each level: the lowest level is unchanged, and each excited level is
// much closer to the eigenvalue of the Hamiltonian built at its own energy
// than with the lowest level's E0. Dense and Davidson agree
TEST_CASE("CI: E0 for each level", "[CI]") {

  const auto t = UnitTest::Li_CI_tables();
  const double E0_sigma2 = -7.4;
  const int num_solutions = 3;

  std::ostringstream discard;
  const auto run = [&](bool BW_each_level, CI::Solver solver) {
    return CI::run_CI(t.orbitals, 1, +1, 3, num_solutions, {}, t.h1, t.qk, t.Bk,
                      t.Sk, true, false, false, discard, "", &t.s1c, t.hk,
                      &t.dSk, E0_sigma2, solver, BW_each_level);
  };
  const auto psi_lowest = run(false, CI::Solver::dense);
  const auto psi_each = run(true, CI::Solver::dense);
  const auto psi_each_davidson = run(true, CI::Solver::Davidson);
  REQUIRE(psi_each.num_solutions() == std::size_t(num_solutions));

  // Energy of level i of the Hamiltonian built with reference energy E
  const auto level_energy = [&t, E0_sigma2](std::size_t i, double E) {
    auto s1c_E = t.s1c;
    s1c_E.E0 = E;
    CI::PsiJPi psi{1, +1, t.orbitals, 3};
    const auto Hci = CI::construct_Hci(psi, t.h1, t.qk, &t.Bk, &t.Sk, &s1c_E,
                                       t.hk, &t.dSk, E - E0_sigma2);
    psi.solve(Hci, int(i) + 1);
    return psi.energy(i);
  };

  const auto dE_lowest = std::abs(psi_each.energy(0) - psi_lowest.energy(0));
  double worst_each = 0.0;
  double worst_lowest = 0.0;
  double worst_davidson = 0.0;
  for (std::size_t i = 0; i < psi_each.num_solutions(); ++i) {
    worst_davidson =
      std::max(worst_davidson,
               std::abs(psi_each_davidson.energy(i) - psi_each.energy(i)));
    if (i == 0)
      continue;
    const auto E_each = psi_each.energy(i);
    const auto E_lowest = psi_lowest.energy(i);
    const auto residual_each = std::abs(E_each - level_energy(i, E_each));
    const auto residual_lowest = std::abs(E_lowest - level_energy(i, E_lowest));
    REQUIRE(residual_each < residual_lowest);
    worst_each = std::max(worst_each, residual_each);
    worst_lowest = std::max(worst_lowest, residual_lowest);
  }
  fmt::print("Level 0, each vs lowest E0:          diff     = {:.1e}\n",
             dE_lowest);
  fmt::print("E_i vs level i of H(E_i), each:      max diff = {:.1e}\n",
             worst_each);
  fmt::print("E_i vs level i of H(E_i), lowest E0: max diff = {:.1e}\n",
             worst_lowest);
  fmt::print("Each level, Davidson vs dense:       max diff = {:.1e}\n",
             worst_davidson);
  REQUIRE(dE_lowest < 1.0e-12);
  REQUIRE(worst_davidson < 1.0e-9);
}

//==============================================================================
// Three-electron CI, run through configuration_interaction: input parsing
// (fractions, twoJ), the odd-2J blocks, the ci file, and the plumbing
// (block energies reproduce the directly constructed Hamiltonian)
TEST_CASE("CI: three-electron run", "[CI][3VCI][unit]") {

  // Li, all three electrons in the CI space (bare-nucleus basis)
  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"Li", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[]", std::nullopt, 1.0e-5, false);
  auto b_params = SplineBasis::Parameters("3sp", 20, 6, {1.0e-2}, 20.0);
  b_params.verbose = false;
  wf.formBasis(b_params);

  const std::string qk_filename =
    "deleteme_" + qip::random_string(3) + ".qk.abf";

  // Two configuration spaces: every configuration (three full steps), and
  // singles and doubles from given leading configurations (over 1s, 2p), so
  // the test does not depend on the default
  const std::vector<std::string> sd_leading{"1s2.2p", "1s.2p2", "2p3"};
  const std::vector<std::string> sd_excitations{"1", "2"};
  for (const bool full_space : {true, false}) {
    const std::string space_input =
      full_space ? "excitations = 1, 2, 3;" :
                   "leading = 1s2.2p, 1s.2p2, 2p3; excitations = 1, 2;";
    INFO(space_input);

    const std::string ci_filename =
      "deleteme_" + qip::random_string(3) + ".ci.abf";
    const auto files =
      "qk_file = " + qk_filename + "; ci_file = " + ci_filename +
      "; num_solutions = 3; print_details = false;" + space_input;

    const IO::InputBlock input{"CI", "ci_basis = 3sp; num_electrons = 3;"
                                     "J = 1/2, 3/2;" +
                                       files};
    const auto solution = CI::configuration_interaction(input, wf);
    const auto &levels = solution.levels;
    REQUIRE(solution.integrals.num_electrons == 3);
    REQUIRE(levels.size() == 4);

    const auto &basis = solution.integrals.ci_basis;
    const auto space =
      full_space ? CI::ConfigurationSpace{} :
                   CI::ConfigurationSpace{basis, 3, sd_leading, sd_excitations};
    for (const auto &psi : levels) {
      REQUIRE(psi.num_electrons() == 3);
      REQUIRE(psi.twoJ() % 2 == 1);
      // CSFs of exactly the configurations in the space; singles and
      // doubles leave out some of every block
      const auto num_CSFs_all =
        CI::form_CSF3s(psi.twoJ(), psi.parity(), basis).size();
      const auto num_CSFs_space =
        CI::form_CSF3s(psi.twoJ(), psi.parity(), basis, space).size();
      REQUIRE(psi.num_CSFs() == num_CSFs_space);
      if (!full_space) {
        REQUIRE(num_CSFs_space < num_CSFs_all);
      }
      REQUIRE(psi.num_solutions() == 3);
      // Energies reproduce the directly constructed Hamiltonian
      const auto Hci = CI::construct_Hci(psi, solution.integrals);
      const auto E = LinAlg::symmhEigensystem(Hci).first;
      for (std::size_t i = 0; i < psi.num_solutions(); ++i) {
        REQUIRE(psi.energy(i) == Approx(E[i]).epsilon(1.0e-10));
      }
      // Every CSF has three electrons, and the config strings reflect that
      REQUIRE(psi.csf(0).num_electrons() == 3);
      REQUIRE(psi.csf(0).parity() == psi.parity());
    }

    // Ground state: 1s^2 ns, J=1/2 even (in the bare-nucleus basis, the 2s
    // and 3s spline orbitals carry similar weight). Bounds: above the exact
    // energy (-7.478 au, non-rel), and below the single-configuration limit
    const auto *gs = &levels.front();
    REQUIRE(gs->twoJ() == 1);
    REQUIRE(gs->parity() == 1);
    REQUIRE(gs->info(0).config.rfind("1s^2 ", 0) == 0);
    REQUIRE(gs->info(0).config.back() == 's');
    REQUIRE(gs->energy(0) > -7.48);
    REQUIRE(gs->energy(0) < -7.3);

    // Term identification: 2S1/2 ground state (gJ = 2), and 1s^2 np 2P1/2 and
    // 2P3/2 (gJ = 2/3, 4/3, in the non-rel limit)
    REQUIRE(gs->info(0).gJ == Approx(2.0).epsilon(1.0e-2));
    REQUIRE(gs->info(0).L2 == Approx(0.0).margin(0.05));
    REQUIRE(gs->info(0).S2 == Approx(0.75).margin(0.05));
    REQUIRE(gs->info(0).twoS == 1.0);
    REQUIRE(gs->info(0).L == 0.0);
    for (const auto &psi : levels) {
      if (psi.parity() != -1)
        continue;
      const auto gJ_P = psi.twoJ() == 1 ? 2.0 / 3.0 : 4.0 / 3.0;
      REQUIRE(psi.info(0).gJ == Approx(gJ_P).epsilon(1.0e-2));
      REQUIRE(psi.info(0).L2 == Approx(2.0).margin(0.05));
      REQUIRE(psi.info(0).S2 == Approx(0.75).margin(0.05));
      REQUIRE(psi.info(0).twoS == 1.0);
      REQUIRE(psi.info(0).L == 1.0);
    }

    // Same blocks via twoJ (which supersedes J), read from the ci file
    // (read_only: a failed read would leave the blocks unsolved)
    const IO::InputBlock input_twoJ{"CI",
                                    "ci_basis = 3sp; num_electrons = 3;"
                                    "J = 2.5; twoJ = 1, 3; read_only = true;" +
                                      files};
    const auto levels_twoJ =
      CI::configuration_interaction(input_twoJ, wf).levels;
    REQUIRE(levels_twoJ.size() == levels.size());
    for (std::size_t i = 0; i < levels.size(); ++i) {
      REQUIRE(levels_twoJ[i].twoJ() == levels[i].twoJ());
      REQUIRE(levels_twoJ[i].parity() == levels[i].parity());
      for (std::size_t s = 0; s < levels[i].num_solutions(); ++s) {
        REQUIRE(levels_twoJ[i].energy(s) == Approx(levels[i].energy(s)));
        for (std::size_t k = 0; k < levels[i].num_CSFs(); ++k) {
          REQUIRE(levels_twoJ[i].coef(s, k) == Approx(levels[i].coef(s, k)));
        }
      }
    }

    // Integer J is not possible for three electrons: no blocks
    std::cout << "\nExpect one warning below (deliberate): J = 1 is not "
                 "possible for 3 electrons\n";
    const IO::InputBlock input_bad{"CI", "ci_basis = 3sp; num_electrons = 3;"
                                         "J = 1;" +
                                           files};
    REQUIRE(CI::configuration_interaction(input_bad, wf).levels.empty());

    // Parity-specific list supersedes the plain one, for its parity only:
    // odd parity takes twoJ- (rejected), even parity takes J
    std::cout << "\nExpect one warning below (deliberate): twoJ- = 2\n";
    const IO::InputBlock input_pi{"CI",
                                  "ci_basis = 3sp; num_electrons = 3;"
                                  "J = 1/2; twoJ- = 2; read_only = true;" +
                                    files};
    const auto levels_pi = CI::configuration_interaction(input_pi, wf).levels;
    REQUIRE(levels_pi.size() == 1);
    REQUIRE(levels_pi[0].twoJ() == 1);
    REQUIRE(levels_pi[0].parity() == 1);
    REQUIRE(levels_pi[0].num_solutions() == 3);
  }
}

//==============================================================================
// One valence electron: the CI Hamiltonian of each (J, parity) is the one-body
// h1 between the orbitals of that kappa. Without Sigma_1 it is diagonal, so
// each level is one orbital: its energy, g-factor, and term. Matrix elements
// between levels are the single-particle ones, and the mixed states are the
// sums over the orbitals. With Sigma_1, h1 has off-diagonal elements, which
// the Hamiltonian must carry
TEST_CASE("CI: one-electron run", "[CI][unit]") {

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"Li", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[He]", std::nullopt, 1.0e-5, false);
  auto b_params = SplineBasis::Parameters("4spd", 20, 6, {1.0e-2}, 20.0);
  b_params.verbose = false;
  wf.formBasis(b_params);

  const std::string qk_filename =
    "deleteme_" + qip::random_string(3) + ".qk.abf";
  const std::string ci_filename =
    "deleteme_" + qip::random_string(3) + ".ci.abf";
  const auto files = "qk_file = " + qk_filename + "; ci_file = " + ci_filename +
                     "; print_details = false;";

  const IO::InputBlock input{"CI", "ci_basis = 4spd; num_electrons = 1;"
                                   "J = 1/2, 3/2; num_solutions = 2;" +
                                     files};
  const auto solution = CI::configuration_interaction(input, wf);
  const auto &levels = solution.levels;
  const auto &basis = solution.integrals.ci_basis;
  REQUIRE(solution.integrals.num_electrons == 1);
  // 1/2+ (s), 3/2+ (d3/2), 1/2- (p1/2), 3/2- (p3/2)
  REQUIRE(levels.size() == 4);

  // The CI orbitals of one kappa, in energy order
  const auto orbitals_of = [&basis](int twoJ, int parity) {
    std::vector<DiracSpinor> orbitals;
    for (const auto &Fn : basis) {
      if (Fn.twoj() == twoJ && Fn.parity() == parity) {
        orbitals.push_back(Fn);
      }
    }
    std::sort(orbitals.begin(), orbitals.end(),
              [](const auto &a, const auto &b) { return a.en() < b.en(); });
    return orbitals;
  };
  const DiracOperator::M1 m1{wf.grid(), PhysConst::alpha, 0.0};

  for (const auto &psi : levels) {
    REQUIRE(psi.num_electrons() == 1);
    const auto twoJ = psi.twoJ();
    const auto orbitals = orbitals_of(twoJ, psi.parity());
    REQUIRE(orbitals.size() >= 2);
    REQUIRE(psi.num_CSFs() == orbitals.size());
    REQUIRE(psi.num_solutions() == 2);
    for (std::size_t i = 0; i < psi.num_solutions(); ++i) {
      const auto &Fa = orbitals[i];
      // Diagonal Hamiltonian: each level is one orbital
      REQUIRE(psi.energy(i) == Approx(Fa.en()).epsilon(1.0e-10));
      REQUIRE(psi.csf(i).config(true) == Fa.shortSymbol());
      REQUIRE(std::abs(psi.coef(i, i)) == Approx(1.0).epsilon(1.0e-10));
      auto nonrel_symbol = Fa.shortSymbol();
      nonrel_symbol.pop_back();
      REQUIRE(psi.info(i).config == nonrel_symbol);
      REQUIRE(psi.info(i).ci2 == Approx(1.0).epsilon(1.0e-10));
      // Term and g-factor of the orbital: gJ = 3j * <a||M1||a> / J
      const auto l = Fa.l();
      REQUIRE(psi.info(i).L2 == Approx(l * (l + 1.0)).epsilon(1.0e-12));
      REQUIRE(psi.info(i).S2 == Approx(0.75).epsilon(1.0e-12));
      REQUIRE(psi.info(i).L == 1.0 * l);
      REQUIRE(psi.info(i).twoS == 1.0);
      const auto gJ = Angular::threej_2(twoJ, twoJ, 2, twoJ, -twoJ, 0) *
                      m1.reducedME(Fa, Fa) / (0.5 * twoJ);
      REQUIRE(psi.info(i).gJ == Approx(gJ).epsilon(1.0e-10));
    }
  }

  // Matrix elements between levels: the single-particle reduced matrix
  // elements (up to the sign of the CI eigenvectors)
  const auto level = [&levels](int twoJ, int parity) {
    return *std::find_if(levels.cbegin(), levels.cend(),
                         [twoJ, parity](const auto &psi) {
                           return psi.twoJ() == twoJ && psi.parity() == parity;
                         });
  };
  const auto &s_half = level(1, +1);
  const auto &p_half = level(1, -1);
  const auto &p_three = level(3, -1);
  const auto s_orbitals = orbitals_of(1, +1);
  const auto p_orbitals = orbitals_of(1, -1);
  const auto pp_orbitals = orbitals_of(3, -1);
  const auto &F2s = s_orbitals.front();
  const auto &F2p = p_orbitals.front();
  const auto &F2pp = pp_orbitals.front();

  const DiracOperator::E1 e1{wf.grid()};
  const auto e1_tab = Amplitudes::me_table(basis, &e1);
  const auto m1_tab = Amplitudes::me_table(basis, &m1);
  const auto d_sp = CI::ReducedME(s_half, 0, p_half, 0, e1_tab, 1, -1);
  const auto d_spp = CI::ReducedME(s_half, 0, p_three, 0, e1_tab, 1, -1);
  const auto m_ppp = CI::ReducedME(p_half, 0, p_three, 0, m1_tab, 1, +1);
  REQUIRE(d_sp != 0.0);
  REQUIRE(std::abs(d_sp) ==
          Approx(std::abs(e1.reducedME(F2s, F2p))).epsilon(1.0e-10));
  REQUIRE(std::abs(d_spp) ==
          Approx(std::abs(e1.reducedME(F2s, F2pp))).epsilon(1.0e-10));
  REQUIRE(std::abs(m_ppp) ==
          Approx(std::abs(m1.reducedME(F2p, F2pp))).epsilon(1.0e-10));
  // Parity selection: E1 within a parity, M1 across, vanish
  REQUIRE(CI::ReducedME(p_half, 0, p_three, 0, e1_tab, 1, -1) == 0.0);
  REQUIRE(CI::ReducedME(s_half, 0, p_half, 0, m1_tab, 1, +1) == 0.0);

  // Mixed state of 2s under E1, in the p1/2 block: with a diagonal
  // Hamiltonian, c_n = -<n||d||2s> c_2s / (e_n - e_2s - omega)
  const auto Hci = CI::construct_Hci(p_half, solution.integrals);
  for (const double omega : {0.0, 0.05}) {
    const auto dPsi =
      CI::solve_mixed_state(s_half, 0, p_half, Hci, e1_tab, 1, omega);
    REQUIRE(dPsi.num_electrons() == 1);
    REQUIRE(dPsi.num_CSFs() == p_half.num_CSFs());
    const auto c_2s = s_half.coef(0, 0);
    for (std::size_t n = 0; n < p_half.num_CSFs(); ++n) {
      const auto &Fn = p_orbitals[n];
      REQUIRE(p_half.csf(n).config(true) == Fn.shortSymbol());
      const auto expected =
        -e1.reducedME(Fn, F2s) * c_2s / (Fn.en() - F2s.en() - omega);
      REQUIRE(dPsi.coef(0, n) == Approx(expected).epsilon(1.0e-10));
    }
  }

  // With Sigma_1 (second order, from the Coulomb integrals), h1 has
  // off-diagonal elements between the orbitals of a kappa: the Hamiltonian is
  // that table, and the levels are its eigenvalues
  const IO::InputBlock input_s1{"CI",
                                "ci_basis = 4spd; num_electrons = 1; J = 1/2;"
                                "sigma1 = true; s1_basis = 6spd;"
                                "num_solutions = 0; ci_file = false;"
                                "qk_file = " +
                                  qk_filename + "; print_details = false;"};
  const auto solution_s1 = CI::configuration_interaction(input_s1, wf);
  REQUIRE(solution_s1.levels.size() == 2);
  for (const auto &psi : solution_s1.levels) {
    const auto N = psi.num_CSFs();
    REQUIRE(N >= 2);
    REQUIRE(psi.num_solutions() == N);
    const auto H = CI::construct_Hci(psi, solution_s1.integrals);
    double largest_off_diagonal = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < N; ++j) {
        const auto h1_ij = solution_s1.integrals.h1.getv(psi.csf(i).state(0),
                                                         psi.csf(j).state(0));
        REQUIRE(H(i, j) == h1_ij);
        if (i != j) {
          largest_off_diagonal =
            std::max(largest_off_diagonal, std::abs(h1_ij));
        }
      }
    }
    REQUIRE(largest_off_diagonal > 1.0e-6);
    const auto E = LinAlg::symmhEigensystem(H).first;
    for (std::size_t i = 0; i < N; ++i) {
      REQUIRE(psi.energy(i) == Approx(E[i]).epsilon(1.0e-10));
    }
    // Correlations lower the lowest level below its orbital energy
    const auto e_orbital = orbitals_of(psi.twoJ(), psi.parity()).front().en();
    REQUIRE(psi.energy(0) < e_orbital);
  }
}

//==============================================================================
// One valence electron, Na (core [Ne]), against single-particle calculations
// on the grid. Without MBPT, the CI levels are the Hartree-Fock valence
// energies. With Sigma_1 (second order, from the Coulomb integrals), the CI
// diagonalises h_HF + Sigma in the CI basis: the levels are the second-order
// Brueckner energies, for which h_HF + Sigma is solved on the grid (with
// Sigma from the Goldstone method, as a radial matrix). Agreement is limited
// by the completeness of the CI basis (a 12spd CI basis misses over 1% of the
// 3s correction), and by the different evaluation of Sigma (radial matrix on
// a sub-grid, vs sums over Coulomb integrals)
TEST_CASE("CI: one-electron Na", "[CI][integration]") {

  Wavefunction wf({1500, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear"},
                  {"Na", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-8, false);
  wf.solve_valence("4sp3d", false);
  auto b_params = SplineBasis::Parameters("30spd", 35, 7, {1.0e-5}, 50.0);
  b_params.verbose = false;
  wf.formBasis(b_params);

  // CI energies of the valence states: the level of the (J, parity) of each,
  // indexed by the number of lower valence states of the same kappa
  const auto ci_energies = [&wf](const std::string &options) {
    const IO::InputBlock input{
      "CI", "ci_basis = 20spd; num_electrons = 1; J+ = 1/2, 3/2, 5/2;"
            "J- = 1/2, 3/2; num_solutions = 2; qk_file = false;"
            "ci_file = false; print_details = false;" +
              options};
    const auto levels = CI::configuration_interaction(input, wf).levels;
    std::vector<double> energies;
    for (const auto &Fv : wf.valence()) {
      const auto level =
        std::find_if(levels.cbegin(), levels.cend(), [&Fv](const auto &psi) {
          return psi.twoJ() == Fv.twoj() && psi.parity() == Fv.parity();
        });
      REQUIRE(level != levels.cend());
      const auto index = std::count_if(
        wf.valence().cbegin(), wf.valence().cend(), [&Fv](const auto &Fw) {
          return Fw.kappa() == Fv.kappa() && Fw.n() < Fv.n();
        });
      energies.push_back(level->energy(std::size_t(index)));
    }
    return energies;
  };

  const auto E_CI = ci_energies("");
  // Second-order Sigma_1, calculated within CI (there is no correlation
  // potential yet)
  const auto E_CI_Sigma1 = ci_energies("sigma1 = true;");

  std::vector<double> E_HF;
  for (const auto &Fv : wf.valence()) {
    E_HF.push_back(Fv.en());
  }

  // Second-order Brueckner: Goldstone Sigma (including the small components,
  // as the Coulomb integrals do), at the lowest valence energy of each kappa
  wf.formSigma(1, 1.0e-4, 30.0, 6, false, true, false, 0, {}, {}, {}, false);
  wf.hartreeFockBrueckner(false);

  std::cout << "\n";
  fmt::print("CI (one valence electron) vs Hartree-Fock and Brueckner, Na "
             "(/cm)\n");
  fmt::print("       HF         CI     diff  Brueckner   CI+S1     diff  "
             "diff/dE_Sigma\n");
  // CI vs HF in /cm; CI+Sigma_1 vs Brueckner relative to the correlation
  // correction, dE_Sigma = E_Brueckner - E_HF
  double max_diff_HF = 0.0;
  double max_relative_diff_Br = 0.0;
  for (std::size_t i = 0; i < E_HF.size(); ++i) {
    const auto &Fv = wf.valence()[i];
    const auto E_Br = Fv.en();
    const auto diff_HF = E_CI[i] - E_HF[i];
    const auto diff_Br = E_CI_Sigma1[i] - E_Br;
    const auto relative_diff_Br = diff_Br / (E_Br - E_HF[i]);
    fmt::print(
      "{:5s} {:9.2f} {:9.2f} {:6.3f} {:9.2f} {:9.2f} {:6.3f}  {:+.1e}\n",
      Fv.shortSymbol(), E_HF[i] * PhysConst::Hartree_invcm,
      E_CI[i] * PhysConst::Hartree_invcm, diff_HF * PhysConst::Hartree_invcm,
      E_Br * PhysConst::Hartree_invcm,
      E_CI_Sigma1[i] * PhysConst::Hartree_invcm,
      diff_Br * PhysConst::Hartree_invcm, relative_diff_Br);
    max_diff_HF = std::max(max_diff_HF, std::abs(diff_HF));
    max_relative_diff_Br =
      std::max(max_relative_diff_Br, std::abs(relative_diff_Br));
  }
  const auto max_diff_HF_invcm = max_diff_HF * PhysConst::Hartree_invcm;
  fmt::print("CI vs HF: max diff = {:.3f} /cm\n", max_diff_HF_invcm);
  fmt::print("CI+Sigma_1 vs Brueckner: max diff = {:.1e} of dE_Sigma\n",
             max_relative_diff_Br);
  REQUIRE(max_diff_HF_invcm < 1.0);
  REQUIRE(max_relative_diff_Br < 0.01);
}

//==============================================================================
// The CI solutions file: blocks with different numbers of solutions; a new
// block appended; a block re-written with more solutions (the other block
// kept exactly); an incomplete final block dropped; the file size; and a file
// without the format tag (blocks holding all num_CSFs solutions), which is
// not read, and is replaced on write
TEST_CASE("CI: solutions file", "[CI][unit]") {

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"He", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[]", std::nullopt, 1.0e-5, false);
  auto b_params = SplineBasis::Parameters("4sp", 20, 6, {1.0e-2}, 20.0);
  b_params.verbose = false;
  wf.formBasis(b_params);

  // Solutions of a random symmetric matrix: distinct, non-trivial
  // coefficients
  std::mt19937_64 rng{123456789};
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);
  const auto solved = [&](int twoJ, int num_solutions) {
    CI::PsiJPi psi{twoJ, +1, wf.basis()};
    const auto N = psi.num_CSFs();
    LinAlg::Matrix<double> H(N, N);
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        H(i, j) = H(j, i) = uniform(rng);
      }
    }
    psi.solve(H, num_solutions);
    return psi;
  };
  const auto A2 = solved(0, 2);
  const auto B2 = solved(2, 2);
  const auto A3 = solved(0, 3);
  const auto N_A = A2.num_CSFs();
  const auto N_B = B2.num_CSFs();

  // Format tag, then per block: 2J, parity, N, num_solutions, E, coefficients
  const auto block_bytes = [](std::size_t N, std::size_t num_solutions) {
    return 2 * sizeof(int) + 2 * sizeof(std::size_t) +
           (num_solutions + num_solutions * N) * sizeof(double);
  };

  const auto fname = "deleteme_" + qip::random_string(4) + ".ci.abf";
  auto write_A2 = A2;
  auto write_B2 = B2;
  auto write_A3 = A3;
  REQUIRE(write_A2.read_write(fname, IO::FRW::write));
  // New block: appended
  REQUIRE(write_B2.read_write(fname, IO::FRW::write));
  REQUIRE(std::filesystem::file_size(fname) ==
          8 + block_bytes(N_A, 2) + block_bytes(N_B, 2));
  // A grows: the file is re-written
  REQUIRE(write_A3.read_write(fname, IO::FRW::write));

  CI::PsiJPi read_A{0, +1, wf.basis()};
  CI::PsiJPi read_B{2, +1, wf.basis()};
  REQUIRE(read_A.read_write(fname, IO::FRW::read));
  REQUIRE(read_B.read_write(fname, IO::FRW::read));
  REQUIRE(read_A.num_solutions() == 3);
  REQUIRE(read_B.num_solutions() == 2);
  // Binary round trip: exact
  for (std::size_t i = 0; i < 3; ++i) {
    REQUIRE(read_A.energy(i) == A3.energy(i));
    for (std::size_t j = 0; j < N_A; ++j) {
      REQUIRE(read_A.coef(i, j) == A3.coef(i, j));
    }
  }
  for (std::size_t i = 0; i < 2; ++i) {
    REQUIRE(read_B.energy(i) == B2.energy(i));
    for (std::size_t j = 0; j < N_B; ++j) {
      REQUIRE(read_B.coef(i, j) == B2.coef(i, j));
    }
  }
  REQUIRE(std::filesystem::file_size(fname) ==
          8 + block_bytes(N_A, 3) + block_bytes(N_B, 2));

  // Incomplete final block (A, written last), e.g. from an interrupted
  // write: a new block (C) is not appended after it; the file is re-written
  // without it
  std::filesystem::resize_file(fname, std::filesystem::file_size(fname) - 5);
  auto write_C2 = solved(4, 2);
  const auto N_C = write_C2.num_CSFs();
  REQUIRE(write_C2.read_write(fname, IO::FRW::write));
  REQUIRE(std::filesystem::file_size(fname) ==
          8 + block_bytes(N_B, 2) + block_bytes(N_C, 2));
  CI::PsiJPi read_B_again{2, +1, wf.basis()};
  CI::PsiJPi read_C{4, +1, wf.basis()};
  REQUIRE(read_B_again.read_write(fname, IO::FRW::read));
  REQUIRE(read_C.read_write(fname, IO::FRW::read));
  REQUIRE(read_B_again.coef(1, N_B - 1) == B2.coef(1, N_B - 1));
  REQUIRE(read_C.coef(1, N_C - 1) == write_C2.coef(1, N_C - 1));
  REQUIRE(!read_A.read_write(fname, IO::FRW::read));

  // No format tag: twoJ, pi, N, num_solutions, E[N], M[N x N]
  const auto fname_old = "deleteme_" + qip::random_string(4) + ".ci.abf";
  {
    std::fstream f;
    IO::FRW::open_binary(f, fname_old, IO::FRW::write);
    int twoJ = 0, pi = 1;
    std::size_t N = N_A, num_solutions = 2;
    IO::FRW::rw_binary(f, IO::FRW::write, twoJ, pi, N, num_solutions);
    std::vector<double> E(N, 0.0), M(N * N, 0.0);
    for (std::size_t i = 0; i < 2; ++i) {
      E[i] = A2.energy(i);
      for (std::size_t j = 0; j < N; ++j) {
        M[i * N + j] = A2.coef(i, j);
      }
    }
    f.write(reinterpret_cast<const char *>(E.data()),
            std::streamsize(E.size() * sizeof(double)));
    f.write(reinterpret_cast<const char *>(M.data()),
            std::streamsize(M.size() * sizeof(double)));
  }
  std::cout << "Expect a message below (deliberate): file without the format "
               "tag is not read\n";
  CI::PsiJPi read_old{0, +1, wf.basis()};
  REQUIRE(!read_old.read_write(fname_old, IO::FRW::read));
  REQUIRE(write_A2.read_write(fname_old, IO::FRW::write));
  REQUIRE(read_old.read_write(fname_old, IO::FRW::read));
  REQUIRE(read_old.num_solutions() == 2);
  REQUIRE(read_old.energy(1) == A2.energy(1));
  REQUIRE(std::filesystem::file_size(fname_old) == 8 + block_bytes(N_A, 2));
}

//==============================================================================
// Three-electron reduced matrix elements and <L^2>, <S^2>, against exact
// results:
//  1. The angular momentum operator J itself: <A||J||B> = delta_AB
//     sqrt(J(J+1)(2J+1)) for every CSF of a configuration (the sums over
//     parents between different CSFs of one configuration must cancel).
//  2. Closed 1s^2 pair plus z: <[1s^2]z||t^k||[1s^2]z'> = <z||t^k||z'>, for
//     odd (E1) and even (E2) rank.
//  3. General case (different j of the active orbital, different J of the two
//     CSFs): against the explicit expansion of each CSF over single-particle
//     m states, with the Wigner-Eckart theorem applied to each orbital.
//  4. <L^2>, <S^2> of 1s^2 nl states: l(l+1), 3/4; the S^2 and L^2 spectra
//     of 2s 2p^2 J=1/2.
//  5. Mixed states: against the sum over the full spectrum.
TEST_CASE("CI: three-electron matrix elements", "[CI][3VCI][unit]") {

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"Li", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[He]", std::nullopt, 1.0e-5, false);
  auto b_params = SplineBasis::Parameters("4spd", 20, 6, {1.0e-2}, 20.0);
  b_params.verbose = false;
  wf.formBasis(b_params);

  const auto &a = wf.core().front();
  std::vector<DiracSpinor> orbitals{a};
  for (const auto &Fn : wf.basis()) {
    if (Fn.n() > 1 && Fn.l() <= 2) {
      orbitals.push_back(Fn);
    }
  }

  const auto view = [](const std::vector<double> &c) {
    return LinAlg::View<const double>{c.data(), 0, c.size(), 1};
  };
  const auto unit = [](std::size_t n, std::size_t i) {
    std::vector<double> c(n, 0.0);
    c[i] = 1.0;
    return c;
  };

  //----------------------------------------------------------------------------
  // 1. J operator: single-particle <z||j||z> = sqrt(j(j+1)(2j+1)); rank 1, even
  Coulomb::meTable<double> j_tab;
  for (const auto &z : orbitals) {
    const auto j = 0.5 * z.twoj();
    j_tab.add(z, z, std::sqrt(j * (j + 1.0) * (2.0 * j + 1.0)));
  }

  const auto index = [&orbitals](int n, int kappa) {
    return std::find_if(orbitals.cbegin(), orbitals.cend(),
                        [n, kappa](const auto &Fn) {
                          return Fn.n() == n && Fn.kappa() == kappa;
                        })
      ->nk_index();
  };
  const std::vector<std::array<DiracSpinor::Index, 3>> configs{
    {index(2, -1), index(2, 1), index(3, 2)},
    {index(2, 1), index(3, 1), index(4, 1)},
    {index(2, -2), index(2, -2), index(3, -1)},
    {index(2, -1), index(2, -2), index(2, -2)},
    {index(2, -2), index(2, -2), index(2, -2)},
    {index(3, -3), index(3, -3), index(3, -3)},
    {index(3, -3), index(3, -3), index(4, -3)}};

  // Convention of the single-particle table: the same identity must hold for
  // two-electron CSFs (which fixes what "the" convention is)
  double worst_J2 = 0.0;
  for (const int twoJ : {0, 2, 4}) {
    const auto csfs = CI::form_CSF2s(twoJ, +1, orbitals);
    const auto J = 0.5 * twoJ;
    const auto expected = std::sqrt(J * (J + 1.0) * (2.0 * J + 1.0));
    for (std::size_t A = 0; A < csfs.size(); ++A) {
      for (std::size_t B = A; B < csfs.size(); ++B) {
        const auto cA = unit(csfs.size(), A);
        const auto cB = unit(csfs.size(), B);
        const auto rme = CI::ReducedME(view(cA), csfs, twoJ, view(cB), csfs,
                                       twoJ, j_tab, 1, +1);
        worst_J2 =
          std::max(worst_J2, std::abs(rme - (A == B ? expected : 0.0)));
      }
    }
  }

  double worst_J3 = 0.0;
  for (const auto &config : configs) {
    for (int twoJ = 1; twoJ <= 15; twoJ += 2) {
      const auto csfs_config = CI::form_CSF3s(config, twoJ);
      const std::vector<CI::CSF> csfs(csfs_config.cbegin(), csfs_config.cend());
      const auto J = 0.5 * twoJ;
      const auto expected = std::sqrt(J * (J + 1.0) * (2.0 * J + 1.0));
      for (std::size_t A = 0; A < csfs.size(); ++A) {
        for (std::size_t B = 0; B < csfs.size(); ++B) {
          const auto cA = unit(csfs.size(), A);
          const auto cB = unit(csfs.size(), B);
          const auto rme = CI::ReducedME(view(cA), csfs, twoJ, view(cB), csfs,
                                         twoJ, j_tab, 1, +1);
          worst_J3 =
            std::max(worst_J3, std::abs(rme - (A == B ? expected : 0.0)));
        }
      }
    }
  }

  //----------------------------------------------------------------------------
  // 2. Closed pair plus one electron: E1 (rank 1, odd), E2 (rank 2, even)
  const DiracOperator::E1 d{wf.grid()};
  const DiracOperator::E2 q{wf.grid()};
  const auto d_tab = Amplitudes::me_table(orbitals, &d);
  const auto q_tab = Amplitudes::me_table(orbitals, &q);

  double worst_closed = 0.0;
  const auto ia = a.nk_index();
  for (const auto &z : orbitals) {
    if (z == a)
      continue;
    const auto csf_z = CI::form_CSF3s({ia, ia, z.nk_index()}, z.twoj());
    REQUIRE(csf_z.size() == 1);
    for (const auto &zp : orbitals) {
      if (zp == a)
        continue;
      const auto csf_zp = CI::form_CSF3s({ia, ia, zp.nk_index()}, zp.twoj());
      for (const auto &[tab, op] :
           {std::pair{&d_tab,
                      &static_cast<const DiracOperator::TensorOperator &>(d)},
            std::pair{
              &q_tab,
              &static_cast<const DiracOperator::TensorOperator &>(q)}}) {
        const auto rme = CI::RME_CSF3(csf_z[0], z.twoj(), csf_zp[0], zp.twoj(),
                                      *tab, op->rank());
        worst_closed = std::max(worst_closed, std::abs(rme - tab->getv(z, zp)));
      }
    }
  }

  //----------------------------------------------------------------------------
  // 3. General case, against the m-scheme expansion of the CSFs.
  // Amplitude of |o1 m1>|o2 m2>|o3 m3> (electrons 1, 2, 3 in definite
  // orbitals) in a CSF at M = J, keyed by {o1, 2m1, o2, 2m2, o3, 2m3}
  const auto expand_m = [](const CI::CSF3 &csf, int twoJ) {
    std::map<std::array<int, 6>, double> amplitude;
    for (const auto &parent : csf.parents()) {
      const auto [x, y] = parent.pair;
      const auto z = parent.spectator;
      const auto tjx = Angular::nkindex_to_twoj(x);
      const auto tjy = Angular::nkindex_to_twoj(y);
      const auto tjz = Angular::nkindex_to_twoj(z);
      // The parent pair, antisymmetric in electrons 1 and 2 (as CSF2)
      const auto eta = x == y ? 1.0 : 1.0 / std::sqrt(2.0);
      std::vector<std::tuple<int, int, double>> pair_terms{{x, y, eta}};
      if (x != y) {
        const auto s = Angular::neg1pow_2(tjx + tjy - parent.twoJ2);
        pair_terms.push_back({y, x, -s * eta});
      }
      for (const auto &[o1, o2, c] : pair_terms) {
        const auto tj1 = Angular::nkindex_to_twoj(std::uint64_t(o1));
        const auto tj2 = Angular::nkindex_to_twoj(std::uint64_t(o2));
        for (int m1 = -tj1; m1 <= tj1; m1 += 2) {
          for (int m2 = -tj2; m2 <= tj2; m2 += 2) {
            const auto twoM2 = m1 + m2;
            const auto m3 = twoJ - twoM2;
            if (std::abs(m3) > tjz)
              continue;
            const auto cg =
              Angular::cg_2(tj1, m1, tj2, m2, parent.twoJ2, twoM2) *
              Angular::cg_2(parent.twoJ2, twoM2, tjz, m3, twoJ, twoJ);
            if (cg != 0.0) {
              amplitude[{o1, m1, o2, m2, int(z), m3}] += parent.coef * c * cg;
            }
          }
        }
      }
    }
    return amplitude;
  };

  // <X J_X J_X| sum_i t^k_q |V J_V J_V>, q = J_X - J_V, with the Wigner-Eckart
  // theorem for each orbital; then the reduced matrix element by the same
  // theorem for the whole state
  const auto brute_force_rme =
    [&expand_m](const CI::CSF3 &X, int twoJX, const CI::CSF3 &V, int twoJV,
                const Coulomb::meTable<double> &tab, int K) {
      const auto twoq = twoJX - twoJV;
      const auto amp_X = expand_m(X, twoJX);
      const auto amp_V = expand_m(V, twoJV);
      double me = 0.0;
      for (const auto &[key_V, c_V] : amp_V) {
        for (std::size_t i = 0; i < 3; ++i) {
          const auto o = std::uint64_t(key_V[2 * i]);
          const auto m = key_V[2 * i + 1];
          const auto tj = Angular::nkindex_to_twoj(o);
          // Act with t on electron i: orbital o -> op (same other electrons)
          for (const auto &[key_X, c_X] : amp_X) {
            bool others_same = true;
            for (std::size_t l = 0; l < 3; ++l) {
              if (l != i && (key_X[2 * l] != key_V[2 * l] ||
                             key_X[2 * l + 1] != key_V[2 * l + 1])) {
                others_same = false;
              }
            }
            if (!others_same)
              continue;
            const auto op = std::uint64_t(key_X[2 * i]);
            const auto mp = key_X[2 * i + 1];
            if (mp != m + twoq)
              continue;
            const auto t =
              tab.getv(DiracSpinor::Index(op), DiracSpinor::Index(o));
            if (t == 0.0)
              continue;
            const auto tjp = Angular::nkindex_to_twoj(op);
            me += c_X * c_V * Angular::neg1pow_2(tjp - mp) *
                  Angular::threej_2(tjp, 2 * K, tj, -mp, twoq, m) * t;
          }
        }
      }
      const auto wigner_eckart =
        Angular::neg1pow_2(twoJX - twoJX) *
        Angular::threej_2(twoJX, 2 * K, twoJV, -twoJX, twoq, twoJV);
      return wigner_eckart == 0.0 ? 0.0 : me / wigner_eckart;
    };

  // Pairs of configurations differing by one orbital of different j:
  // E1 (odd): 2s 2p- 3d- <-> 2s 2p- 2p+ ; E2 (even): 2s 2p- 2p+ <-> 2s 2p- 3p-
  double worst_general = 0.0;
  int num_general = 0;
  for (const auto &[config_X, config_V, tab, K] :
       std::vector<std::tuple<std::array<DiracSpinor::Index, 3>,
                              std::array<DiracSpinor::Index, 3>,
                              const Coulomb::meTable<double> *, int>>{
         {{index(2, -1), index(2, 1), index(3, 2)},
          {index(2, -1), index(2, 1), index(2, -2)},
          &d_tab,
          1},
         {{index(2, -1), index(2, 1), index(2, -2)},
          {index(2, -1), index(2, 1), index(3, 1)},
          &q_tab,
          2}}) {
    for (int twoJX = 1; twoJX <= 9; twoJX += 2) {
      for (int twoJV = 1; twoJV <= 9; twoJV += 2) {
        if (std::abs(twoJX - twoJV) > 2 * K)
          continue;
        for (const auto &X : CI::form_CSF3s(config_X, twoJX)) {
          for (const auto &V : CI::form_CSF3s(config_V, twoJV)) {
            const auto rme = CI::RME_CSF3(X, twoJX, V, twoJV, *tab, K);
            const auto bf = brute_force_rme(X, twoJX, V, twoJV, *tab, K);
            worst_general = std::max(worst_general, std::abs(rme - bf));
            ++num_general;
          }
        }
      }
    }
  }
  REQUIRE(num_general > 20);

  //----------------------------------------------------------------------------
  // 4. <L^2>, <S^2>: 1s^2 nl is a pure doublet with L = l
  double worst_LS = 0.0;
  for (const auto &z : orbitals) {
    if (z == a)
      continue;
    const auto csfs_z = CI::form_CSF3s({ia, ia, z.nk_index()}, z.twoj());
    const std::vector<CI::CSF> csfs(csfs_z.cbegin(), csfs_z.cend());
    const std::vector<double> c{1.0};
    const auto [L2, S2] = CI::expectation_L2S2(view(c), csfs, z.twoj());
    worst_LS = std::max(worst_LS, std::abs(L2 - z.l() * (z.l() + 1.0)));
    worst_LS = std::max(worst_LS, std::abs(S2 - 0.75));
  }
  // 2s 2p^2, J=1/2: the jj CSFs (over p-^2, p-p+, p+^2) mix doublets and
  // quartets, but the S^2 matrix over them must have eigenvalues 3/4 or 15/4
  // only (the states are 4P1/2, 2P1/2, 2S1/2, in the non-rel limit)
  {
    std::vector<CI::CSF> csfs;
    for (const auto &[p1, p2] :
         std::vector<std::pair<int, int>>{{1, 1}, {1, -2}, {-2, -2}}) {
      for (auto &csf :
           CI::form_CSF3s({index(2, -1), index(2, p1), index(2, p2)}, 1)) {
        csfs.push_back(std::move(csf));
      }
    }
    const auto n = csfs.size();
    // S^2 matrix in this basis, from expectation values of basis combinations
    LinAlg::Matrix<double> S2(n, n);
    for (std::size_t A = 0; A < n; ++A) {
      for (std::size_t B = 0; B < n; ++B) {
        // <A|S^2|B> = (<A+B|S^2|A+B> - <A|S^2|A> - <B|S^2|B>)/2
        auto cAB = unit(n, A);
        cAB[B] += 1.0;
        const auto SAB = CI::expectation_L2S2(view(cAB), csfs, 1).second;
        const auto SA = CI::expectation_L2S2(view(unit(n, A)), csfs, 1).second;
        const auto SB = CI::expectation_L2S2(view(unit(n, B)), csfs, 1).second;
        S2(A, B) = A == B ? SA : 0.5 * (SAB - SA - SB);
      }
    }
    // and L^2 likewise
    LinAlg::Matrix<double> L2(n, n);
    for (std::size_t A = 0; A < n; ++A) {
      for (std::size_t B = 0; B < n; ++B) {
        auto cAB = unit(n, A);
        cAB[B] += 1.0;
        const auto LAB = CI::expectation_L2S2(view(cAB), csfs, 1).first;
        const auto LA = CI::expectation_L2S2(view(unit(n, A)), csfs, 1).first;
        const auto LB = CI::expectation_L2S2(view(unit(n, B)), csfs, 1).first;
        L2(A, B) = A == B ? LA : 0.5 * (LAB - LA - LB);
      }
    }
    // Terms with J=1/2: 2S1/2, 2P1/2, 4P1/2 (ascending)
    const std::vector<double> expected_S2{0.75, 0.75, 3.75};
    const std::vector<double> expected_L2{0.0, 2.0, 2.0};
    const auto eigen_S2 = LinAlg::symmhEigensystem(S2).first;
    const auto eigen_L2 = LinAlg::symmhEigensystem(L2).first;
    REQUIRE(n == 3);
    for (std::size_t i = 0; i < n; ++i) {
      worst_LS = std::max(worst_LS, std::abs(eigen_S2[i] - expected_S2[i]));
      worst_LS = std::max(worst_LS, std::abs(eigen_L2[i] - expected_L2[i]));
    }
  }

  REQUIRE(worst_J2 < 1.0e-12);
  REQUIRE(worst_J3 < 1.0e-12);
  REQUIRE(worst_closed < 1.0e-12);
  REQUIRE(worst_general < 1.0e-12);
  REQUIRE(worst_LS < 1.0e-12);

  //----------------------------------------------------------------------------
  // 5. Mixed states vs sum over states: Li, all three electrons active
  {
    Wavefunction wf3({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                     {"Li", -1, "pointlike"}, 1.0);
    wf3.solve_core("HartreeFock", "[]", std::nullopt, 1.0e-5, false);
    auto b_params3 = SplineBasis::Parameters("3sp", 20, 6, {1.0e-2}, 20.0);
    b_params3.verbose = false;
    wf3.formBasis(b_params3);
    const auto ci_basis = CI::basis_subset(wf3.basis(), "3sp");
    Coulomb::QkTable qk;
    qk.fill(ci_basis, Coulomb::YkTable(ci_basis), 8, false);
    const auto h1 = CI::calculate_h1_table(ci_basis, {}, {}, qk, false);

    CI::PsiJPi Psi0(1, +1, ci_basis, 3);
    Psi0.solve(CI::construct_Hci(Psi0, h1, qk));
    const auto E0 = Psi0.energy(0);

    const auto test_mixed_state = [&](int twoJ, int parity,
                                      const DiracOperator::TensorOperator *h) {
      const auto h_tab = Amplitudes::me_table(ci_basis, h);
      CI::PsiJPi target(twoJ, parity, ci_basis, 3);
      const auto Hci = CI::construct_Hci(target, h1, qk);
      target.solve(Hci);
      const auto dPsi =
        CI::solve_mixed_state(Psi0, 0, target, Hci, h_tab, h->rank());
      const auto dc = dPsi.coefs(0);
      // T_I = sum_K <I||t||K> c_K: projected onto a level A, it is <A||t||0>
      const auto T =
        CI::TPsi_reduced(target.CSFs(), twoJ, Psi0, 0, h_tab, h->rank());
      double worst = 0.0;
      double max_me = 0.0;
      for (std::size_t iA = 0; iA < target.num_solutions(); ++iA) {
        const auto dE = E0 - target.energy(iA);
        const auto me =
          CI::ReducedME(target, iA, Psi0, 0, h_tab, h->rank(), h->parity());
        const auto cA = target.coefs(iA);
        double overlap = 0.0;
        double cA_T = 0.0;
        for (std::size_t i = 0; i < cA.size(); ++i) {
          overlap += cA[i] * dc[i];
          cA_T += cA[i] * T[i];
        }
        REQUIRE(cA_T == Approx(me).margin(1.0e-12));
        if (std::abs(dE) < 1.0e-10) {
          REQUIRE(std::abs(overlap) < 1.0e-10);
          continue;
        }
        max_me = std::max(max_me, std::abs(me));
        worst = std::max(worst, std::abs(dE * overlap - me));
      }
      return worst / max_me;
    };

    const DiracOperator::E1 d3{wf3.grid()};
    REQUIRE(test_mixed_state(1, -1, &d3) < 1.0e-10);
    REQUIRE(test_mixed_state(3, -1, &d3) < 1.0e-10);
    // Scalar, same J/Pi as Psi0 (singular case)
    const DiracOperator::RadialF r2{wf3.grid(), 2.0};
    REQUIRE(test_mixed_state(1, +1, &r2) < 1.0e-10);
  }
}

//==============================================================================
UnitTest::CITables UnitTest::Li_CI_tables() {
  CITables t;

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"Li", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[He]", std::nullopt, 1.0e-5, false);
  auto b_params = SplineBasis::Parameters("5spd", 20, 6, {1.0e-2}, 20.0);
  b_params.verbose = false;
  wf.formBasis(b_params);

  t.orbitals.push_back(wf.core().front());
  for (const auto &Fn : wf.basis()) {
    if (Fn.n() > 1 && Fn.l() <= 2) {
      t.orbitals.push_back(Fn);
    }
  }
  std::vector<DiracSpinor> orbitals_small; // Sk tabulated here only
  for (const auto &Fn : t.orbitals) {
    if (Fn.n() <= 3) {
      orbitals_small.push_back(Fn);
    }
  }

  t.qk.fill(t.orbitals, Coulomb::YkTable(t.orbitals), -1, false);

  for (const auto &v : t.orbitals) {
    for (const auto &w : t.orbitals) {
      if (v.kappa() == w.kappa()) {
        t.h1.add(v, w, v == w ? v.en() : 0.01 / (v.n() + w.n()));
      }
    }
  }

  const auto all_k = [](int, const DiracSpinor &, const DiracSpinor &,
                        const DiracSpinor &,
                        const DiracSpinor &) { return true; };
  const auto &qk = t.qk;
  t.Bk.fill(
    t.orbitals,
    [&qk](int k, const DiracSpinor &v, const DiracSpinor &w,
          const DiracSpinor &x,
          const DiracSpinor &y) { return 0.3 * qk.Q(k, v, w, x, y); },
    all_k, -1, false);
  t.Sk.fill(
    orbitals_small,
    [&qk](int k, const DiracSpinor &v, const DiracSpinor &w,
          const DiracSpinor &x,
          const DiracSpinor &y) { return 0.2 * qk.Q(k, v, w, x, y); },
    all_k, -1, false);
  t.dSk.fill(
    orbitals_small,
    [&qk](int k, const DiracSpinor &v, const DiracSpinor &w,
          const DiracSpinor &x,
          const DiracSpinor &y) { return 0.05 * qk.Q(k, v, w, x, y); },
    all_k, -1, false);
  t.hk = {0.1, 0.12, 0.14, 0.16};

  for (const auto &v : t.orbitals) {
    t.s1c.en[v.nk_index()] = v.en();
    t.s1c.e_sigma.insert({v.kappa(), v.en()});
    for (const auto &w : t.orbitals) {
      if (v.kappa() == w.kappa()) {
        t.s1c.S1.add(v, w, -0.02 / (v.n() * w.n()));
        t.s1c.dS1.add(v, w, -0.005 / (v.n() + w.n()));
      }
    }
  }
  t.s1c.E0 = -7.3;
  t.dE0 = 0.03;

  return t;
}

//==============================================================================
std::string UnitTest::Be_sigma2_input() {
  return "ci_basis = 6spd; J+ = 0; J- = 1; num_solutions = 2; "
         "sigma2 = true; cis2_basis = 3sp; print_details = false; "
         "qk_file = false; sk_file = false; dsk_file = false; "
         "ds1_file = false; ci_file = false;";
}

//==============================================================================
UnitTest::HoleAtom
UnitTest::hole_atom(const std::string &name, const std::string &core,
                    const std::vector<std::string> &hole_shells,
                    const std::string &basis) {
  HoleAtom atom{Wavefunction({500, 1.0e-4, 30.0, 0.33 * 30.0, "loglinear"},
                             {name, -1, "pointlike"}, 1.0)};
  atom.wf.solve_core("HartreeFock", core, std::nullopt, 1.0e-6, false);
  auto b_params = SplineBasis::Parameters(basis, 20, 6, {1.0e-2}, 20.0);
  b_params.verbose = false;
  atom.wf.formBasis(b_params);
  const auto core_shells = AtomData::core_parser(atom.wf.coreConfiguration());
  for (const auto &Fn : atom.wf.basis()) {
    const auto in_core = std::any_of(
      core_shells.cbegin(), core_shells.cend(),
      [&Fn](const auto &c) { return c.n == Fn.n() && c.l == Fn.l(); });
    (in_core ? atom.core : atom.excited).push_back(Fn);
    for (const auto &shell : hole_shells) {
      const auto s = CI::parse_configuration(shell).front();
      if (in_core && s.n == Fn.n() && s.l == Fn.l()) {
        atom.holes.push_back(Fn);
      }
    }
  }
  atom.all = qip::merge(atom.core, atom.excited);
  atom.qk.fill(atom.all, Coulomb::YkTable(atom.all), -1, false);
  HF::Breit breit{};
  breit.fill_gb(atom.all);
  const auto Bk_function = [&breit](int k, const DiracSpinor &p,
                                    const DiracSpinor &q, const DiracSpinor &r,
                                    const DiracSpinor &t) {
    return breit.Bk_abcd_2(k, p, q, r, t);
  };
  atom.Bk.fill(atom.all, Bk_function, HF::Breit::Bk_SR, 99, false);
  return atom;
}

std::vector<UnitTest::HoleAtomInput> UnitTest::hole_atoms() {
  return {{"Li", "[He]", {"1s"}, "3spd", {{1, -1}, {1, +1}, {3, +1}}},
          {"Be", "[Be]", {"1s", "2s"}, "3sp", {{1, +1}, {1, -1}, {3, -1}}},
          {"Na", "[Ne]", {"2s", "2p"}, "3spd", {{1, +1}, {3, -1}}}};
}

std::vector<UnitTest::SpinOrbital>
UnitTest::spin_orbitals(const std::vector<DiracSpinor> &orbitals) {
  std::vector<SpinOrbital> so;
  for (const auto &Fn : orbitals) {
    for (int m = -Fn.twoj(); m <= Fn.twoj(); m += 2) {
      so.push_back({Fn.nk_index(), m});
    }
  }
  return so;
}

std::size_t UnitTest::spin_orbital(const std::vector<SpinOrbital> &so,
                                   DiracSpinor::Index nk, int twom) {
  const auto it = std::find(so.cbegin(), so.cend(), SpinOrbital{nk, twom});
  assert(it != so.cend());
  return std::size_t(it - so.cbegin());
}

std::vector<double>
UnitTest::energies(const std::vector<SpinOrbital> &so,
                   const std::vector<DiracSpinor> &orbitals) {
  std::vector<double> e;
  for (const auto &s_i : so) {
    const auto nk = s_i.first;
    const auto Fn =
      std::find_if(orbitals.cbegin(), orbitals.cend(),
                   [nk](const auto &F) { return F.nk_index() == nk; });
    e.push_back(Fn->en());
  }
  return e;
}

std::pair<double, UnitTest::Determinant> UnitTest::create(std::size_t i,
                                                          Determinant D) {
  const auto pos = std::lower_bound(D.begin(), D.end(), i);
  if (pos != D.end() && *pos == i)
    return {0.0, {}};
  const auto sign = (pos - D.begin()) % 2 == 0 ? 1.0 : -1.0;
  D.insert(pos, i);
  return {sign, D};
}

std::pair<double, UnitTest::Determinant> UnitTest::annihilate(std::size_t i,
                                                              Determinant D) {
  const auto pos = std::lower_bound(D.begin(), D.end(), i);
  if (pos == D.end() || *pos != i)
    return {0.0, {}};
  const auto sign = (pos - D.begin()) % 2 == 0 ? 1.0 : -1.0;
  D.erase(pos);
  return {sign, D};
}

LinAlg::Matrix<double>
UnitTest::g_table(const std::vector<SpinOrbital> &so,
                  const std::vector<DiracSpinor> &orbitals,
                  const Coulomb::QkTable &qk, const Coulomb::WkTable *Bk) {
  const auto n = so.size();
  std::vector<const DiracSpinor *> orbital;
  for (const auto &s_i : so) {
    const auto nk = s_i.first;
    orbital.push_back(
      &*std::find_if(orbitals.cbegin(), orbitals.cend(),
                     [nk](const auto &F) { return F.nk_index() == nk; }));
  }
  LinAlg::Matrix<double> g(n * n, n * n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t l = 0; l < n; ++l) {
          // m is conserved
          if (so[i].second + so[j].second != so[k].second + so[l].second)
            continue;
          const auto g_ijkl =
            qk.g(*orbital[i], *orbital[j], *orbital[k], *orbital[l],
                 so[i].second, so[j].second, so[k].second, so[l].second);
          const auto b_ijkl =
            Bk ? Bk->g(*orbital[i], *orbital[j], *orbital[k], *orbital[l],
                       so[i].second, so[j].second, so[k].second, so[l].second) :
                 0.0;
          g(i * n + j, k * n + l) = g_ijkl + b_ijkl;
        }
      }
    }
  }
  return g;
}

double UnitTest::H_det(const Determinant &Dp, const Determinant &D,
                       const Determinant &core, const std::vector<double> &e,
                       const LinAlg::Matrix<double> &g) {
  const auto n = std::size_t(std::lround(std::sqrt(double(g.rows()))));
  // antisymmetrised <ij|g|kl> - <ij|g|lk>
  const auto gt = [&g, n](std::size_t i, std::size_t j, std::size_t k,
                          std::size_t l) {
    return g(i * n + j, k * n + l) - g(i * n + j, l * n + k);
  };
  const auto h0 = [&](std::size_t i, std::size_t j) {
    double V_HF = 0.0;
    for (const auto c : core) {
      V_HF += gt(i, c, j, c);
    }
    return (i == j ? e[i] : 0.0) - V_HF;
  };
  Determinant added, removed;
  std::set_difference(Dp.begin(), Dp.end(), D.begin(), D.end(),
                      std::back_inserter(added));
  std::set_difference(D.begin(), D.end(), Dp.begin(), Dp.end(),
                      std::back_inserter(removed));
  if (added.size() > 2)
    return 0.0;
  if (added.empty()) {
    double E = 0.0;
    for (std::size_t t = 0; t < D.size(); ++t) {
      E += h0(D[t], D[t]);
      for (std::size_t u = 0; u < t; ++u) {
        E += gt(D[t], D[u], D[t], D[u]);
      }
    }
    return E;
  }
  if (added.size() == 1) {
    const auto p = added[0];
    const auto r = removed[0];
    const auto [s1, D1] = annihilate(r, D);
    const auto [s2, D2] = create(p, D1);
    assert(D2 == Dp);
    double two_body = 0.0;
    for (const auto t : D) {
      if (t != r) {
        two_body += gt(p, t, r, t);
      }
    }
    return s1 * s2 * (h0(p, r) + two_body);
  }
  const auto [s1, D1] = annihilate(removed[0], D);
  const auto [s2, D2] = annihilate(removed[1], D1);
  const auto [s3, D3] = create(added[1], D2);
  const auto [s4, D4] = create(added[0], D3);
  assert(D4 == Dp);
  return s1 * s2 * s3 * s4 * gt(added[0], added[1], removed[0], removed[1]);
}

UnitTest::Expansion UnitTest::expand(const CI::CSF &csf, int twoJ,
                                     const Determinant &core,
                                     const std::vector<SpinOrbital> &so) {
  Expansion amplitude;
  if (csf.num_holes() == 0) {
    const auto [sign, D] = create(spin_orbital(so, csf.state(0), twoJ), core);
    amplitude[D] += sign;
    return amplitude;
  }
  const auto &c = std::get<CI::CSF2p1h>(csf.variant());
  const auto v = c.state(0);
  const auto w = c.state(1);
  const auto a = c.hole;
  const auto tjv = Angular::nkindex_to_twoj(v);
  const auto tjw = Angular::nkindex_to_twoj(w);
  const auto tja = Angular::nkindex_to_twoj(a);
  const auto eta = v == w ? 1.0 / std::sqrt(2.0) : 1.0;
  for (int mv = -tjv; mv <= tjv; mv += 2) {
    for (int mw = -tjw; mw <= tjw; mw += 2) {
      const auto M12 = mv + mw;
      const auto m = twoJ - M12;
      if (std::abs(m) > tja)
        continue;
      const auto cg = Angular::cg_2(tjv, mv, tjw, mw, c.twoJ12(), M12) *
                      Angular::cg_2(c.twoJ12(), M12, tja, m, twoJ, twoJ);
      if (cg == 0.0)
        continue;
      // h+_{a m} = (-1)^(ja - m) a_{a,-m}; then a+_w, then a+_v
      const auto [s1, D1] = annihilate(spin_orbital(so, a, -m), core);
      const auto [s2, D2] = create(spin_orbital(so, w, mw), D1);
      const auto [s3, D3] = create(spin_orbital(so, v, mv), D2);
      const auto sign = s1 * s2 * s3;
      if (sign != 0.0) {
        amplitude[D3] += eta * cg * Angular::neg1pow_2(tja - m) * sign;
      }
    }
  }
  return amplitude;
}

UnitTest::Expansion
UnitTest::particle_hole(DiracSpinor::Index w, DiracSpinor::Index a, int twoK,
                        const Determinant &core,
                        const std::vector<SpinOrbital> &so) {
  const auto tjw = Angular::nkindex_to_twoj(w);
  const auto tja = Angular::nkindex_to_twoj(a);
  Expansion amplitude;
  for (int mw = -tjw; mw <= tjw; mw += 2) {
    const auto m = twoK - mw;
    if (std::abs(m) > tja)
      continue;
    const auto cg = Angular::cg_2(tjw, mw, tja, m, twoK, twoK);
    if (cg == 0.0)
      continue;
    const auto [s1, D1] = annihilate(spin_orbital(so, a, -m), core);
    const auto [s2, D2] = create(spin_orbital(so, w, mw), D1);
    if (s1 * s2 != 0.0) {
      amplitude[D2] += cg * Angular::neg1pow_2(tja - m) * s1 * s2;
    }
  }
  return amplitude;
}

std::vector<UnitTest::Determinant>
UnitTest::Q_space(const std::vector<SpinOrbital> &so, const Determinant &core,
                  const std::vector<DiracSpinor> &holes, int twoM, int parity,
                  int num_electrons) {
  // Core spin-orbitals are the first ones; the rest are excited
  const auto num_core = core.size();
  const auto num_excited = so.size() - num_core;
  const auto in_holes = [&](std::size_t i) {
    const auto nk = Angular::index_to_nk(so[i].first);
    const auto n = nk.first;
    const auto kappa = nk.second;
    return std::any_of(holes.cbegin(), holes.cend(),
                       [n, kappa](const auto &Fa) {
                         return Fa.n() == n && Fa.kappa() == kappa;
                       });
  };
  const auto parity_of = [&](std::size_t i) {
    return Angular::parity_k(Angular::index_to_nk(so[i].first).second);
  };
  // Each r-combination of n items, in order
  const auto for_each_combination = [](std::size_t n, std::size_t r,
                                       const auto &fn) {
    if (r > n)
      return;
    std::vector<std::size_t> idx(r);
    for (std::size_t i = 0; i < r; ++i) {
      idx[i] = i;
    }
    while (true) {
      fn(idx);
      std::size_t i = r;
      while (i > 0 && idx[i - 1] == n - r + i - 1) {
        --i;
      }
      if (i == 0)
        break;
      ++idx[i - 1];
      for (std::size_t j = i; j < r; ++j) {
        idx[j] = idx[j - 1] + 1;
      }
    }
  };
  std::vector<Determinant> Q;
  for (std::size_t num_holes = 1; num_holes <= 3; ++num_holes) {
    const auto num_particles = num_holes + std::size_t(num_electrons);
    for_each_combination(
      num_core, num_holes, [&](const std::vector<std::size_t> &ih) {
        if (std::all_of(ih.cbegin(), ih.cend(), in_holes))
          return;
        for_each_combination(num_excited, num_particles,
                             [&](const std::vector<std::size_t> &ip) {
                               int twoM_D = 0;
                               int parity_D = 1;
                               Determinant D = core;
                               for (const auto i : ih) {
                                 twoM_D -= so[i].second;
                                 parity_D *= parity_of(i);
                                 D.erase(std::find(D.begin(), D.end(), i));
                               }
                               for (const auto i : ip) {
                                 twoM_D += so[num_core + i].second;
                                 parity_D *= parity_of(num_core + i);
                                 D.push_back(num_core + i);
                               }
                               if (twoM_D != twoM || parity_D != parity)
                                 return;
                               std::sort(D.begin(), D.end());
                               Q.push_back(D);
                             });
      });
  }
  return Q;
}

LinAlg::Matrix<double>
UnitTest::second_order(const std::vector<Expansion> &states,
                       const std::vector<Determinant> &Q,
                       const Determinant &core, const std::vector<double> &e,
                       const LinAlg::Matrix<double> &g) {
  const auto E0 = [&e](const Determinant &D) {
    double E = 0.0;
    for (const auto i : D) {
      E += e[i];
    }
    return E;
  };
  const auto N = states.size();
  const auto NQ = Q.size();
  // <A|H|Q>: between different determinants, H is the two-body interaction
  LinAlg::Matrix<double> V(N, NQ);
  std::vector<double> E_A(N), E_Q(NQ);
  for (std::size_t q = 0; q < NQ; ++q) {
    E_Q[q] = E0(Q[q]);
  }
  for (std::size_t i = 0; i < N; ++i) {
    E_A[i] = E0(states[i].cbegin()->first);
    for (std::size_t q = 0; q < NQ; ++q) {
      double sum = 0.0;
      for (const auto &[D, amp] : states[i]) {
        sum += amp * H_det(D, Q[q], core, e, g);
      }
      V(i, q) = sum;
    }
  }
  LinAlg::Matrix<double> H2(N, N);
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      double sum = 0.0;
      for (std::size_t q = 0; q < NQ; ++q) {
        sum += V(i, q) * V(j, q) * 0.5 *
               (1.0 / (E_A[i] - E_Q[q]) + 1.0 / (E_A[j] - E_Q[q]));
      }
      H2(i, j) = sum;
    }
  }
  return H2;
}

//==============================================================================
// Hole CSFs (two particles, one hole) with one-electron CSFs in one block,
// against an explicit expansion of the CSFs over Slater determinants: the
// determinants are formed by acting with the creation (and hole creation)
// operators on the closed-shell core, keeping the fermion signs; the
// Hamiltonian between determinants follows the Slater-Condon rules, from
// the single-particle Coulomb integrals in the m-scheme and the one-body
// h0 = e - V_HF (so that h0 + V_HF is diagonal, as the V^{N-M} CI assumes).
// Checks the CSFs (normalised, orthogonal), every matrix element of the
// block Hamiltonian (one-body, pair Coulomb, particle-hole Coulomb, and the
// one-electron / hole-CSF coupling), the reduced matrix elements of one-body
// tensor operators (E1 across blocks, M1 within), and the CI driver with the
// 'holes' option. Li (hole in 1s) and Na (holes in 2s, 2p: different hole
// orbitals in one block)
// construct_Hci with hole CSFs (two particles and a hole), with Coulomb and
// with Breit, against the Hamiltonian between Slater determinants
TEST_CASE("CI: hole CSF Hamiltonian", "[CI][unit][ciholes]") {
  for (const auto &input : UnitTest::hole_atoms()) {
    const auto atom = UnitTest::hole_atom(input.atom, input.core,
                                          input.hole_shells, input.basis);
    const auto so = UnitTest::spin_orbitals(atom.all);
    const auto e = UnitTest::energies(so, atom.all);
    const auto g = UnitTest::g_table(so, atom.all, atom.qk);
    const auto g_Breit = UnitTest::g_table(so, atom.all, atom.qk, &atom.Bk);
    // The core: its spin-orbitals come first
    UnitTest::Determinant core(UnitTest::spin_orbitals(atom.core).size());
    std::iota(core.begin(), core.end(), 0);
    const auto E_core = UnitTest::H_det(core, core, core, e, g);
    const auto E_core_Breit = UnitTest::H_det(core, core, core, e, g_Breit);
    // The CI tables: one-body = orbital energies
    Coulomb::meTable<double> h1;
    for (const auto &Fn : qip::merge(atom.excited, atom.holes)) {
      h1.add(Fn, Fn, Fn.en());
    }

    double worst = 0.0, worst_Breit = 0.0, largest_Breit = 0.0;
    for (const auto &[twoJ, parity] : input.blocks) {
      const CI::PsiJPi psi(twoJ, parity, atom.excited, 1, {}, atom.holes);
      REQUIRE(psi.has_holes());
      std::vector<UnitTest::Expansion> csfs;
      for (const auto &csf : psi.CSFs()) {
        csfs.push_back(UnitTest::expand(csf, twoJ, core, so));
      }
      const auto Hci = CI::construct_Hci(psi, h1, atom.qk);
      const auto Hci_Breit = CI::construct_Hci(psi, h1, atom.qk, &atom.Bk);
      for (std::size_t i = 0; i < psi.num_CSFs(); ++i) {
        for (std::size_t j = 0; j < psi.num_CSFs(); ++j) {
          double H = 0.0, H_Breit = 0.0;
          for (const auto &[Dp, amp_i] : csfs[i]) {
            for (const auto &[D, amp_j] : csfs[j]) {
              H += amp_i * amp_j * UnitTest::H_det(Dp, D, core, e, g);
              H_Breit +=
                amp_i * amp_j * UnitTest::H_det(Dp, D, core, e, g_Breit);
            }
          }
          if (i == j) {
            H -= E_core;
            H_Breit -= E_core_Breit;
          }
          worst = std::max(worst, std::abs(Hci(i, j) - H));
          worst_Breit =
            std::max(worst_Breit, std::abs(Hci_Breit(i, j) - H_Breit));
          largest_Breit = std::max(largest_Breit, std::abs(H_Breit - H));
        }
      }
    }
    fmt::print("{}: Hci vs determinants {:.1e}; with Breit {:.1e} (Breit "
               "terms up to {:.1e})\n",
               input.atom, worst, worst_Breit, largest_Breit);
    REQUIRE(worst < 1.0e-10);
    REQUIRE(worst_Breit < 1.0e-10);
    REQUIRE(largest_Breit > 1.0e-6);
  }
}

//==============================================================================
// L^2 and S^2 over a block with hole CSFs (L2S2_element): the eigenvalues are
// L(L+1) and S(S+1), with integer L and half-integer S (one electron), and
// expectation_L2S2 of each eigenvector (which groups the CSFs by
// configuration) returns its eigenvalue
TEST_CASE("CI: hole CSF L2 and S2", "[CI][unit][ciholes]") {
  for (const auto &input : UnitTest::hole_atoms()) {
    const auto atom = UnitTest::hole_atom(input.atom, input.core,
                                          input.hole_shells, input.basis);
    double worst_eigenvalue = 0.0, worst_expectation = 0.0;
    for (const auto &[twoJ, parity] : input.blocks) {
      CI::PsiJPi psi(twoJ, parity, atom.excited, 1, {}, atom.holes);
      const auto N = psi.num_CSFs();
      LinAlg::Matrix<double> L2(N, N), S2(N, N);
      for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
          std::tie(L2(i, j), S2(i, j)) =
            CI::L2S2_element(psi.csf(i), psi.csf(j), twoJ);
        }
      }
      for (const auto &[J2, is_S] :
           {std::pair{&L2, false}, std::pair{&S2, true}}) {
        psi.solve(*J2, 0);
        for (std::size_t n = 0; n < psi.num_solutions(); ++n) {
          // lambda = X(X+1): X integer (L), or half-integer (S)
          const auto lambda = psi.energy(n);
          const auto twoX = std::round(std::sqrt(1.0 + 4.0 * lambda) - 1.0);
          REQUIRE(std::fmod(twoX, 2.0) == (is_S ? 1.0 : 0.0));
          worst_eigenvalue = std::max(
            worst_eigenvalue, std::abs(lambda - 0.25 * twoX * (twoX + 2.0)));
          const auto [L2_c, S2_c] =
            CI::expectation_L2S2(psi.coefs(n), psi.CSFs(), twoJ);
          worst_expectation = std::max(worst_expectation,
                                       std::abs((is_S ? S2_c : L2_c) - lambda));
        }
      }
    }
    fmt::print("{}: L^2, S^2 eigenvalues vs X(X+1) {:.1e}; expectation_L2S2 "
               "{:.1e}\n",
               input.atom, worst_eigenvalue, worst_expectation);
    REQUIRE(worst_eigenvalue < 1.0e-10);
    REQUIRE(worst_expectation < 1.0e-10);
  }
}

//==============================================================================
// Reduced matrix elements between hole CSFs, and between one-electron and
// hole CSFs (RME_CSF; E1 and M1): <A||T||B> = (-1)^(J_A - J_B) <B||T||A>
TEST_CASE("CI: hole CSF RMEs", "[CI][unit][ciholes]") {
  for (const auto &input : UnitTest::hole_atoms()) {
    const auto atom = UnitTest::hole_atom(input.atom, input.core,
                                          input.hole_shells, input.basis);
    std::vector<CI::PsiJPi> blocks;
    for (const auto &[twoJ, parity] : input.blocks) {
      blocks.emplace_back(twoJ, parity, atom.excited, 1,
                          CI::ConfigurationSpace{}, atom.holes);
    }
    const auto ci_orbitals = qip::merge(atom.excited, atom.holes);
    const DiracOperator::E1 e1{atom.wf.grid()};
    const DiracOperator::M1 m1{atom.wf.grid(), PhysConst::alpha, 0.0};
    double worst = 0.0, largest = 0.0;
    std::size_t num_rme = 0;
    for (const auto *op :
         {static_cast<const DiracOperator::TensorOperator *>(&e1),
          static_cast<const DiracOperator::TensorOperator *>(&m1)}) {
      const auto tab = Amplitudes::me_table(ci_orbitals, op);
      const auto K = op->rank();
      for (const auto &A : blocks) {
        for (const auto &B : blocks) {
          if (A.parity() * B.parity() != op->parity() ||
              std::abs(A.twoJ() - B.twoJ()) > 2 * K)
            continue;
          const auto sign = Angular::neg1pow_2(A.twoJ() - B.twoJ());
          for (std::size_t i = 0; i < A.num_CSFs(); ++i) {
            for (std::size_t j = 0; j < B.num_CSFs(); ++j) {
              const auto AB =
                CI::RME_CSF(A.csf(i), A.twoJ(), B.csf(j), B.twoJ(), tab, K);
              const auto BA =
                CI::RME_CSF(B.csf(j), B.twoJ(), A.csf(i), A.twoJ(), tab, K);
              worst = std::max(worst, std::abs(AB - sign * BA));
              largest = std::max(largest, std::abs(AB));
              ++num_rme;
            }
          }
        }
      }
    }
    fmt::print("{}: {} RMEs, bra-ket symmetry {:.1e} (largest {:.1e})\n",
               input.atom, num_rme, worst, largest);
    REQUIRE(num_rme > 100);
    REQUIRE(largest > 0.1);
    REQUIRE(worst < 1.0e-12);
  }
}

//==============================================================================
// The CI driver with hole CSFs: its energies against the eigenvalues of
// construct_Hci over the same CSFs
TEST_CASE("CI: hole CSF driver", "[CI][unit][ciholes]") {
  for (const auto &input : UnitTest::hole_atoms()) {
    const auto atom = UnitTest::hole_atom(input.atom, input.core,
                                          input.hole_shells, input.basis);
    Coulomb::meTable<double> h1;
    for (const auto &Fn : qip::merge(atom.excited, atom.holes)) {
      h1.add(Fn, Fn, Fn.en());
    }
    std::string J_even, J_odd;
    for (const auto &[twoJ, parity] : input.blocks) {
      auto &J_list = parity == 1 ? J_even : J_odd;
      J_list += (J_list.empty() ? "" : ",") + CI::J_string(twoJ);
    }
    std::string holes;
    for (const auto &shell : input.hole_shells) {
      holes += (holes.empty() ? "" : ",") + shell;
    }
    const IO::InputBlock ci_input{
      "CI", "ci_basis = " + input.basis + "; num_electrons = 1; holes = " +
              holes + "; J+ = " + J_even + "; J- = " + J_odd +
              "; num_solutions = 0; qk_file = false; ci_file = false; "
              "print_details = false;"};
    const auto solution = CI::configuration_interaction(ci_input, atom.wf);
    REQUIRE(solution.levels.size() == input.blocks.size());
    REQUIRE(solution.integrals.holes.size() == atom.holes.size());
    double worst = 0.0;
    for (const auto &level : solution.levels) {
      REQUIRE(level.has_holes());
      CI::PsiJPi psi(level.twoJ(), level.parity(), atom.excited, 1, {},
                     atom.holes);
      REQUIRE(psi.num_CSFs() == level.num_CSFs());
      psi.solve(CI::construct_Hci(psi, h1, atom.qk), 0);
      for (std::size_t i = 0; i < level.num_solutions(); ++i) {
        worst = std::max(worst, std::abs(level.energy(i) - psi.energy(i)));
      }
    }
    fmt::print("{}: driver vs construct_Hci {:.1e}\n", input.atom, worst);
    REQUIRE(worst < 1.0e-9);
  }
}

//==============================================================================
// Second order with one hole leg: the one-electron / hole CSF elements,
// Sigma2_AB (the S^k with a hole leg, the terms with every internal line in
// the CI orbitals left out) plus one_body_1p_2p1h with <w|Sigma|a> between
// the other particle and the hole (Sigma_1 less the same terms, at the hole
// and the particle energies, averaged), against the second-order sum over
// the determinants outside the CI orbitals. Every excited orbital is a CI
// orbital, so those are the determinants with a hole outside the hole
// shells. RS denominators: exact
TEST_CASE("CI: Sigma_2 with one hole leg", "[CI][unit][ciholes]") {
  // The hole shells are some of the core shells: with every core shell a hole
  // shell (and every excited orbital a CI orbital), every term is left out
  const std::vector<UnitTest::HoleAtomInput> inputs{
    {"Be", "[Be]", {"2s"}, "3sp", {{1, +1}, {1, -1}, {3, -1}}},
    {"Be", "[Be]", {"1s"}, "3sp", {{1, +1}, {1, -1}, {3, -1}}},
    {"Na", "[Ne]", {"2p"}, "3sp", {{1, +1}, {3, -1}}}};
  for (const auto &input : inputs) {
    const auto atom = UnitTest::hole_atom(input.atom, input.core,
                                          input.hole_shells, input.basis);
    const auto so = UnitTest::spin_orbitals(atom.all);
    const auto e = UnitTest::energies(so, atom.all);
    const auto g = UnitTest::g_table(so, atom.all, atom.qk);
    UnitTest::Determinant core(UnitTest::spin_orbitals(atom.core).size());
    std::iota(core.begin(), core.end(), 0);

    // S^k over the CI orbitals, with hole legs, less the terms within the CI
    // orbitals (the hole shells and every excited orbital)
    const auto &holes = atom.holes;
    const MBPT::CIOrbitals ci_orbitals(holes, atom.excited);
    const Angular::SixJTable sixj(DiracSpinor::max_tj(atom.all));
    const auto Sk_function = [&](int k, const DiracSpinor &p,
                                 const DiracSpinor &q, const DiracSpinor &r,
                                 const DiracSpinor &t) {
      return MBPT::Sk_vwxy(k, p, q, r, t, atom.qk, atom.core, atom.excited,
                           sixj, MBPT::Denominators::RS, {}, 0.0, nullptr,
                           &ci_orbitals);
    };
    Coulomb::WkTable Sk;
    Sk.fill(qip::merge(atom.excited, atom.holes), Sk_function, MBPT::Sk_vwxy_SR,
            99, false);
    // <w|Sigma|a> between the excited and hole orbitals of one kappa, less
    // its terms within the CI orbitals (core lines in the hole shells)
    const auto Sigma_wa = [&](const DiracSpinor &Fw, const DiracSpinor &Fa,
                              double E) {
      return MBPT::Sigma_vw(Fw, Fa, atom.qk, atom.core, atom.excited, 99, E) -
             MBPT::Sigma_vw(Fw, Fa, atom.qk, holes, atom.excited, 99, E);
    };
    Coulomb::meTable<double> h1_Sigma;
    for (const auto &Fw : atom.excited) {
      for (const auto &Fa : holes) {
        if (Fw.kappa() != Fa.kappa())
          continue;
        h1_Sigma.add(Fw, Fa,
                     0.5 *
                       (Sigma_wa(Fw, Fa, Fa.en()) + Sigma_wa(Fw, Fa, Fw.en())));
      }
    }

    double worst = 0.0, largest = 0.0;
    for (const auto &[twoJ, parity] : input.blocks) {
      const CI::PsiJPi psi(twoJ, parity, atom.excited, 1, {}, atom.holes);
      std::vector<UnitTest::Expansion> csfs;
      for (const auto &csf : psi.CSFs()) {
        csfs.push_back(UnitTest::expand(csf, twoJ, core, so));
      }
      const auto Q = UnitTest::Q_space(so, core, holes, twoJ, parity, 1);
      REQUIRE(!Q.empty());
      const auto H2 = UnitTest::second_order(csfs, Q, core, e, g);
      for (std::size_t i = 0; i < psi.num_CSFs(); ++i) {
        for (std::size_t j = 0; j < psi.num_CSFs(); ++j) {
          if (psi.csf(i).num_holes() != 0 || psi.csf(j).num_holes() != 1)
            continue;
          const auto &x = std::get<CI::CSF1>(psi.csf(i).variant());
          const auto &vwa = std::get<CI::CSF2p1h>(psi.csf(j).variant());
          const auto Sigma2 = CI::Sigma2_AB(x, vwa, twoJ, Sk) +
                              CI::one_body_1p_2p1h(h1_Sigma, x, vwa, twoJ);
          worst = std::max(worst, std::abs(Sigma2 - H2(i, j)));
          largest = std::max(largest, std::abs(H2(i, j)));
        }
      }
    }
    fmt::print("{}, holes {}: Sigma_2 (one hole leg) vs second-order "
               "determinants {:.1e} (largest {:.1e})\n",
               input.atom, input.hole_shells.front(), worst, largest);
    REQUIRE(worst < 1.0e-10);
    REQUIRE(largest > 1.0e-7);
  }
}

//==============================================================================
// Second order with two hole legs: particle_hole_element with the Sigma_2
// pair elements, for the particle-hole pairs of the closed-shell core,
// |(w a~) K>, with different particles and different holes (then the element
// is the two-body part alone, and no spectator), against the second-order sum
// over the determinants outside the CI orbitals (with a hole in 1s or 2s).
// First order (the Coulomb element) checks the pair states. Na: the holes are
// the two 2p orbitals; 1s and 2s are not hole shells
TEST_CASE("CI: Sigma_2 with two hole legs", "[CI][unit][ciholes]") {
  const auto atom = UnitTest::hole_atom("Na", "[Ne]", {"2p"}, "3spd");
  const auto so = UnitTest::spin_orbitals(atom.all);
  const auto e = UnitTest::energies(so, atom.all);
  const auto g = UnitTest::g_table(so, atom.all, atom.qk);
  UnitTest::Determinant core(UnitTest::spin_orbitals(atom.core).size());
  std::iota(core.begin(), core.end(), 0);

  const auto &holes = atom.holes;
  const MBPT::CIOrbitals ci_orbitals(holes, atom.excited);
  const Angular::SixJTable sixj(DiracSpinor::max_tj(atom.all));
  const auto Sk_function = [&](int k, const DiracSpinor &p,
                               const DiracSpinor &q, const DiracSpinor &r,
                               const DiracSpinor &t) {
    return MBPT::Sk_vwxy(k, p, q, r, t, atom.qk, atom.core, atom.excited, sixj,
                         MBPT::Denominators::RS, {}, 0.0, nullptr,
                         &ci_orbitals);
  };
  Coulomb::WkTable Sk;
  Sk.fill(qip::merge(atom.excited, atom.holes), Sk_function, MBPT::Sk_vwxy_SR,
          99, false);
  const auto coulomb = CI::coulomb_pair_element(atom.qk);
  const auto sigma2 = CI::sigma2_pair_element(Sk);

  double worst_first = 0.0, worst_second = 0.0, largest = 0.0;
  std::size_t num_pairs = 0;
  for (const auto &Fw : atom.excited) {
    for (const auto &Fa : holes) {
      for (int twoK = std::abs(Fw.twoj() - Fa.twoj());
           twoK <= Fw.twoj() + Fa.twoj(); twoK += 2) {
        // The pairs (y b~) of the same K and parity, with y != w and b != a
        const auto parity = Fw.parity() * Fa.parity();
        std::vector<std::pair<DiracSpinor::Index, DiracSpinor::Index>> pairs;
        std::vector<UnitTest::Expansion> states;
        for (const auto &Fy : atom.excited) {
          for (const auto &Fb : holes) {
            if (Fy == Fw || Fb == Fa || Fy.parity() * Fb.parity() != parity ||
                std::abs(Fy.twoj() - Fb.twoj()) > twoK ||
                Fy.twoj() + Fb.twoj() < twoK)
              continue;
            pairs.push_back({Fy.nk_index(), Fb.nk_index()});
            states.push_back(UnitTest::particle_hole(
              Fy.nk_index(), Fb.nk_index(), twoK, core, so));
          }
        }
        if (pairs.empty())
          continue;
        const auto ket =
          UnitTest::particle_hole(Fw.nk_index(), Fa.nk_index(), twoK, core, so);
        const auto Q = UnitTest::Q_space(so, core, holes, twoK, parity, 0);
        const auto H2 =
          UnitTest::second_order(qip::merge(states, {ket}), Q, core, e, g);
        for (std::size_t i = 0; i < pairs.size(); ++i) {
          const auto [y, b] = pairs[i];
          double H1 = 0.0;
          for (const auto &[Dp, amp_i] : states[i]) {
            for (const auto &[D, amp_j] : ket) {
              H1 += amp_i * amp_j * UnitTest::H_det(Dp, D, core, e, g);
            }
          }
          const auto first = CI::particle_hole_element(
            coulomb, y, b, Fw.nk_index(), Fa.nk_index(), twoK);
          const auto second = CI::particle_hole_element(
            sigma2, y, b, Fw.nk_index(), Fa.nk_index(), twoK);
          worst_first = std::max(worst_first, std::abs(first - H1));
          worst_second =
            std::max(worst_second, std::abs(second - H2(i, pairs.size())));
          largest = std::max(largest, std::abs(H2(i, pairs.size())));
          ++num_pairs;
        }
      }
    }
  }
  fmt::print(
    "Na, core particle-hole pairs ({}): Coulomb vs determinants {:.1e}; "
    "Sigma_2 (two hole legs) vs second order {:.1e} (largest {:.1e})\n",
    num_pairs, worst_first, worst_second, largest);
  REQUIRE(num_pairs > 10);
  REQUIRE(worst_first < 1.0e-10);
  REQUIRE(worst_second < 1.0e-10);
  REQUIRE(largest > 1.0e-7);
}

//==============================================================================
// Sigma_1 with hole CSFs in the block: the driver's h1, for the excited and
// the hole orbitals alike, is the Sigma_1 of the kappa at the energy of the
// lowest CI orbital of that kappa, less its terms with every internal line in
// the CI orbitals (core lines in the hole shells, excited lines in
// ci_basis); its Q^k table holds every integral that needs (against a table
// of all integrals). The CI basis is smaller than the basis, so part of
// Sigma_1 remains
TEST_CASE("CI: Sigma_1 with hole CSFs", "[CI][unit][ciholes]") {
  Wavefunction wf({500, 1.0e-4, 30.0, 10.0, "loglinear"},
                  {"Li", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[He]", std::nullopt, 1.0e-6, false);
  auto b_params = SplineBasis::Parameters("5spd", 20, 6, {1.0e-2}, 20.0);
  b_params.verbose = false;
  wf.formBasis(b_params);

  const auto [core, excited] =
    DiracSpinor::split_by_core(wf.basis(), wf.core());
  const auto &holes = core;
  Coulomb::QkTable qk;
  qk.fill(wf.basis(), Coulomb::YkTable(wf.basis()), -1, false);

  const IO::InputBlock input{
    "CI", "ci_basis = 4spd; num_electrons = 1; holes = 1s; "
          "ci_hole_basis = 3sp; sigma1 = true; s1_basis = 5spd; J- = 1/2; "
          "num_solutions = 0; qk_file = false; ci_file = false; "
          "print_details = false;"};
  const auto solution = CI::configuration_interaction(input, wf);
  REQUIRE(solution.levels.size() == 1);
  const auto &ints = solution.integrals;

  const auto kappa_energy = [&ints](int kappa) {
    for (const auto &Fn : ints.ci_basis) {
      if (Fn.kappa() == kappa)
        return Fn.en();
    }
    return 0.0;
  };
  const auto full_Sigma = [&, &core = core, &excited = excited](
                            const DiracSpinor &v, const DiracSpinor &w) {
    return MBPT::Sigma_vw(v, w, qk, core, excited, 99, kappa_energy(v.kappa()));
  };
  const auto Sigma_in_CI = [&](const DiracSpinor &v, const DiracSpinor &w) {
    return MBPT::Sigma_vw(v, w, qk, holes, ints.ci_basis, 99,
                          kappa_energy(v.kappa()));
  };

  const auto orbitals = qip::merge(ints.ci_basis, ints.holes);
  double worst = 0.0, largest_cross = 0.0;
  for (const auto &v : orbitals) {
    for (const auto &w : orbitals) {
      if (w.kappa() != v.kappa())
        continue;
      const auto expected =
        (v == w ? v.en() : 0.0) + full_Sigma(v, w) - Sigma_in_CI(v, w);
      worst = std::max(worst, std::abs(ints.h1.getv(v, w) - expected));
      if (v.n() > 1 && w == holes.front()) {
        largest_cross = std::max(largest_cross, std::abs(expected));
      }
    }
  }

  const auto &a = holes.front();
  const auto Sigma_hole = full_Sigma(a, a) - Sigma_in_CI(a, a);
  fmt::print("h1 (excited and hole orbitals), driver vs all integrals: {:.1e}; "
             "hole Sigma_1 {:.6e} (full {:.6e}); <w|Sigma|a> up to {:.1e}\n",
             worst, Sigma_hole, full_Sigma(a, a), largest_cross);
  REQUIRE(worst < 1.0e-10);
  // For the hole: part is left out, and part remains
  REQUIRE(std::abs(Sigma_in_CI(a, a)) > 1.0e-4);
  REQUIRE(std::abs(Sigma_hole) > 1.0e-4);
  REQUIRE(largest_cross > 1.0e-5);
}

//==============================================================================
// The driver with sigma2 and hole CSFs: the S^k table holds the integrals
// with a hole leg (at most one per pair) that the block's hole CSFs need, and
// every integral (with or without a hole leg) leaves out the terms with every
// internal line in the CI orbitals (the hole shells and ci_basis); and the
// levels are those of the Hamiltonian built from the stored tables. Li, 1s
// hole, hole CSFs = single excitations from 2s only, particles in 3sp; the CI
// basis is smaller than the basis, so the terms with an internal line beyond
// it remain
TEST_CASE("CI: Sigma_2 driver with hole CSFs", "[CI][unit][ciholes]") {
  Wavefunction wf({500, 1.0e-4, 30.0, 10.0, "loglinear"},
                  {"Li", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[He]", std::nullopt, 1.0e-6, false);
  auto b_params = SplineBasis::Parameters("5spd", 20, 6, {1.0e-2}, 20.0);
  b_params.verbose = false;
  wf.formBasis(b_params);
  const auto [core, excited] =
    DiracSpinor::split_by_core(wf.basis(), wf.core());
  Coulomb::QkTable qk;
  qk.fill(wf.basis(), Coulomb::YkTable(wf.basis()), -1, false);
  const Angular::SixJTable sixj(DiracSpinor::max_tj(wf.basis()));

  const IO::InputBlock input{
    "CI", "ci_basis = 4spd; num_electrons = 1; holes = 1s; "
          "ci_hole_basis = 3sp; leading = 2s; excitations = 1; "
          "sigma1 = true; sigma2 = true; s1_basis = 5spd; cis2_basis = 4spd; "
          "denominators = RS; J+ = 1/2; J- = 1/2; num_solutions = 0; "
          "qk_file = false; sk_file = false; ci_file = false; "
          "print_details = false;"};
  const auto solution = CI::configuration_interaction(input, wf);
  const auto &ints = solution.integrals;
  REQUIRE(solution.levels.size() == 2);
  REQUIRE(ints.holes.size() == 1);
  REQUIRE(!ints.Sk.emptyQ());
  REQUIRE(!ints.config_space.unrestricted());
  REQUIRE(ints.ci_basis.size() < excited.size());
  const auto &hole = ints.holes.front();
  const auto &hole_basis = ints.holes_basis;

  // The hole CSFs of the block, directly from the configuration space (which
  // integrals with a hole leg are stored)
  const auto in = [](const std::vector<DiracSpinor> &set,
                     const DiracSpinor &n) {
    return std::find(set.cbegin(), set.cend(), n) != set.cend();
  };
  const auto in_hole_basis = [&](const DiracSpinor &n) {
    return in(hole_basis, n);
  };
  const MBPT::InCISpace hole_CSF =
    [&](const DiracSpinor &a, const DiracSpinor &n, const DiracSpinor *p) {
      if (a != hole || !in_hole_basis(n))
        return false;
      if (p != nullptr) {
        return in_hole_basis(*p) &&
               ints.config_space.included({n.nk_index(), p->nk_index()},
                                          {a.nk_index()});
      }
      return std::any_of(hole_basis.cbegin(), hole_basis.cend(),
                         [&](const DiracSpinor &q) {
                           return ints.config_space.included(
                             {n.nk_index(), q.nk_index()}, {a.nk_index()});
                         });
    };
  const MBPT::CIOrbitals ci_orbitals(ints.holes, ints.ci_basis);

  // (1) Which integrals are stored (the legs and pair basis is the whole CI
  // basis, so every integral without a hole leg is; with a hole leg, those of
  // the hole CSFs: one hole leg, the other pair its particles; one in each
  // pair, each hole with the particle of the other pair), and their values
  const auto legs = qip::merge(ints.ci_basis, ints.holes);
  std::size_t num_hole_leg = 0, num_wrong = 0, num_not_finite = 0;
  double worst = 0.0, largest = 0.0, left_out = 0.0;
  for (const auto &v : legs) {
    for (const auto &w : legs) {
      for (const auto &x : legs) {
        for (const auto &y : legs) {
          const auto holes_vw = int(v == hole) + int(w == hole);
          const auto holes_xy = int(x == hole) + int(y == hole);
          bool hole_CSFs = holes_vw <= 1 && holes_xy <= 1;
          if (holes_vw == 1 && holes_xy == 1) {
            hole_CSFs = hole_CSF(hole, v == hole ? w : v, nullptr) &&
                        hole_CSF(hole, x == hole ? y : x, nullptr);
          } else if (holes_vw == 1) {
            hole_CSFs = hole_CSF(hole, x, &y);
          } else if (holes_xy == 1) {
            hole_CSFs = hole_CSF(hole, v, &w);
          }
          const auto [k0, k1] = MBPT::k_minmax_S(v, w, x, y);
          for (int k = k0; k <= k1; ++k) {
            const auto selected = hole_CSFs && MBPT::Sk_vwxy_SR(k, v, w, x, y);
            if (ints.Sk.contains(k, v, w, x, y) != selected) {
              ++num_wrong;
            }
            if (!selected)
              continue;
            const auto stored = ints.Sk.Q(k, v, w, x, y);
            if (!std::isfinite(stored)) {
              ++num_not_finite;
            }
            const auto expected = MBPT::Sk_vwxy(
              k, v, w, x, y, qk, core, excited, sixj, MBPT::Denominators::RS,
              {}, 0.0, nullptr, &ci_orbitals);
            worst = std::max(worst, std::abs(stored - expected));
            if (holes_vw + holes_xy > 0) {
              ++num_hole_leg;
              largest = std::max(largest, std::abs(expected));
            } else {
              // Without a hole leg the terms within the CI orbitals are left
              // out too
              const auto all_terms = MBPT::Sk_vwxy(
                k, v, w, x, y, qk, core, excited, sixj, MBPT::Denominators::RS,
                {}, 0.0, nullptr, nullptr);
              left_out = std::max(left_out, std::abs(all_terms - expected));
            }
          }
        }
      }
    }
  }
  fmt::print("S^k: {} with hole legs stored ({} not finite), {} wrongly "
             "(un)stored; driver vs direct {:.1e} (hole legs up to {:.1e}; "
             "terms left out without a hole leg up to {:.1e})\n",
             num_hole_leg, num_not_finite, num_wrong, worst, largest, left_out);
  REQUIRE(num_hole_leg > 0);
  REQUIRE(num_not_finite == 0);
  REQUIRE(num_wrong == 0);
  REQUIRE(worst < 1.0e-14);
  REQUIRE(largest > 1.0e-6);
  REQUIRE(left_out > 1.0e-8);

  // (2) The levels are those of the Hamiltonian from the stored tables
  double worst_level = 0.0;
  for (const auto &level : solution.levels) {
    REQUIRE(level.has_holes());
    CI::PsiJPi psi(level.twoJ(), level.parity(), ints.ci_basis, 1,
                   ints.config_space, ints.holes, ints.holes_basis);
    REQUIRE(psi.num_CSFs() == level.num_CSFs());
    psi.solve(CI::construct_Hci(psi, ints), 0);
    for (std::size_t i = 0; i < level.num_solutions(); ++i) {
      worst_level =
        std::max(worst_level, std::abs(level.energy(i) - psi.energy(i)));
    }
  }
  fmt::print("Levels: driver vs construct_Hci {:.1e}\n", worst_level);
  REQUIRE(worst_level < 1.0e-9);
}
