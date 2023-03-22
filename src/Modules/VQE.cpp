#include "Modules/VQE.hpp"
#include "Coulomb/Coulomb.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include "qip/Vector.hpp"
#include <random>
#include <vector>

namespace Module {

void VQE(const IO::InputBlock &input, const Wavefunction &wf) {

  // Check input options:
  input.check({{"basis", "Basis used for CI expansion; must be a sub-set of "
                         "MBPT basis [default: 5spd]"},
               {"twoJ", "2*J angular symmetry for CSFs [0]"},
               {"parity", "parity (+/-1) symmetry for CSFs [+1]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  using namespace std::string_literals;

  // Determine the sub-set of basis to use in CI:
  const auto basis_string = input.get("basis", "5spd"s);
  // Select from wf.basis() [MBPT basis], those which match input 'basis_string'
  std::vector<DiracSpinor> ci_sp_basis;
  const auto nk_list = AtomData::n_kappa_list(basis_string);
  for (const auto &nk : nk_list) {
    const auto select_orbital = [&wf, nk](const auto &Fn) {
      const auto [max_n, kappa] = nk;
      return Fn.n() <= max_n && Fn.kappa() == kappa &&
             !wf.isInCore(Fn.n(), Fn.kappa());
    };
    qip::insert_into_if(wf.basis(), &ci_sp_basis, select_orbital);
  }

  // Print info re: basis to screen:
  std::cout << "\nUsing " << DiracSpinor::state_config(ci_sp_basis) << " = "
            << ci_sp_basis.size()
            << " relativistic single-particle orbitals:\n";
  std::cout << "n  kappa (j   l pi) index  string  E(nk)\n";
  for (const auto &Fn : ci_sp_basis) {
    fmt::print("{:<2}   {:+1}  ({}/2 {} {:+1d})   {:3}  {:4s}    {:<+10.4}\n",
               Fn.n(), Fn.kappa(), Fn.twoj(), Fn.l(), Fn.parity(),
               Fn.nk_index(), Fn.shortSymbol(), Fn.en());
  }
  std::cout << "\n'index' is unique index that maps {n,kappa} to integer\n"
            << "index(n,k) := n^2 - 2n + 1 + 2*|k| - 2   ; for k<0\n"
            << "index(n,k) := n^2 - 2n + 1 + 2*|k| - 1   ; for k>0\n";
  std::cout << "(It is *not* an efficient packing..)\n\n";

  //----------------------------------------------------------------------------

  std::cout << "\nSingle-particle matrix elements:\n";
  const auto sigma = wf.Sigma();

  // We could calculate Sigma^1 matrix elements directly
  // But it's much faster to first calculate the Sigma matrix
  // For Sigma^2 (not currently implemented), we have to do direct calculation
  // MBPT::CorrelationPotential Sigma(wf.vHF(), wf.basis(), {}, {});

  std::cout << "Diagonal:\n";
  std::cout << "v      |h^HF|_vv  |𝛴^1|_vv\n";
  for (const auto &v : ci_sp_basis) {
    const auto h0vv = v.en();
    const auto Svv = sigma ? v * sigma->SigmaFv(v) : 0.0;
    fmt::print("{:4s}  {:+9.6f}  {:+11.8f}\n", v.shortSymbol(), h0vv, Svv);
  }
  std::cout << "\nnb: since our basis is eigenstates of single-particle \n"
               "Hartree-Fock Hamiltonian, |h^HF|_vv = E_v and |h^HF|_vw = 0\n";

  const int max_count = 20;

  // if we don't have Sigma, these are all zero
  if (sigma) {
    fmt::print("\nOff-diagonal (just first ~{}):\n", max_count);
    std::cout << "v     w      |𝛴^1|_vw\n";
    int count = 0;
    for (const auto &v : ci_sp_basis) {
      for (const auto &w : ci_sp_basis) {
        // Symmetric, so only do half:
        if (w <= v)
          continue;
        // Angular selection rule:
        if (w.kappa() != v.kappa())
          continue;
        const auto Svw = sigma ? v * sigma->SigmaFv(w) : 0.0;
        fmt::print("{:4s}  {:4s}   {:+11.8f}\n", v.shortSymbol(),
                   w.shortSymbol(), Svw);
        ++count;
      }
      if (count > max_count)
        break;
    }
  }

  //----------------------------------------------------------------------------

  std::cout << "\nCalculate two-body Coulomb integrals: Q^k_abcd\n";
  Coulomb::QkTable qk;
  const auto qfname =
      wf.identity() + DiracSpinor::state_config(ci_sp_basis) + ".qk";
  const auto read_from_file_ok = qk.read(qfname);
  if (!read_from_file_ok) {
    // if didn't find Qk file to read in, calculate from scratch:
    const Coulomb::YkTable yk(ci_sp_basis);
    qk.fill(ci_sp_basis, yk);
    qk.write(qfname);
  }

  fmt::print("\nPrint a few randomly-chosen Coulomb integrals:\n");
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<std::size_t> rand_index(0,
                                                        ci_sp_basis.size() - 1);
  fmt::print("k  a    b    c    d      R^k        Q^k        W^k\n");
  for (int count = 0; count < max_count;) {
    const auto &a = ci_sp_basis.at(rand_index(gen));
    const auto &b = ci_sp_basis.at(rand_index(gen));
    const auto &c = ci_sp_basis.at(rand_index(gen));
    const auto &d = ci_sp_basis.at(rand_index(gen));
    const auto [kmin, kmax] = Coulomb::k_minmax_W(a, b, c, d);
    // nb: there can be non-zero Ws if Q is zero!
    // can k+=2 for Q, due to parity, but not for W
    for (int k = kmin; k <= kmax; ++k) {
      const auto Rk = qk.R(k, a, b, c, d);
      const auto Qk = qk.Q(k, a, b, c, d);
      const auto Wk = qk.W(k, a, b, c, d);
      fmt::print("{}  {:4s},{:4s},{:4s},{:4s}  {:+10.3e} {:+10.3e} {:+10.3e}\n",
                 k, a.shortSymbol(), b.shortSymbol(), c.shortSymbol(),
                 d.shortSymbol(), Rk, Qk, Wk);
      ++count;
    }
  }
  std::cout << "...\n";

  //----------------------------------------------------------------------------
  // List CSFs.. for now, two-particle only
  const auto twoJ = input.get("twoJ", 0);
  const auto parity = input.get("parity", 1);

  if (std::abs(parity) != 1) {
    std::cout << "parity must be +/-1\n";
    return;
  }
  if (twoJ < 0) {
    std::cout << "twoJ must >=0\n";
    return;
  }
  if (twoJ % 2 != 0) {
    std::cout << "twoJ must be even for two-electron CSF\n";
  }

  std::cout << "\nListing all 2-particle CSFs (only print first few):\n";
  fmt::print("For J={}/2, pi={:+}\n", twoJ, parity);
  int count_CSFs = 0;
  int count_projs = 0;
  for (const auto &v : ci_sp_basis) {
    for (const auto &w : ci_sp_basis) {
      // Symmetric; only show half
      if (w < v)
        continue;
      // Parity symmetry:
      if (v.parity() * w.parity() != parity)
        continue;
      // J triangle:
      if (v.twoj() + w.twoj() < twoJ || std::abs(v.twoj() - w.twoj()) > twoJ)
        continue;

      count_CSFs++;
      if (count_CSFs < max_count)
        std::cout << v << " " << w << "\n";
      double sum = 0.0;
      for (int two_m_v = -v.twoj(); two_m_v <= v.twoj(); two_m_v += 2) {
        const auto two_m_w = twoJ - two_m_v;
        if (std::abs(two_m_w) > w.twoj())
          continue;
        const auto cgc =
            Angular::cg_2(v.twoj(), two_m_v, w.twoj(), two_m_w, twoJ, twoJ);
        if (count_CSFs < max_count)
          fmt::print("  m = ({:+}/2, {:+}/2) , c_d = {:+.5f}\n", two_m_v,
                     two_m_w, cgc);
        sum += cgc * cgc;
        count_projs++;
      }
      if (count_CSFs < max_count)
        std::cout << "Sum c_d^2 = " << sum << "\n";
    }
  }
  std::cout << "...\n";

  std::cout << "\nTotal CSFs: " << count_CSFs << "\n";
  std::cout << "Total Projections: " << count_projs << "\n";
  std::cout << "*We don't usually need *all* CSFs in expansion;\n"
               "often reasonable to include all single-, but fewer double-\n"
               "excitations from a small set of leading configurations...\n";
}

} // namespace Module
