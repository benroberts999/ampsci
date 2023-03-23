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
               {"twoJ", "2*J angular symmetry for CSFs [default: 0]"},
               {"parity", "parity (+/-1) symmetry for CSFs [default: +1]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  // Format for parsing input:
  // input.get("input_option", default_value);
  // If "input_option" is present in the input file, this will return that value
  // If it is not present, will instead return 'default_value'.

  // Determine the sub-set of basis to use in CI:
  const auto basis_string = input.get("basis", std::string{"5spd"});
  // Select from wf.basis() [MBPT basis], those which match input 'basis_string'
  const std::vector<DiracSpinor> ci_sp_basis =
      basis_subset(wf.basis(), basis_string, wf.FermiLevel());

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

  // Correlation potential:
  const auto Sigma = wf.Sigma();
  // We could calculate Sigma^1 matrix elements directly
  // But it's much faster to first calculate the Sigma matrix
  // For Sigma^2 (not currently implemented), we have to do direct calculation
  // MBPT::CorrelationPotential Sigma(wf.vHF(), wf.basis(), {}, {});

  std::cout << "Diagonal:\n";
  std::cout << "v      |h^HF|_vv  |𝛴^1|_vv\n";
  for (const auto &v : ci_sp_basis) {
    const auto h0_vv = v.en();
    const auto Sigma_vv = Sigma ? v * Sigma->SigmaFv(v) : 0.0;
    fmt::print("{:4s}  {:+9.6f}  {:+11.8f}\n", v.shortSymbol(), h0_vv,
               Sigma_vv);
  }
  std::cout << "\nnb: since our basis is eigenstates of single-particle \n"
               "Hartree-Fock Hamiltonian, |h^HF|_vv = E_v and |h^HF|_vw = 0\n";

  const int max_count = 20;

  // if we don't have Sigma, these are all zero
  if (Sigma) {
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
        const auto Sigma_vw = Sigma ? v * Sigma->SigmaFv(w) : 0.0;
        fmt::print("{:4s}  {:4s}   {:+11.8f}\n", v.shortSymbol(),
                   w.shortSymbol(), Sigma_vw);
        ++count;
      }
      if (count > max_count)
        break;
    }
  }

  //----------------------------------------------------------------------------

  std::cout << "\nCalculate two-body Coulomb integrals: Q^k_abcd\n";
  // Lookup table; stores all qk's
  Coulomb::QkTable qk;
  const auto qk_filename =
      wf.identity() + DiracSpinor::state_config(ci_sp_basis) + ".qk";
  // Try to read from disk (may already have calculated Qk)
  const auto read_from_file_ok = qk.read(qk_filename);
  if (!read_from_file_ok) {
    // if didn't find Qk file to read in, calculate from scratch:
    const Coulomb::YkTable yk(ci_sp_basis);
    qk.fill(ci_sp_basis, yk);
    qk.write(qk_filename);
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
  // Count all projections; only print a few o them:
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

      ++count_CSFs;
      if (count_CSFs < max_count)
        std::cout << v << " " << w << "\n";
      // Each indevidual m projection:
      double sum_c2 = 0.0; // just to check: should \sum(c^2)=1
      for (int two_m_v = -v.twoj(); two_m_v <= v.twoj(); two_m_v += 2) {
        const auto two_m_w = twoJ - two_m_v;
        if (std::abs(two_m_w) > w.twoj())
          continue;
        const auto cgc =
            Angular::cg_2(v.twoj(), two_m_v, w.twoj(), two_m_w, twoJ, twoJ);
        if (count_CSFs < max_count) {
          fmt::print("  m = ({:+}/2, {:+}/2) , c_d = {:+.5f}\n", two_m_v,
                     two_m_w, cgc);
        }
        sum_c2 += cgc * cgc;
        ++count_projs;
      }
      if (count_CSFs < max_count) {
        std::cout << "Sum c_d^2 = " << sum_c2 << " (should=1)\n";
      }
    }
  }
  std::cout << "...\n";

  std::cout << "\nTotal CSFs: " << count_CSFs << "\n";
  std::cout << "Total Projections: " << count_projs << "\n";
  std::cout << "*We don't usually need *all* CSFs in expansion;\n"
               "often reasonable to include all single-, but fewer double-\n"
               "excitations from a small set of leading configurations...\n";
}

//==============================================================================
std::vector<DiracSpinor> basis_subset(const std::vector<DiracSpinor> &basis,
                                      const std::string &subset_string,
                                      double eFermi) {
  std::vector<DiracSpinor> subset;
  const auto nk_list = AtomData::n_kappa_list(subset_string);
  for (const auto &nk : nk_list) {
    const auto select_orbital = [eFermi, nk](const auto &Fn) {
      const auto [max_n, kappa] = nk;
      return Fn.n() <= max_n && Fn.kappa() == kappa && Fn.en() > eFermi;
    };
    qip::insert_into_if(basis, &subset, select_orbital);
  }
  return subset;
}

} // namespace Module
