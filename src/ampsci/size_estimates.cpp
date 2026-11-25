#include "size_estimates.hpp"
#include "CI/CI_Integrals.hpp"
#include "CI/CSF.hpp"
#include "CI/ConfigurationInteraction.hpp"
#include "CI/ConfigurationSpace.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/QkTable.hpp"
#include "HF/Breit.hpp"
#include "IO/ChronoTimer.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/format.hpp"
#include "fmt/ostream.hpp"
#include "qip/String.hpp"
#include "qip/omp.hpp"
#include "version/version.hpp"
#include <algorithm>
#include <memory>

//==============================================================================
SizeEstimates estimate_sizes(const IO::InputBlock &input) {
  using namespace std::string_literals;

  SizeEstimates out;
  out.core_string = input.get({"HartreeFock"}, "core", "[]"s);
  out.basis_string = input.get({"Basis"}, "states", ""s);

  // Orbitals as labels only: each n, kappa of the basis, with a placeholder
  // energy that puts the core below the Fermi level
  const auto grid =
    std::make_shared<const Grid>(GridParameters{10, 1.0e-4, 250.0});
  const auto core_configs = AtomData::core_parser(out.core_string);
  const auto is_core = [&core_configs](int n, int l) {
    return std::any_of(core_configs.cbegin(), core_configs.cend(),
                       [n, l](const auto &c) { return c.n == n && c.l == l; });
  };
  constexpr double e_core = -2.0;
  constexpr double e_excited = -1.0;
  constexpr double e_Fermi = -1.5;
  int N_max_core = 0;
  std::vector<DiracSpinor> basis;
  const auto nlmax_in = input.get({"Basis"}, "nlmax", std::vector<int>{});
  const auto nlmax = nlmax_in.size() == 2 ?
                       std::pair{nlmax_in[0], nlmax_in[1]} :
                       std::pair{0, -1};
  for (const auto &nk : SplineBasis::states_list(out.basis_string, nlmax)) {
    const auto n_max = nk.n;
    const auto kappa = nk.k;
    const auto l = Angular::l_k(kappa);
    for (int n = l + 1; n <= n_max; ++n) {
      auto &Fn = basis.emplace_back(n, kappa, grid);
      const auto core = is_core(n, l);
      Fn.en() = core ? e_core : e_excited;
      if (core) {
        ++out.num_core;
        N_max_core = std::max(N_max_core, n);
      }
    }
  }
  out.num_orbitals = basis.size();

  out.Feynman_GB = Feynman_memory_GB(input);

  // The integral counts loop over all sets of four orbitals, and can take a
  // while: each reports as it starts, and the time it took (timer, at the end
  // of its scope)

  // Structure radiation: its own module, or a StructureRadiation{} block
  // within another module. It forms the full Q^k table
  for (const auto &block : input.blocks()) {
    if (qip::ci_wc_compare(block.name(), "Module::StrucRad*") ||
        qip::ci_wc_compare(block.name(), "Module::StructureRad*") ||
        block.getBlock("StructureRadiation")) {
      out.qk_all_required = true;
    }
  }

  // All Q^k over the full basis, all k: counted only if required, since it
  // is the slowest count
  Coulomb::QkTable qk;
  if (out.qk_all_required) {
    std::cout << "Counting all Q^k (" << basis.size() << " orbitals).. "
              << std::flush;
    const IO::ChronoTimer timer;
    out.qk_all = qk.count_non_zero_integrals(basis, 999, e_Fermi)[0];
  }

  // CI options, with the same defaults as configuration_interaction. Without
  // a CI{} block, the CI and MBPT counts are for those defaults, but are not
  // required
  const auto have_CI = input.getBlock("CI") != std::nullopt;
  const auto ci_input = input.getBlock("CI").value_or(IO::InputBlock{"CI"});
  out.ci_basis_string = ci_input.get("ci_basis", "20spdf"s);
  out.num_electrons = ci_input.get("num_electrons", 2);
  const auto include_Sigma1 = ci_input.get("sigma1", false);
  // Sigma_2 and the two-body Breit are two-electron terms
  const auto include_Sigma2 =
    ci_input.get("sigma2", false) && out.num_electrons > 1;
  const auto n_min_core = ci_input.get("n_min_core", 1);
  // As the table fills: a cut of zero or less means all k
  const auto max_k = ci_input.get("max_k", 10);
  const auto k_cut = max_k <= 0 ? 999ul : std::size_t(max_k);
  const auto sigma1_dSdE = ci_input.get("sigma1_dSdE", false);

  const auto ci_orbitals =
    CI::basis_subset(basis, out.ci_basis_string, out.core_string);
  out.num_ci_orbitals = ci_orbitals.size();

  // Q^k among the CI orbitals
  {
    std::cout << "Counting CI Q^k (" << ci_orbitals.size() << " orbitals).. "
              << std::flush;
    const IO::ChronoTimer timer;
    out.qk_CI = qk.count_non_zero_integrals(ci_orbitals, k_cut, e_Fermi)[0];
    out.qk_CI_required = have_CI;
  }

  // Q^k with one or two core orbitals, for Sigma_1 and/or Sigma_2. Sigma_1
  // needs them unless the correlation potential supplies it (and, with the
  // dSigma/dE correction, its derivative)
  const auto have_Sigma = input.getBlock("Correlations") != std::nullopt;
  const auto Sigma_derivative =
    input.get({"Correlations"}, "derivative", false);
  const auto need_s1_qk =
    include_Sigma1 && (!have_Sigma || (sigma1_dSdE && !Sigma_derivative));
  {
    const auto mbpt_string =
      need_s1_qk ? ci_input.get("s1_basis", ""s) :
                   ci_input.get("s2_basis", ci_input.get("s1_basis", ""s));
    std::vector<DiracSpinor> mbpt_orbitals;
    for (const auto &Fn : CI::basis_subset(basis, mbpt_string)) {
      // core below n_min_core is not polarised
      if (Fn.en() == e_core && Fn.n() < n_min_core)
        continue;
      mbpt_orbitals.push_back(Fn);
    }
    std::cout << "Counting MBPT Q^k (" << mbpt_orbitals.size()
              << " orbitals).. " << std::flush;
    const IO::ChronoTimer timer;
    out.qk_MBPT = qk.count_non_zero_integrals(mbpt_orbitals, k_cut, e_Fermi)[2];
    out.qk_MBPT_required = have_CI && (need_s1_qk || include_Sigma2);
  }

  // Count of Coulomb-like integrals, for a table with symmetry of Table, with
  // the given selection rule; k up to the cut (as CoulombTable::fill)
  const auto count_integrals = [k_cut](const auto &table,
                                       const std::vector<DiracSpinor> &orbitals,
                                       const auto &selection_rule) {
    const auto k_max =
      std::min(std::size_t(DiracSpinor::max_tj(orbitals)), k_cut);
    std::size_t count = 0;
    for (std::size_t k = 0; k <= k_max; ++k) {
      for (const auto &a : orbitals) {
        for (const auto &b : orbitals) {
          for (const auto &c : orbitals) {
            for (const auto &d : orbitals) {
              if (table.is_NormalOrdered(a, b, c, d) &&
                  selection_rule(int(k), a, b, c, d)) {
                ++count;
              }
            }
          }
        }
      }
    }
    return count;
  };

  // B^k over Breit_basis, if Breit is included in HF
  const auto Breit_HF = input.get({"HartreeFock"}, "Breit", false) ||
                        input.get({"HartreeFock"}, "Breit", 0.0) != 0.0;
  if (have_CI && Breit_HF && ci_input.get("Breit2", true) &&
      out.num_electrons > 1) {
    const auto Breit_orbitals = CI::basis_subset(
      ci_orbitals,
      ci_input.get("Breit_basis", std::to_string(N_max_core + 6) + "spdf"));
    std::cout << "Counting B^k (" << Breit_orbitals.size() << " orbitals).. "
              << std::flush;
    const IO::ChronoTimer timer;
    out.bk =
      count_integrals(Coulomb::WkTable{}, Breit_orbitals, &HF::Breit::Bk_SR);
  }

  // CI blocks, dense: the matrix, plus a copy while diagonalising. Sparse:
  // the stored elements (lower triangle), at most the pairs of CSFs that can
  // interact; the Davidson vectors are small next to it
  const auto davidson = ci_input.get("davidson", out.num_electrons >= 3);
  const auto J_pi_list = CI::requested_J_pi(ci_input, out.num_electrons);
  const CI::ConfigurationSpace config_space{
    ci_orbitals, out.num_electrons,
    ci_input.get("leading", std::vector<std::string>{}),
    ci_input.get("excitations", std::vector<std::string>{"1", "2"})};
  std::cout << "Forming CSFs (" << J_pi_list.size() << " J/pi).. "
            << std::flush;
  const IO::ChronoTimer timer;
  for (const auto &[twoJ, parity] : J_pi_list) {
    // Number of CSFs, and of pairs that can interact (sparse elements)
    const auto CSFs =
      CI::form_CSFs(twoJ, parity, ci_orbitals, out.num_electrons, config_space);
    const auto num_CSFs = CSFs.size();
    const auto num_pairs = CI::num_interacting_pairs(CSFs);
    const auto N = double(num_CSFs);
    // Sparse: each stored element is a value (double) and a column (32 bit)
    const auto matrix_GB =
      davidson ? 12.0 * double(num_pairs) * 1.0e-9 : 2.0 * 8.0 * N * N * 1.0e-9;
    out.ci_blocks.push_back({twoJ, parity, num_CSFs, davidson, matrix_GB});
  }

  return out;
}

//==============================================================================
void print_size_estimates(const SizeEstimates &estimates, std::ostream &os) {

  const auto &e = estimates;
  fmt::print(os, "Basis: {}, {} orbitals ({} core: {})\n", e.basis_string,
             e.num_orbitals, e.num_core, e.core_string);

  // Tables the run requires are marked, and counted in the total
  double total_GB = 0.0;
  fmt::print(os, "\n{:<22} {:>14} {:>10}\n", "Integrals", "count", "memory");
  const auto print_table = [&os, &total_GB](const std::string &name,
                                            std::size_t count, bool required) {
    const auto memory_GB = Coulomb::estimate_memory_GB(count);
    fmt::print(os, "{:<22} {:>14} {:>7.1f} GB\n", name, count, memory_GB);
    if (required) {
      total_GB += memory_GB;
    }
  };
  if (e.qk_all_required) {
    print_table("Q^k (all; SR)", e.qk_all, true);
  }
  print_table("Q^k (CI)", e.qk_CI, e.qk_CI_required);
  print_table("Q^k (MBPT, 1-2 core)", e.qk_MBPT, e.qk_MBPT_required);
  if (e.bk > 0) {
    print_table("B^k (Breit)", e.bk, true);
  }

  if (e.Feynman_GB > 0.0) {
    fmt::print(os, "{:<22} {:>14} {:>7.1f} GB\n", "Sigma (Feynman)", "",
               e.Feynman_GB);
    total_GB += e.Feynman_GB;
  }

  // CI blocks are solved one at a time: the largest counts in the total
  double largest_block_GB = 0.0;
  if (!e.ci_blocks.empty()) {
    fmt::print(os, "\nCI: {}, {} orbitals, {} valence electrons\n",
               e.ci_basis_string, e.num_ci_orbitals, e.num_electrons);
    fmt::print(os, "{:<5} {:>2} {:>10} {:>14}\n", "J", "pi", "N_CSF",
               "matrix (peak)");
    for (const auto &block : e.ci_blocks) {
      fmt::print(os, "{:<5} {:>2} {:>10} {:>11.1f} GB{}\n",
                 CI::J_string(block.twoJ), block.parity == 1 ? "+" : "-",
                 block.num_CSFs, block.matrix_GB,
                 block.davidson ? "  (sparse)" : "  (dense)");
      largest_block_GB = std::max(largest_block_GB, block.matrix_GB);
    }
  }
  fmt::print(os, "\nTotal: ~{:.1f} GB in memory\n",
             total_GB + largest_block_GB);
}

//==============================================================================

double Feynman_memory_GB(const IO::InputBlock &input) {
  using namespace std::string_literals;

  // Same options and defaults as the run
  const auto all_order = input.get({"Correlations"}, "all_order", false);
  if (!input.getBlock({"Correlations"}) ||
      !input.get({"Correlations"}, "Feynman", all_order)) {
    return 0.0;
  }

  // Sigma sub-grid: every stride-th point of the radial grid
  const auto r0 = input.get({"Grid"}, "r0", 1.0e-6);
  const auto rmax = input.get({"Grid"}, "rmax", 120.0);
  const auto du = input.get<double>({"Grid"}, "du");
  const auto num_points = du ? 0ul : input.get({"Grid"}, "num_points", 2000ul);
  const auto b = input.get({"Grid"}, "b", rmax / 3.0);
  const auto grid_type = (b <= r0 || b >= rmax) ?
                           "logarithmic" :
                           input.get({"Grid"}, "type", "loglinear"s);
  const Grid grid{
    GridParameters{num_points, r0, rmax, b, grid_type, du ? *du : 0}};
  const auto sigma_rmin = input.get({"Correlations"}, "rmin", 1.0e-4);
  const auto sigma_rmax =
    std::min(input.get({"Correlations"}, "rmax", 30.0), grid.rmax());
  const auto default_stride =
    std::max(2, int(grid.getIndex(30.0) - grid.getIndex(1.0e-4)) / 200);
  const auto stride = std::size_t(
    std::max(1, input.get({"Correlations"}, "stride", default_stride)));
  const auto subgrid_points =
    (grid.getIndex(sigma_rmax) - grid.getIndex(sigma_rmin)) / stride + 1;

  // Core orbitals in the polarisation loops; partial waves; multipoles
  const auto n_min_core = input.get({"Correlations"}, "n_min_core", 1);
  std::size_t num_core = 0;
  int l_max_core = 0;
  for (const auto &shell :
       AtomData::core_parser(input.get({"HartreeFock"}, "core", "[]"s))) {
    l_max_core = std::max(l_max_core, shell.l);
    if (shell.n >= n_min_core) {
      num_core += shell.l == 0 ? 1 : 2;
    }
  }
  const auto l_max = input.get({"Correlations"}, "lmax", 6);
  const auto num_kappas = std::size_t(2 * std::max(l_max, l_max_core)) + 1;
  const auto num_k = num_kappas + 1;

  // Frequency grid: w0 * ratio^i, up to twice the deepest core energy
  // included (taken as hydrogen-like), and at least 5000 au
  const auto imag_omega =
    input.get({"Correlations"}, "imag_omega", std::vector{0.01, 1.5});
  const auto w0 = imag_omega.size() == 2 ? imag_omega[0] : 0.01;
  const auto w_ratio = imag_omega.size() == 2 ? imag_omega[1] : 1.5;
  const auto Z = AtomData::atomic_Z(input.get({"Atom"}, "Z", "H"s));
  const auto e_core_deepest =
    0.5 * double(Z * Z) / double(n_min_core * n_min_core);
  const auto w_max = std::max(2.0 * w_ratio * e_core_deepest, 5000.0);
  const auto num_w =
    std::size_t(std::log(w_ratio * w_max / w0) / std::log(w_ratio)) + 2;

  // With the small components (Breit, or include_G), a Green's function has
  // four blocks, and a spinor two components
  const auto include_G = input.get({"Correlations"}, "include_G", false) ||
                         input.get({"HartreeFock"}, "Breit", false) ||
                         input.get({"HartreeFock"}, "Breit", 0.0) != 0.0;
  const auto num_blocks = include_G ? 4ul : 1ul;
  const auto num_components = include_G ? 2ul : 1ul;

  // Number of complex sub-grid matrices held:
  // Q*Pi*Q, and a second table without screening (for the fk)
  auto num_matrices = 2 * num_w * num_k;
  const auto exchange = qip::ci_compare(
    input.get({"Correlations"}, "exchange", "Goldstone"s), "Feynman");
  if (exchange) {
    // Valence-line Green's functions
    num_matrices += num_w * num_kappas * num_blocks;
    // Each thread: two terms, for each partial wave and each multipole that
    // can join it to a core orbital (at most l_max_core + 1)
    num_matrices += std::size_t(omp_get_max_threads()) * 2 * num_components *
                    num_kappas * std::size_t(l_max_core + 1);
    // Stored loop Green's functions: e_a + w and e_a - w
    if (input.get({"Correlations"}, "store_gex", true)) {
      num_matrices += 2 * num_w * num_core * num_kappas * num_blocks;
    }
  }

  return 16.0 * double(subgrid_points * subgrid_points) * double(num_matrices) *
         1.0e-9;
}
