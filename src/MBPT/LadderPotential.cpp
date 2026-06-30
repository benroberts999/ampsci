#include "MBPT/LadderPotential.hpp"
#include "Angular/include.hpp"
#include "CI/CI_Integrals.hpp"
#include "Coulomb/include.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/Ladder.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace MBPT {

//==============================================================================
std::optional<LadderOptions>
parse_ladder_options(const IO::InputBlock &correlations,
                     const std::string &identity) {
  using namespace std::string_literals;

  auto ladder_block = correlations.getBlock("Ladder");

  // If the parent (Correlations) is in help mode, propagate 'help' into the
  // ladder block so its check() documents the options under
  // `ampsci -i Correlations` (creating the block if it isn't present).
  if (correlations.has_option("help")) {
    if (!ladder_block) {
      ladder_block = IO::InputBlock{"Ladder"};
    }
    ladder_block->add("help;");
  }

  if (!ladder_block)
    return std::nullopt;

  ladder_block->check(
    {{"min_n_core", "lowest core n to include [1]"},
     {"basis", "Basis string specifying which states to include; must be a "
               "subset of full "
               "basis), e.g. '30spdf'. Default: include entire excited basis."},
     {"max_k", "maximum k to include in Qk. Put -1 to include all. [8]"},
     {"include_L4", "Inlcude 4th Ladder diagram [false]"},
     {"full_basis",
      "Construct the ladder correlation potential Sigma_L = sum_i L|i><i| "
      "using the full basis for {|i>}, otherwise just single HF |v> "
      "eigenstate. Using the full basis usually requires extending Qk. "
      "[false]"},
     {"Qk_file", "Filename for storing Qk Coulomb integrals. By default, is "
                 "<Identity>.qk. If 'false' will not read or write."},
     {"Lk_file", "Filename for storing Lk ladder integrals. By default, is "
                 "<Identity>.lk. If 'false' will not read or write."},
     {"from_scratch", "If true, don't read existing Qk/Lk files (still "
                      "writes). [false]"},
     {"max_it", "Max # iterations. If zero, will simply read ladder diagrams "
                "in (and then run a symmetry test). [15]"},
     {"damp",
      "Damping factor for iterations, [0,1). 0 means no damping. [0.0]"},
     {"eps_target", "Target for convergance [1.0e-5]"}});

  if (ladder_block->has_option("help"))
    return std::nullopt;

  const auto L4 = ladder_block->get("include_L4", false);
  const auto lk4 = L4 ? "_l4" : "";
  return LadderOptions{
    ladder_block->get("basis", ""s),
    ladder_block->get("min_n_core", 1),
    ladder_block->get("max_k", 8),
    L4,
    ladder_block->get("full_basis", false),
    ladder_block->get("max_it", 15),
    ladder_block->get("damp", 0.0),
    ladder_block->get("eps_target", 1.0e-5),
    ladder_block->get("Qk_file", identity + ".qk.abf"s),
    ladder_block->get("Lk_file", identity + lk4 + ".lk.abf"s),
    ladder_block->get("from_scratch", false)};
}

//==============================================================================
LadderPotential::LadderPotential(const Wavefunction &wf, double r0, double rmax,
                                 std::size_t stride,
                                 const LadderOptions &options)
  : m_r0(r0),
    m_rmax(rmax),
    m_stride(stride),
    m_include_L4(options.include_L4),
    m_full_basis(options.full_basis) {

  const auto min_n_core = options.n_min_core;
  const auto &basis_str = options.basis_str;
  const auto max_k = options.max_k;
  const auto include_L4 = options.include_L4;
  const auto full_basis = options.full_basis;
  const auto max_it = options.max_it;
  const auto a_damp = options.damp;
  const auto eps_target = options.eps_target;
  const auto &Qk_file = options.Qk_file;
  const auto &Lk_file = options.Lk_file;
  const auto from_scratch = options.from_scratch;

  // Sort basis into core/excited/valence
  const auto en_core = wf.FermiLevel();
  const auto holes = qip::select_if(wf.basis(), [&](const auto &Fn) {
    return Fn.en() < en_core && Fn.n() >= min_n_core;
  });
  const auto excited =
    CI::basis_subset(wf.basis(), basis_str, wf.coreConfiguration());

  // nb: use _basis_ version of valence states in iterations
  std::vector<DiracSpinor> valence;
  for (const auto &Fv : wf.valence()) {
    const auto pFv = std::find(wf.basis().cbegin(), wf.basis().cend(), Fv);
    if (pFv == wf.basis().cend()) {
      std::cout << "Warning: Basis missing valence state: " << Fv << "\n";
      continue;
    }
    valence.push_back(*pFv);
  }

  std::cout << "\n";
  std::cout << "basis        = " << DiracSpinor::state_config(excited) << "\n";
  std::cout << "min_n (core) = " << min_n_core << "\n";
  std::cout << std::boolalpha;
  std::cout << "include_L4   = " << include_L4 << "\n";
  std::cout << "full_basis   = " << full_basis << "\n";
  std::cout << "max_k        = " << max_k << "\n";
  std::cout << "max_it       = " << max_it << "\n";
  std::cout << "damp         = " << a_damp << "\n";
  std::cout << "eps_target   = " << eps_target << "\n";

  // Create the "total" basis, which has core+excited, but only those states
  // actually included (i.e., [n_min, n_max]). This is used to calculate Qk.
  // Reduces size of Qk; should also reduce lookup time.
  const auto both = qip::merge(holes, excited);
  // the "i" orbitals: core + valence
  const auto core_and_val = qip::merge(holes, valence);

  // Fill Yk table (used to fill Qk table)
  const Coulomb::YkTable yk(both);
  // Store reference to Yk's 6J table (save typing):
  const auto &sjt = yk.SixJ();

  // Form the Qk table.
  Coulomb::QkTable qk;
  std::cout << "\nFill (or read) Qk table:\n";
  if (!from_scratch)
    qk.read(Qk_file);
  qk.fill(both, yk, max_k);
  qk.write(Qk_file);

  std::cout << "\nMBPT(2) energy shifts, using HF vs. spline legs" << std::endl;
  fmt::print("{:<5s} {:>13s} {:>13s} {:>13s}   {}\n", "state", "de(HF)",
             "de(basis)", "de(Qk)", "eps");
  const auto dec_1 = MBPT::de_core(yk, yk, wf.core(), excited);
  const auto dec_2 = MBPT::de_core(yk, yk, holes, excited);
  const auto dec_3 = MBPT::de_core(qk, qk, holes, excited);
  const auto dec_eps = std::abs(2.0 * (dec_1 - dec_2) / (dec_1 + dec_2));
  fmt::print("{:<5s} {:+13.7f} {:+13.7f} {:+13.7f}   {:.0e}\n", "core", dec_1,
             dec_2, dec_3, dec_eps);
  for (const auto &sv : valence) {
    const auto pv = wf.getState(sv.n(), sv.kappa());
    assert(pv != nullptr); // valence is subset of wf.valence()
    const auto de_1 = MBPT::de_valence(*pv, yk, yk, holes, excited);
    const auto de_2 = MBPT::de_valence(sv, yk, yk, holes, excited);
    const auto de_3 = MBPT::de_valence(sv, qk, qk, holes, excited);
    const auto de_eps = std::abs(2.0 * (de_1 - de_2) / (de_1 + de_2));
    fmt::print("{:<5s} {:+13.7f} {:+13.7f} {:+13.7f}   {:.0e}\n",
               sv.shortSymbol(), de_1, de_2, de_3, de_eps);
  }
  std::cout << "\n";

  // Fill Lk table:
  Coulomb::LkTable lk;
  Coulomb::LkTable lk_next;
  // Each run will continue from where last left off
  const bool did_read_Lk = !from_scratch ? lk_next.read(Lk_file) : false;
  std::cout << (did_read_Lk ? "\nRe-starting using existing ladder diagrams\n" :
                              "\nCalculating ladder diagrams from scratch\n");

  //----------------------------------------------------------------------------
  // Iterate Lk equations

  // "First iteration": fill_Lk_mnib keeps existing (read) entries and adds only
  // newly-required diagrams (or fills from scratch if not read). Subsequent
  // iterations below re-iterate the existing table.
  std::cout << "\nInitial fill of Lk table: core + valence\n" << std::flush;
  MBPT::fill_Lk_mnib(&lk_next, qk, excited, holes, core_and_val, include_L4,
                     sjt, max_k, true);
  lk_next.write(Lk_file);
  lk = lk_next;

  // convert to inverse cm
  const auto icm = PhysConst::Hartree_invcm;

  // initial corrections
  double de_c0 = MBPT::de_core(qk, lk, holes, excited);
  std::vector<double> de_0;
  std::cout << "\nInitial Ladder Corrections:\n";
  fmt::print("de_l({:4}): {:+12.9f} au  {:+12.5f} cm^-1\n", "core", de_c0,
             de_c0 * icm);
  for (const auto &v : valence) {
    const auto de_v = MBPT::de_valence(v, qk, lk, holes, excited);
    de_0.push_back(de_v);
    fmt::print("de_l({:>4}): {:+12.9f} au  {:+12.5f} cm^-1\n", v.shortSymbol(),
               de_v, de_v * icm);
  }

  // Iterate for core states:
  for (int it = 1; it <= max_it; ++it) {
    std::cout << "\ncore it:" << it << "\n";
    MBPT::update_Lk_mnib(&lk_next, qk, excited, holes, holes, include_L4, sjt,
                         &lk, a_damp, true);
    // for the damping, don't swap; lk and lk_prev must match before next iter
    lk = lk_next;
    lk.write(Lk_file);
    const auto de_c = MBPT::de_core(qk, lk, holes, excited);
    const auto eps_c = std::abs(2.0 * (de_c - de_c0) / (de_c + de_c0));
    de_c0 = de_c;
    fmt::print("de_l({:>4}): {:+12.9f} au  {:+12.5f} cm^-1  {:.1e}\n", "core",
               de_c, de_c * icm, eps_c);
    if (eps_c < eps_target)
      break;
  }

  // Print energy corrections after core convergence; refresh de_0 for valence
  std::cout << "\nAfter core convergence:\n";
  fmt::print("de_l({:4}): {:+12.9f} au  {:+12.5f} cm^-1\n", "core", de_c0,
             de_c0 * icm);
  for (std::size_t i = 0; i < valence.size(); ++i) {
    const auto de_v = MBPT::de_valence(valence[i], qk, lk, holes, excited);
    de_0[i] = de_v;
    fmt::print("de_l({:>4}): {:+12.9f} au  {:+12.5f} cm^-1\n",
               valence[i].shortSymbol(), de_v, de_v * icm);
  }

  // Iterate for each required valence state; drop converged states each iter
  std::vector<DiracSpinor> unconverged_valence = valence;
  for (int it = 1; it <= max_it; ++it) {
    std::cout << "\nvalence it:" << it << "\n";
    MBPT::update_Lk_mnib(&lk_next, qk, excited, holes, unconverged_valence,
                         include_L4, sjt, &lk, a_damp, true);
    lk = lk_next;
    lk.write(Lk_file);
    for (std::size_t i = 0; i < valence.size(); ++i) {
      const auto &v = valence[i];
      const auto de_v = MBPT::de_valence(v, qk, lk, holes, excited);
      const auto eps_v = std::abs(2.0 * (de_v - de_0[i]) / (de_v + de_0[i]));
      fmt::print("de_l({:>4}): {:+12.9f} au  {:+12.5f} cm^-1  {:.1e}\n",
                 v.shortSymbol(), de_v, de_v * icm, eps_v);
      auto it2 =
        std::find(unconverged_valence.begin(), unconverged_valence.end(), v);
      if (it2 != unconverged_valence.end()) {
        de_0[i] = de_v;
        if (eps_v < eps_target) {
          // if converged: remove from list (speed up convergance)
          unconverged_valence.erase(it2);
        }
      }
    }
    if (unconverged_valence.empty())
      break;
  }
  std::cout << "\n";

  //----------------------------------------------------------------------------

  // Extend Qk table to the FULL basis. With full_basis, Sigma_ladder projects
  // the ladder onto the full basis of each kappa_v, so Lkmnij() needs Qk
  // integrals involving those (high-n) projection states - which are absent if
  // Qk was filled only over the [n_min,n_max] subset 'both'. Not always
  // strictly required (only when projection states fall outside 'both'), but we
  // just refill over the full basis to be safe. Skipped for single-state |v>.
  if (full_basis) {
    std::cout << "\nExtending Qk table to full basis (for Sigma_ladder):\n"
              << std::flush;
    const Coulomb::YkTable yk_full(wf.basis());
    qk.fill(wf.basis(), yk_full, max_k);
    qk.write(Qk_file);
  }

  // Store for later Sigma_L construction (see Sigma_ladder):
  m_holes = holes;
  m_excited = excited;
  m_valence = valence;
  m_qk = std::move(qk);
  m_lk = std::move(lk);
  m_sjt = sjt;
  if (full_basis)
    m_proj_basis = wf.basis();
}

//==============================================================================
GMatrix LadderPotential::Sigma_ladder(int kappa, double ev,
                                      const DiracSpinor *Fv,
                                      bool include_G) const {
  // Projection set (Sigma_ladder filters by kappa internally): the full basis
  // if full_basis was set, otherwise the single valence state |v>.
  const auto single = Fv ? std::vector<DiracSpinor>{*Fv} : m_excited;
  const auto &proj = m_full_basis ? m_proj_basis : single;
  return MBPT::Sigma_ladder(kappa, ev, m_holes, m_excited, proj, m_qk, &m_lk,
                            m_sjt, m_include_L4, m_r0, m_rmax, m_stride,
                            include_G);
}

} // namespace MBPT
