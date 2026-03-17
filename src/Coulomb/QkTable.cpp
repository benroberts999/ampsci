#include "QkTable.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp"
#include <memory>

namespace Coulomb {

void estimate_memory_usage(const std::string &basis_string,
                           const std::string &core_string, int k_cut) {

  std::cout << "Estimating memory usage (nb: may be very rough)\n";

  std::cout << "For basis         : " << basis_string << "\n" << std::flush;
  const auto nk_list = AtomData::n_kappa_list(basis_string);

  const auto core_configs = AtomData::core_parser(core_string);

  // dummy grid:
  const auto radial_grid =
    std::make_shared<const Grid>(GridParameters{10, 1.0e-4, 250.0});

  // build set of H-like orbitals, one n for each kappa up to l=lmax
  std::size_t orbital_count{0};
  std::vector<DiracSpinor> orbs;
  std::vector<DiracSpinor> core;
  for (const auto [n_max, kappa] : nk_list) {
    auto l = Angular::l_k(kappa);
    int n_min = l + 1;
    for (int n = n_min; n <= n_max; ++n) {
      const auto &t =
        orbs.emplace_back(DiracSpinor::exactHlike(n, kappa, radial_grid, 1.0));
      ++orbital_count;
      // if(any c = core_configs.at(i) have c.l()=t.() && c.n()==t.n()){}
      const bool is_core = std::any_of(
        core_configs.begin(), core_configs.end(),
        [&t](const auto &c) { return c.l == t.l() && c.n == t.n(); });
      if (is_core) {
        core.push_back(t);
      }
    }
  }
  assert(orbital_count == orbs.size());

  const auto e_Fermi = core.empty() ? 0.0 : DiracSpinor::max_En(core) + 1.0e-6;

  std::cout << "Total orbitals    : " << orbs.size() << "\n" << std::flush;
  if (!core_string.empty())
    std::cout << "Orbitals in core  : " << core.size() << "\n" << std::flush;

  std::cout << "Counting integrals..." << std::flush;
  QkTable qk;
  const auto [integrals, integrals_Val, integrals_CI, max_k] =
    qk.count_non_zero_integrals(orbs, std::size_t(k_cut), e_Fermi);

  std::cout << "\n";
  std::cout << "With maximum k    : " << max_k << "\n";
  std::cout << "Total integrals   : " << integrals << "\n" << std::flush;
  std::cout << "Excited-only      : " << integrals_Val << "\n" << std::flush;
  std::cout << "Core-Excited      : " << integrals_CI << "\n" << std::flush;

  const auto memory_GB1 = estimate_memory_GB(integrals);
  const auto memory_GB2 = estimate_memory_GB(integrals_Val);
  const auto memory_GB3 = estimate_memory_GB(integrals_CI);

  std::cout << "Estimated Qk size :\n";
  std::cout << "All possible Qk   : " << std::round(memory_GB1)
            << " Gb  [SR]\n";
  std::cout << "Excited-only      : " << std::round(memory_GB2)
            << " Gb  [CI]\n";
  std::cout << "Core-Excited      : " << std::round(memory_GB3)
            << " Gb  (1 or 2 core states allowed)  [MBPT]\n";
}

//==============================================================================

double estimate_memory_GB(std::size_t number_of_integrals, double load_factor) {

  // quite a rough guess

  const auto bucket_count = double(number_of_integrals) / load_factor;

  const auto node_size = sizeof(nk4Index) + // key
                         sizeof(double) +   // value
                         sizeof(void *);    // next pointer (chaining)

  // one pointer per bucket
  const auto bucket_array = bucket_count * double(sizeof(void *));

  return (double(number_of_integrals * node_size) + bucket_array) * 1.0e-9;
}
} // namespace Coulomb