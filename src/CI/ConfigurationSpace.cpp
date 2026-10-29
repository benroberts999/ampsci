#include "ConfigurationSpace.hpp"
#include "Angular/Wigner369j.hpp"
#include "CI_Integrals.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/color.hpp"
#include "fmt/format.hpp"
#include "qip/String.hpp"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace CI {

//==============================================================================
ConfigurationSpace::ConfigurationSpace(
  const std::vector<DiracSpinor> &ci_basis, int num_electrons,
  const std::vector<std::string> &leading,
  const std::vector<std::string> &excitations,
  const std::vector<DiracSpinor> &holes,
  const std::vector<DiracSpinor> &hole_basis)
  : m_num_electrons(num_electrons),
    m_basis_string(DiracSpinor::state_config(ci_basis)),
    m_hole_shells(shells_of(holes)),
    m_hole_basis_shells(
      holes.empty() ? std::vector<Shell>{} :
                      shells_of(hole_basis.empty() ? ci_basis : hole_basis)) {

  // Shells of the CI basis, with energies averaged over both j (weight 2j+1):
  // shell -> {sum of (2j+1) * energy, sum of (2j+1)}
  std::map<Shell, std::pair<double, double>> energy_sums;
  for (const auto &Fn : ci_basis) {
    auto &[weighted_energy, weight] = energy_sums[{Fn.n(), Fn.l()}];
    weighted_energy += (Fn.twoj() + 1) * Fn.en();
    weight += Fn.twoj() + 1;
  }
  for (const auto &[shell, sums] : energy_sums) {
    m_shells.push_back(shell);
  }

  for (const auto &text : leading) {
    const auto configuration = parse_configuration(IO::removeSpaces(text));
    // A hole counts as minus one electron
    int count = 0;
    for (const auto &shell : configuration) {
      count += shell.num;
    }
    if (configuration.empty() || count != num_electrons) {
      fmt2::warning();
      fmt::print(": leading configuration '{}' not understood, or not {} "
                 "electrons (form: 6s2.6p, or 6s-1.6p2 with a hole); "
                 "skipped\n",
                 text, num_electrons);
      continue;
    }
    // At most one hole, in one of the hole shells
    const auto configuration_hole_shells = configuration_holes(configuration);
    const auto hole_ok =
      configuration_hole_shells.empty() ||
      (configuration_hole_shells.size() == 1 &&
       std::find(m_hole_shells.cbegin(), m_hole_shells.cend(),
                 configuration_hole_shells.front()) != m_hole_shells.cend());
    if (!hole_ok) {
      fmt2::warning();
      fmt::print(": leading configuration '{}': at most one hole, in a shell "
                 "of 'holes'; skipped\n",
                 text);
      continue;
    }
    m_leading.push_back(configuration);
  }

  // Default: the lowest configurations of each parity, by orbital energy;
  // none with a hole
  m_default_leading = m_leading.empty();
  if (m_default_leading) {
    std::vector<std::pair<double, Configuration>> even, odd;
    for (auto &configuration : electron_configurations()) {
      double energy = 0.0;
      int sum_l = 0;
      for (const auto &shell : configuration) {
        const auto &[weighted_energy, weight] =
          energy_sums.at({shell.n, shell.l});
        energy += shell.num * weighted_energy / weight;
        sum_l += shell.num * shell.l;
      }
      auto &same_parity = sum_l % 2 == 0 ? even : odd;
      same_parity.emplace_back(energy, std::move(configuration));
    }
    for (auto *same_parity : {&even, &odd}) {
      std::stable_sort(
        same_parity->begin(), same_parity->end(),
        [](const auto &a, const auto &b) { return a.first < b.first; });
      const auto num =
        std::min(default_leading_per_parity, same_parity->size());
      for (std::size_t i = 0; i < num; ++i) {
        m_leading.push_back((*same_parity)[i].second);
      }
    }
  }

  // Excitation steps: "k" or "k:basis"; the largest k is the number of steps
  std::map<int, std::string> step_basis_strings;
  int num_steps = 0;
  for (const auto &entry : excitations) {
    // Spaces are kept in quoted options
    const auto text = IO::removeSpaces(entry);
    const auto colon = text.find(':');
    const auto step = text.substr(0, colon);
    if (!qip::string_is_integer(step)) {
      fmt2::warning();
      fmt::print(": excitation '{}' not understood (form: 2 or 2:16spdf); "
                 "skipped\n",
                 text);
      continue;
    }
    num_steps = std::max(num_steps, std::stoi(step));
    if (colon != std::string::npos) {
      step_basis_strings[std::stoi(step)] = text.substr(colon + 1);
    }
  }
  for (int k = 1; k <= num_steps; ++k) {
    const auto &basis_string = step_basis_strings[k];
    m_step_shells.push_back(shells_of(basis_subset(ci_basis, basis_string)));
    m_step_strings.push_back(basis_string);
  }

  // With as many full-basis steps as particles, every configuration can be
  // reached from a leading configuration without a hole (the configurations
  // with a hole have one more particle than electrons)
  const auto num_full_steps =
    std::count(m_step_shells.cbegin(), m_step_shells.cend(), m_shells);
  const auto num_particles =
    m_hole_shells.empty() ? num_electrons : num_electrons + 1;
  const auto hole_free_leading = std::any_of(
    m_leading.cbegin(), m_leading.cend(), [](const auto &configuration) {
      return configuration_holes(configuration).empty();
    });
  m_unrestricted = hole_free_leading && num_full_steps >= num_particles;
}

//==============================================================================
bool ConfigurationSpace::included(
  const std::vector<DiracSpinor::Index> &orbitals,
  const std::vector<DiracSpinor::Index> &holes) const {
  if (m_unrestricted)
    return true;

  // Both j of a shell count together; a hole counts as minus one electron
  Configuration configuration;
  for (const auto *list : {&orbitals, &holes}) {
    const int num = list == &orbitals ? 1 : -1;
    for (const auto index : *list) {
      const auto [n, kappa] = Angular::index_to_nk(index);
      const AtomData::NonRelConfig electron{n, Angular::l_k(kappa), num};
      const auto shell =
        std::find(configuration.begin(), configuration.end(), electron);
      if (shell == configuration.end()) {
        configuration.push_back(electron);
      } else {
        *shell += electron;
      }
    }
  }
  return included_configuration(configuration);
}

//==============================================================================
bool ConfigurationSpace::included_configuration(
  const Configuration &configuration) const {
  if (m_unrestricted)
    return true;

  const auto holes = configuration_holes(configuration);
  for (const auto &leading_configuration : m_leading) {
    // The hole is created by an excitation from a leading configuration
    // without one, or is that of the leading configuration: it is never
    // moved or filled (that would put an electron into the core)
    const auto leading_holes = configuration_holes(leading_configuration);
    if (!leading_holes.empty() && leading_holes != holes)
      continue;

    // Electrons beyond the leading configuration's occupancy of each shell
    // (the particles; the hole is not an excited electron)
    std::vector<Shell> excited;
    for (const auto &shell : configuration) {
      if (shell.num < 0)
        continue;
      const auto in_leading = std::find(leading_configuration.cbegin(),
                                        leading_configuration.cend(), shell);
      const auto num_leading =
        in_leading == leading_configuration.cend() ? 0 : in_leading->num;
      for (int i = num_leading; i < shell.num; ++i) {
        excited.emplace_back(shell.n, shell.l);
      }
    }
    if (excited.size() > m_step_shells.size())
      continue;

    // Each excited electron needs a step of its own, whose basis contains its
    // shell: try the steps in every order
    std::vector<std::size_t> steps(m_step_shells.size());
    std::iota(steps.begin(), steps.end(), 0);
    do {
      std::size_t num_placed = 0;
      while (num_placed < excited.size() &&
             std::binary_search(m_step_shells[steps[num_placed]].cbegin(),
                                m_step_shells[steps[num_placed]].cend(),
                                excited[num_placed])) {
        ++num_placed;
      }
      if (num_placed == excited.size())
        return true;
    } while (std::next_permutation(steps.begin(), steps.end()));
  }
  return false;
}

//==============================================================================
std::vector<ConfigurationSpace::Configuration>
ConfigurationSpace::particle_configurations(const std::vector<Shell> &shells,
                                            int num_particles) const {
  std::vector<Configuration> configurations;
  if (shells.empty() || num_particles <= 0)
    return configurations;

  // Each configuration is a non-decreasing list of shell indices, one for
  // each particle; step through them in order
  std::vector<std::size_t> chosen(std::size_t(num_particles), 0);
  while (true) {
    Configuration configuration;
    for (const auto index : chosen) {
      const auto &[n, l] = shells[index];
      if (!configuration.empty() && configuration.back().n == n &&
          configuration.back().l == l) {
        ++configuration.back().num;
      } else {
        configuration.emplace_back(n, l, 1);
      }
    }
    if (std::all_of(configuration.cbegin(), configuration.cend(),
                    [](const auto &shell) { return shell.ok(); })) {
      configurations.push_back(std::move(configuration));
    }

    // Next list: increase the last index that is not at the last shell, and
    // set all after it equal to it
    auto position = chosen.size();
    while (position > 0 && chosen[position - 1] + 1 == shells.size()) {
      --position;
    }
    if (position == 0)
      break;
    std::fill(chosen.begin() + std::ptrdiff_t(position - 1), chosen.end(),
              chosen[position - 1] + 1);
  }
  return configurations;
}

//==============================================================================
std::vector<ConfigurationSpace::Configuration>
ConfigurationSpace::electron_configurations() const {
  return particle_configurations(m_shells, m_num_electrons);
}

//==============================================================================
std::vector<ConfigurationSpace::Configuration>
ConfigurationSpace::hole_configurations() const {
  std::vector<Configuration> configurations;
  if (m_hole_shells.empty())
    return configurations;

  // One hole, and one more particle than electrons
  const auto particles =
    particle_configurations(m_hole_basis_shells, m_num_electrons + 1);
  for (const auto &[n, l] : m_hole_shells) {
    for (auto configuration : particles) {
      configuration.emplace_back(n, l, -1);
      std::sort(configuration.begin(), configuration.end());
      configurations.push_back(std::move(configuration));
    }
  }
  return configurations;
}

//==============================================================================
std::vector<ConfigurationSpace::Configuration>
ConfigurationSpace::all_configurations() const {
  auto configurations = electron_configurations();
  for (auto &configuration : hole_configurations()) {
    configurations.push_back(std::move(configuration));
  }
  return configurations;
}

//==============================================================================
std::string ConfigurationSpace::summary() const {
  const auto configurations = all_configurations();
  std::size_t num_included = 0;
  std::size_t num_holes = 0;
  std::size_t num_holes_included = 0;
  for (const auto &configuration : configurations) {
    const auto included = included_configuration(configuration);
    const auto has_hole = !configuration_holes(configuration).empty();
    num_included += included ? 1 : 0;
    num_holes += has_hole ? 1 : 0;
    num_holes_included += included && has_hole ? 1 : 0;
  }
  // Of which, those with a hole (if there are hole CSFs)
  std::string hole_text;
  if (!m_hole_shells.empty()) {
    hole_text = m_unrestricted ? fmt::format(" ({} with a hole)", num_holes) :
                                 fmt::format(" ({} of {} with a hole)",
                                             num_holes_included, num_holes);
  }

  std::string out;
  if (m_unrestricted) {
    out += fmt::format("Including all configurations: {} of {}{}\n",
                       num_included, configurations.size(), hole_text);
  } else {
    // In the form of the 'leading' input option, to be copied and edited;
    // continuation lines indented under the list
    std::string leading_text;
    for (std::size_t i = 0; i < m_leading.size(); ++i) {
      leading_text += (i == 0 ? "" : ", ") + configuration_string(m_leading[i]);
    }
    auto leading_lines = qip::wrap(leading_text, 80, std::string(9, ' '));
    leading_lines.replace(0, 9, "Leading: ");
    out += leading_lines + "\n";

    std::string steps_text;
    for (std::size_t step = 0; step < m_step_strings.size(); ++step) {
      steps_text +=
        (step == 0 ? "" : ", ") + std::to_string(step + 1) + " with " +
        (m_step_strings[step].empty() ? m_basis_string : m_step_strings[step]);
    }
    out += fmt::format("Including excitations: {}\n",
                       steps_text.empty() ? "none" : steps_text);
    out += fmt::format("Configurations: {} of {}{}\n", num_included,
                       configurations.size(), hole_text);
    if (m_default_leading) {
      out += "* NOTE: leading configs auto generated: check they are OK for "
             "your problem!\n";
    }
  }

  return out;
}

//==============================================================================
ConfigurationSpace::Configuration parse_configuration(const std::string &text) {
  ConfigurationSpace::Configuration configuration;
  for (auto term : qip::split(text, '.')) {
    // The number of electrons may be left off: one
    if (!term.empty() &&
        std::isalpha(static_cast<unsigned char>(term.back()))) {
      term += "1";
    }
    // n, the l symbol, then the number of electrons: negative for a hole
    const auto symbol = term.find_first_not_of("0123456789");
    if (symbol == std::string::npos ||
        !qip::string_is_integer(term.substr(symbol + 1))) {
      return {};
    }
    const auto shell = AtomData::term_parser(term);
    const AtomData::NonRelConfig filled{shell.n, shell.l, std::abs(shell.num)};
    if (shell.num == 0 || !filled.ok())
      return {};
    const auto existing =
      std::find(configuration.begin(), configuration.end(), shell);
    if (existing == configuration.end()) {
      configuration.push_back(shell);
    } else {
      *existing += shell;
    }
  }
  for (const auto &shell : configuration) {
    const AtomData::NonRelConfig filled{shell.n, shell.l, std::abs(shell.num)};
    if (shell.num == 0 || !filled.ok())
      return {};
  }
  std::sort(configuration.begin(), configuration.end());
  return configuration;
}

//==============================================================================
std::vector<ConfigurationSpace::Shell>
shells_of(const std::vector<DiracSpinor> &orbitals) {
  std::vector<ConfigurationSpace::Shell> shells;
  for (const auto &Fn : orbitals) {
    shells.emplace_back(Fn.n(), Fn.l());
  }
  std::sort(shells.begin(), shells.end());
  shells.erase(std::unique(shells.begin(), shells.end()), shells.end());
  return shells;
}

//==============================================================================
std::vector<ConfigurationSpace::Shell>
configuration_holes(const ConfigurationSpace::Configuration &configuration) {
  std::vector<ConfigurationSpace::Shell> holes;
  for (const auto &shell : configuration) {
    for (int i = shell.num; i < 0; ++i) {
      holes.emplace_back(shell.n, shell.l);
    }
  }
  return holes;
}

//==============================================================================
std::string
configuration_string(const ConfigurationSpace::Configuration &configuration) {
  auto text = AtomData::configs_to_string(configuration);
  std::replace(text.begin(), text.end(), ',', '.');
  return text;
}

} // namespace CI
