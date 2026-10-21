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
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace CI {

//==============================================================================
ConfigurationSpace::ConfigurationSpace(
  const std::vector<DiracSpinor> &ci_basis, int num_electrons,
  const std::vector<std::string> &leading,
  const std::vector<std::string> &excitations)
  : m_num_electrons(num_electrons),
    m_basis_string(DiracSpinor::state_config(ci_basis)) {

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
    int count = 0;
    for (const auto &shell : configuration) {
      count += shell.num;
    }
    if (count != num_electrons) {
      fmt2::warning();
      fmt::print(": leading configuration '{}' not understood, or not {} "
                 "electrons (form: 6s2.6p); skipped\n",
                 text, num_electrons);
      continue;
    }
    m_leading.push_back(configuration);
  }

  // Default: the lowest configurations of each parity, by orbital energy
  m_default_leading = m_leading.empty();
  if (m_default_leading) {
    std::vector<std::pair<double, Configuration>> even, odd;
    for (auto &configuration : all_configurations()) {
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
    std::vector<Shell> shells;
    for (const auto &Fn : basis_subset(ci_basis, basis_string)) {
      shells.emplace_back(Fn.n(), Fn.l());
    }
    std::sort(shells.begin(), shells.end());
    shells.erase(std::unique(shells.begin(), shells.end()), shells.end());
    m_step_shells.push_back(std::move(shells));
    m_step_strings.push_back(basis_string);
  }

  // With as many full-basis steps as electrons, every configuration can be
  // reached from any leading configuration
  const auto num_full_steps =
    std::count(m_step_shells.cbegin(), m_step_shells.cend(), m_shells);
  m_unrestricted = !m_leading.empty() && num_full_steps >= num_electrons;
}

//==============================================================================
bool ConfigurationSpace::included(
  const std::vector<DiracSpinor::Index> &orbitals) const {
  if (m_unrestricted)
    return true;

  // Both j of a shell count together
  Configuration configuration;
  for (const auto index : orbitals) {
    const auto [n, kappa] = Angular::index_to_nk(index);
    const AtomData::NonRelConfig electron{n, Angular::l_k(kappa), 1};
    const auto shell =
      std::find(configuration.begin(), configuration.end(), electron);
    if (shell == configuration.end()) {
      configuration.push_back(electron);
    } else {
      *shell += electron;
    }
  }
  return included_configuration(configuration);
}

//==============================================================================
bool ConfigurationSpace::included_configuration(
  const Configuration &configuration) const {
  if (m_unrestricted)
    return true;

  for (const auto &leading_configuration : m_leading) {
    // Electrons beyond the leading configuration's occupancy of each shell
    std::vector<Shell> excited;
    for (const auto &shell : configuration) {
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
ConfigurationSpace::all_configurations() const {
  std::vector<Configuration> configurations;
  if (m_shells.empty() || m_num_electrons <= 0)
    return configurations;

  // Each configuration is a non-decreasing list of shell indices, one for
  // each electron; step through them in order
  std::vector<std::size_t> chosen(std::size_t(m_num_electrons), 0);
  while (true) {
    Configuration configuration;
    for (const auto index : chosen) {
      const auto &[n, l] = m_shells[index];
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
    while (position > 0 && chosen[position - 1] + 1 == m_shells.size()) {
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
std::string ConfigurationSpace::summary() const {
  const auto configurations = all_configurations();
  const auto num_included =
    std::count_if(configurations.cbegin(), configurations.cend(),
                  [this](const auto &configuration) {
                    return included_configuration(configuration);
                  });

  std::string out;
  if (m_unrestricted) {
    out += fmt::format("Including all configurations: {} of {}\n", num_included,
                       configurations.size());
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
    out += fmt::format("Configurations: {} of {}\n", num_included,
                       configurations.size());
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
    const auto shell = AtomData::term_parser(term);
    if (!shell.ok() || shell.num < 1)
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
    if (!shell.ok())
      return {};
  }
  std::sort(configuration.begin(), configuration.end());
  return configuration;
}

//==============================================================================
std::string
configuration_string(const ConfigurationSpace::Configuration &configuration) {
  auto text = AtomData::configs_to_string(configuration);
  std::replace(text.begin(), text.end(), ',', '.');
  return text;
}

} // namespace CI
