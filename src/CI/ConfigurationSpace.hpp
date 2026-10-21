#pragma once
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <string>
#include <utility>
#include <vector>

namespace CI {

/*!
  @brief Number of configurations of each parity in the default leading list.
  @details
  More brings in more shells, and so more configurations (and CSFs); fewer
  risks leaving out the dominant configuration of a state of interest. Either
  way, the printed leading list should be checked.
*/
constexpr std::size_t default_leading_per_parity = 5;

/*!
  @brief The configurations included in the CI expansion: leading
  configurations, plus excitations from them.

  @details
  Configurations are non-relativistic (nl): each stands for all of its
  relativistic configurations (both j). The excited electrons of a
  configuration, relative to a leading one, are those beyond the leading
  one's occupancy of each shell.

  Excitations are formed stepwise: each step moves one electron into that
  step's basis, and any of the steps may be used. A configuration is
  included if, for some leading configuration, its excited electrons can each
  be given a different step whose basis contains its shell. E.g., for two
  electrons with leading 4s2 and steps '1:25spdf, 2:16spdf', 25s16s is
  included but 25s2 is not.

  Default leading configurations: the default_leading_per_parity lowest of each
  parity, by the sum of orbital energies. A default-constructed space, or one
  with as many full-basis steps as electrons, includes every configuration.
*/
class ConfigurationSpace {
public:
  //! Non-relativistic configuration: occupancy of each nl shell
  using Configuration = std::vector<AtomData::NonRelConfig>;
  //! Non-relativistic shell: {n, l}
  using Shell = std::pair<int, int>;

private:
  int m_num_electrons{0};
  std::string m_basis_string{};
  std::vector<Shell> m_shells{};
  std::vector<Configuration> m_leading{};
  // True if the leading configurations were generated, not given
  bool m_default_leading{false};
  // Shells, and basis string (empty for ci_basis), of each excitation step
  std::vector<std::vector<Shell>> m_step_shells{};
  std::vector<std::string> m_step_strings{};
  bool m_unrestricted{true};

public:
  //! Includes every configuration
  ConfigurationSpace() = default;

  /*!
    @param ci_basis       The CI single-particle basis.
    @param num_electrons  Number of valence electrons.
    @param leading        Leading configurations, e.g. {"6s2.6p", "6s.6p2"}
                          (see parse_configuration()); empty for the default.
                          Invalid ones are skipped, with a warning.
    @param excitations    Excitation steps: "k", or "k:basis" to restrict step
                          k to that basis (default: ci_basis). The largest k
                          is the number of steps: {"3:6spdf"} is the same as
                          {"1", "2", "3:6spdf"}, and {"0"} gives the leading
                          configurations only.
  */
  ConfigurationSpace(const std::vector<DiracSpinor> &ci_basis,
                     int num_electrons, const std::vector<std::string> &leading,
                     const std::vector<std::string> &excitations);

  //! True if every configuration is included
  bool unrestricted() const { return m_unrestricted; }

  //! True if the configuration of these orbitals (nk_index) is included
  bool included(const std::vector<DiracSpinor::Index> &orbitals) const;

  //! The leading configurations (in the form of the 'leading' input option),
  //! the excitation steps, and the number of configurations included; or one
  //! line, if every configuration is included. Identifies the space.
  std::string summary() const;

private:
  bool included_configuration(const Configuration &configuration) const;

  // Every configuration of the electrons in the CI basis shells
  std::vector<Configuration> all_configurations() const;
};

/*!
  @brief Parses a configuration, e.g. "6s2.6p": shells separated by '.', each
  with its number of electrons (default 1).
  @return Empty if not a valid configuration.
*/
ConfigurationSpace::Configuration parse_configuration(const std::string &text);

//! Configuration as text, e.g. "6s2.6p1"
std::string
configuration_string(const ConfigurationSpace::Configuration &configuration);

} // namespace CI
