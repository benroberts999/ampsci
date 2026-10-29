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

  With hole CSFs (see CSF2p1h), the configurations with a hole (one more
  particle than electrons, and a hole in one of the hole shells; written with
  a negative occupancy, e.g. 6s-1.6p2) are included too. The hole is created
  by an excitation from a leading configuration without one (6s-1.6p2 and
  6s-1.6p.7p are single excitations from 6p), or is that of the leading
  configuration; it is never moved or filled, since that would put an
  electron into the core. There is at most one hole. The default leading
  configurations have no holes. With a hole-free leading configuration and as
  many full-basis steps as particles, every configuration is included.
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
  // Shells that may hold a hole, and shells of the particles of the
  // configurations with a hole (empty if no holes)
  std::vector<Shell> m_hole_shells{};
  std::vector<Shell> m_hole_basis_shells{};

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
    @param holes          Core orbitals that may hold a hole (see CSF2p1h);
                          none by default.
    @param hole_basis     Basis of the particles of the configurations with a
                          hole; @p ci_basis if empty.
  */
  ConfigurationSpace(const std::vector<DiracSpinor> &ci_basis,
                     int num_electrons, const std::vector<std::string> &leading,
                     const std::vector<std::string> &excitations,
                     const std::vector<DiracSpinor> &holes = {},
                     const std::vector<DiracSpinor> &hole_basis = {});

  //! True if every configuration is included
  bool unrestricted() const { return m_unrestricted; }

  //! True if the configuration of these orbitals (nk_index), with holes in
  //! these core orbitals, is included
  bool included(const std::vector<DiracSpinor::Index> &orbitals,
                const std::vector<DiracSpinor::Index> &holes = {}) const;

  //! The leading configurations (in the form of the 'leading' input option),
  //! the excitation steps, and the number of configurations included; or one
  //! line, if every configuration is included. Identifies the space.
  std::string summary() const;

private:
  bool included_configuration(const Configuration &configuration) const;

  // Every configuration of num_particles particles in the given shells
  std::vector<Configuration>
  particle_configurations(const std::vector<Shell> &shells,
                          int num_particles) const;
  // Every configuration of the electrons in the CI basis shells
  std::vector<Configuration> electron_configurations() const;
  // Every configuration with a hole: one hole shell, and one more particle
  // than electrons in the hole basis shells
  std::vector<Configuration> hole_configurations() const;
  // Both of the above
  std::vector<Configuration> all_configurations() const;
};

//! Non-relativistic shells {n, l} of the orbitals: sorted, each once
std::vector<ConfigurationSpace::Shell>
shells_of(const std::vector<DiracSpinor> &orbitals);

//! The holes of a configuration (shells with negative occupancy), one entry
//! per hole
std::vector<ConfigurationSpace::Shell>
configuration_holes(const ConfigurationSpace::Configuration &configuration);

/*!
  @brief Parses a configuration, e.g. "6s2.6p": shells separated by '.', each
  with its number of electrons (default 1); negative for a hole, e.g.
  "6s-1.6p2".
  @return Empty if not a valid configuration.
*/
ConfigurationSpace::Configuration parse_configuration(const std::string &text);

//! Configuration as text, e.g. "6s2.6p1" (a hole as "6s-1")
std::string
configuration_string(const ConfigurationSpace::Configuration &configuration);

} // namespace CI
