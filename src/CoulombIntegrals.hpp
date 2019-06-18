#pragma once
#include "DiracSpinor.hpp"
#include "Grid.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include "Wigner_369j.hpp"
#include <vector>

// * Add ability to just update a single integral
// * (all integrals involving given orbital)

class Coulomb {

public: // constructor + static functions
  Coulomb(const std::vector<DiracSpinor> &in_core,
          const std::vector<DiracSpinor> &in_valence);
  Coulomb(const std::vector<DiracSpinor> &in_orbitals);

  static void calculate_y_ijk(const DiracSpinor &phi_a,
                              const DiracSpinor &phi_b, const int k,
                              std::vector<double> &vabk);

public: // functions
  void calculate_core_core();
  void calculate_valence_valence();
  void calculate_core_valence();

  // NOTE: these take kappa-index! not kappa!
  double get_angular_C_kiakibk(const DiracSpinor &phi_a,
                               const DiracSpinor &phi_b, int k) const;
  double get_angular_L_kiakibk(const DiracSpinor &phi_a,
                               const DiracSpinor &phi_b, int k) const;

  const std::vector<double> &get_y_ijk(const DiracSpinor &phi_i,
                                       const DiracSpinor &phi_j, int k) const;

private: // functions
  void initialise_core_core();
  void initialise_core_valence();
  void initialise_valence_valence();
  void calculate_angular(int ki);

  std::size_t find_valence_index(const DiracSpinor &phi) const;
  std::size_t find_core_index(const DiracSpinor &phi) const;

  const std::vector<std::vector<double>> &get_y_abk(std::size_t a,
                                                    std::size_t b) const;
  const std::vector<std::vector<double>> &get_y_vck(std::size_t a,
                                                    std::size_t b) const;
  const std::vector<std::vector<double>> &get_y_vwk(std::size_t a,
                                                    std::size_t b) const;

  const std::vector<std::vector<double>> &
  get_y_ijk(const DiracSpinor &phi_i, const DiracSpinor &phi_j) const;

  const std::vector<double> &get_angular_C_kiakib_k(int kia, int kib) const;
  const std::vector<double> &get_angular_L_kiakib_k(int kia, int kib) const;

private: // data
  std::size_t num_initialised_vv = 0;
  std::size_t num_initialised_vc = 0;

  const std::vector<DiracSpinor> *const c_orbs_ptr;
  const std::vector<DiracSpinor> *const v_orbs_ptr;
  const Grid *const rgrid_ptr;

  int m_largest_ki = -1; //-1 not valid, but need 0>x to hold true first time

  std::vector<std::vector<std::vector<std::vector<double>>>> m_y_abkr;
  std::vector<std::vector<std::vector<std::vector<double>>>> m_y_vckr;
  std::vector<std::vector<std::vector<std::vector<double>>>> m_y_vwkr;

  std::vector<std::vector<std::vector<double>>> m_angular_L_kakbk;
  std::vector<std::vector<std::vector<double>>> m_angular_C_kakbk;
};
