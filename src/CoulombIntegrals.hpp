#pragma once
#include "DiracSpinor.hpp"
#include "Grid.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include "Wigner_369j.hpp"
#include <vector>

// * Add ability to just update a single integral
// * (all integrals involving given orbital)

/*
Definitions:
y^k_ij(r)   := Int_0^inf [r_min^k/r_max^(k+1)]*rho(f') dr'
rho(r')     := fi(r')*fj(r') + gi(r')gj(r')
Lambda^k_ij := 3js((ji,jj,k),(-1/2,1/2,0))^2 * parity(li+lj+k)

m_C_kakbk "C" (just parity + [j] + 3js, no sign term!)
so, C =
C^k_ij = Sqrt([ji][jj]) * 3js((ji,jj,k),(-1/2,1/2,0)) * parity(li+lj+k)

*/

class Coulomb {

public: // constructor + static functions
  Coulomb(const std::vector<DiracSpinor> &in_core,
          const std::vector<DiracSpinor> &in_valence);
  Coulomb(const Grid &in_grid, const std::vector<DiracSpinor> &in_core,
          const std::vector<DiracSpinor> &in_valence);

  static void calculate_y_ijk(const DiracSpinor &phi_a,
                              const DiracSpinor &phi_b, const int k,
                              std::vector<double> &vabk);

public: // functions
  void form_core_core();
  void form_valence_valence();
  void form_core_valence();

  // MUST calculate values first!
  std::vector<double> calculate_R_abcd_k(const DiracSpinor &psi_a,
                                         const DiracSpinor &psi_b,
                                         const DiracSpinor &psi_c,
                                         const DiracSpinor &psi_d) const;

  // getters
  double get_angular_C_kiakibk(const DiracSpinor &phi_a,
                               const DiracSpinor &phi_b, int k) const;
  double get_angular_L_kiakibk(const DiracSpinor &phi_a,
                               const DiracSpinor &phi_b, int k) const;

  const std::vector<double> &get_y_ijk(const DiracSpinor &phi_i,
                                       const DiracSpinor &phi_j, int k) const;

public: // functions
  void initialise_core_core();

private: // functions
  void initialise_core_valence();
  void initialise_valence_valence();
  void calculate_angular(int ki);

  // write another that returns pair <int, bool> = <index, val?> ?
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
  const std::vector<DiracSpinor> *const c_orbs_ptr;
  const std::vector<DiracSpinor> *const v_orbs_ptr;
  const Grid *const rgrid_ptr;

  std::size_t num_initialised_vv = 0;
  std::size_t num_initialised_vc = 0;
  int m_largest_ki = -1; //-1 not valid, but need 0>x to hold true first time

  // Arrays to store Coulomb integrals. Note: highly non-rectangular
  // Make use of symmety (where appropriate)
  std::vector<std::vector<std::vector<std::vector<double>>>> m_y_abkr;
  std::vector<std::vector<std::vector<std::vector<double>>>> m_y_vckr;
  std::vector<std::vector<std::vector<std::vector<double>>>> m_y_vwkr;
  // Angular coeficients:
  std::vector<std::vector<std::vector<double>>> m_L_kakbk;
  std::vector<std::vector<std::vector<double>>> m_C_kakbk;
};
