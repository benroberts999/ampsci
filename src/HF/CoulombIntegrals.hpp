#pragma once
#include <tuple>
#include <vector>
class DiracSpinor;
class Grid;

// * Add ability to just update a single integral
// * (all integrals involving given orbital)

//! @brief Calculates + stores Coulomb integrals, with \f$y_{ab}^k\f$ look-up
/*! @details
  - Needs a good clean-up, lots of dead code etc.

Definitions:
  - y^k_ij(r)   := Int_0^inf [r_min^k/r_max^(k+1)]*rho(f') dr'
  - rho(r')     := fi(r')*fj(r') + gi(r')gj(r')
  - Lambda^k_ij := 3js((ji,jj,k),(-1/2,1/2,0))^2 * parity(li+lj+k)

  - m_C_kakbk "C" (just parity + [j] + 3js, no sign term!)
  - so, C =
C^k_ij = sqrt([ji][jj]) * 3js((ji,jj,k),(-1/2,1/2,0)) * parity(li+lj+k)
*/
class Coulomb {

public: // constructor + static functions
  Coulomb(const std::vector<DiracSpinor> &in_core,
          const std::vector<DiracSpinor> &in_valence);
  Coulomb(const Grid &in_grid, const std::vector<DiracSpinor> &in_core,
          const std::vector<DiracSpinor> &in_valence);

  Coulomb &operator=(const Coulomb &) = delete; // copy assignment
  Coulomb(const Coulomb &) = default;           // copy constructor
  ~Coulomb() = default;

  //! Calculates Hartree Screening functions \f$y^k_{ab}(r)\f$
  static void calculate_y_ijk(const DiracSpinor &Fa, const DiracSpinor &Fb,
                              const int k, std::vector<double> &vabk,
                              const std::size_t maxi = 0);

public:
  //! Forms/updtaes yk_ab matrix for all core a and b
  void form_core_core();
  //! Updates yk_ab matrix: only terms involving Fa (Fa in core)
  void form_core_core(const DiracSpinor &Fa);
  //! Forms/updtaes yk_vw matrix for all valence v and w
  void form_valence_valence();
  //! Forms/updtaes yk_va matrix for all valence v and core a
  void form_core_valence();
  //! Forms/updtaes yk_va matrix, only terms involving Fn
  void form_core_valence(const DiracSpinor &Fn);

  // MUST calculate values first!
  //! Calculates R^k_abcd for all k. y^k_ab must already exist!
  std::vector<double> calculate_R_abcd_k(const DiracSpinor &Fa,
                                         const DiracSpinor &Fb,
                                         const DiracSpinor &Fc,
                                         const DiracSpinor &Fd) const;

  //! Calculates R^k_abcd for given k. From scratch (calculates y)
  static double Rk_abcd_any(const DiracSpinor &Fa, const DiracSpinor &Fb,
                            const DiracSpinor &Fc, const DiracSpinor &Fd,
                            const int k);
  //! Calculates Q^k_abcd for given k. From scratch (calculates y)
  static double Qk_abcd_any(const DiracSpinor &Fa, const DiracSpinor &Fb,
                            const DiracSpinor &Fc, const DiracSpinor &Fd,
                            const int k);
  //! Calculates Z^k_abcd for given k. From scratch (calculates y)
  static double Zk_abcd_any(const DiracSpinor &Fa, const DiracSpinor &Fb,
                            const DiracSpinor &Fc, const DiracSpinor &Fd,
                            const int k);

  static DiracSpinor Qk_abcd_rhs(const DiracSpinor &Fa, const DiracSpinor &Fb,
                                 const DiracSpinor &Fc, const DiracSpinor &Fd,
                                 const int k);
  static DiracSpinor Rk_abcd_rhs(const DiracSpinor &Fa, const DiracSpinor &Fb,
                                 const DiracSpinor &Fc, const DiracSpinor &Fd,
                                 const int k);

  // getters

  //! Look-up for angular factor. This should be replaced with Angular version!
  const std::vector<double> &get_angular_C_kiakib_k(int kia, int kib) const;
  //! Look-up for angular factor. This should be replaced with Angular version!
  const std::vector<double> &get_angular_L_kiakib_k(int kia, int kib) const;
  //! Look-up for y^k_ij (core or valence, given k); must exist already
  const std::vector<double> &get_y_ijk(const DiracSpinor &Fi,
                                       const DiracSpinor &Fj, int k) const;
  //! Look-up for y^k_ij (core or valence). Returns all k. Must exist already
  const std::vector<std::vector<double>> &
  get_y_ijk(const DiracSpinor &Fi, const DiracSpinor &Fj) const;

public:
  //! Sizes the relvant y and angular matrix for core-core integrals
  void initialise_core_core();
  //! Sizes the relvant y and angular matrix for core-valence integrals
  void initialise_core_valence();

private: // functions
  void initialise_valence_valence();
  void calculate_angular(int ki);

  // write another that returns pair <int, bool> = <index, val?> ?
  std::size_t find_valence_index(const DiracSpinor &Fa) const;
  std::size_t find_core_index(const DiracSpinor &Fa) const;
  std::size_t find_either_index(const DiracSpinor &Fa, bool &valenceQ) const;

  const std::vector<std::vector<double>> &get_y_abk(std::size_t a,
                                                    std::size_t b) const;
  const std::vector<std::vector<double>> &get_y_vck(std::size_t a,
                                                    std::size_t b) const;
  const std::vector<std::vector<double>> &get_y_vwk(std::size_t a,
                                                    std::size_t b) const;

private: // data
  const std::vector<DiracSpinor> *const c_orbs_ptr;
  const std::vector<DiracSpinor> *const v_orbs_ptr;
  const Grid *const rgrid_ptr;

  std::size_t num_initialised_vv = 0;
  std::size_t num_initialised_vc = 0;
  int m_largest_ki = -1; //-1 not valid, but need 0>x to hold true first time

  // Arrays to store Coulomb integrals. Note: highly non-rectangular
  // Make use of symmety (where appropriate)
  std::vector<std::vector<std::vector<std::vector<double>>>> m_y_abkr = {};
  std::vector<std::vector<std::vector<std::vector<double>>>> m_y_vckr = {};
  std::vector<std::vector<std::vector<std::vector<double>>>> m_y_vwkr = {};
  // Angular coeficients:
  std::vector<std::vector<std::vector<double>>> m_L_kakbk = {};
  std::vector<std::vector<std::vector<double>>> m_C_kakbk = {};
};
