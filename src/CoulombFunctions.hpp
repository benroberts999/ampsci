#pragma once
#include "DiracSpinor.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include <vector>

class Coulomb {
public:
  static void calculate_v_abk(const DiracSpinor &phi_a,
                              const DiracSpinor &phi_b, const int k,
                              std::vector<double> &vabk);

  // static double calculate_R_abcdk(); // Core-core? Core-valence?

  void initialise_core_v_abk_r();
  void form_core_v_abk_r(const std::vector<DiracSpinor> &c_orbitals);
  void form_coreval_v_amk_r(const std::vector<DiracSpinor> &c_orbitals,
                            const std::vector<DiracSpinor> &v_orbitals);

  // const std::vector<double> &get_v_aa0(std::size_t a) const;
  const std::vector<std::vector<double>> &get_v_abk(std::size_t a,
                                                    std::size_t b) const;

  const std::vector<double> &get_v_ab_k(std::size_t a, std::size_t b,
                                        int k) const;

private: // data
  // std::vector<std::vector<std::vector<std::vector<double>>>> m_arr_v_abk_r;
  // core-core: a-b symmetry [only calc/store for a>=b]
  std::vector<std::vector<std::vector<std::vector<double>>>> m_core_v_abk_r;

  // core-valence: [core-virtual] a must be core, m must be virtual
  // v[a][m][k]; v[CORE][NUM_VAL][MAX_K]
  std::vector<std::vector<std::vector<std::vector<double>>>> m_coreval_v_amk_r;

}; // namespace Coulomb
