#pragma once
#include "DiracSpinor.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include "Wigner_369j.hpp"
#include <vector>

class Coulomb {
public:
  static void calculate_v_abk(const DiracSpinor &phi_a,
                              const DiracSpinor &phi_b, const int k,
                              std::vector<double> &vabk);

  struct State {
    const int n;
    const int k;
    State(int in_n, int in_k) : n(in_n), k(in_k) {}
    bool operator==(const State &other) const {
      return n == other.n && k == other.k;
    }
    bool operator!=(const State &other) const { return !(*this == other); }
  };

  // hold a list of all states currently calculated, give them unique index
  // Index is equal to index of this list!
  std::vector<State> nka_list;
  std::vector<State> nkb_list;
  // if valence-core, a MUST be the valence state! OK????

  std::vector<int> ki_list;
  std::vector<int> twoj_list;
  int m_largest_ki = -1; //-1 not valid, but need 0>x to hold true first time

  std::vector<std::vector<std::vector<std::vector<double>>>> m_v_abkr;
  std::vector<std::vector<std::vector<double>>> m_angular_L_kakbk;
  std::vector<std::vector<std::vector<double>>> m_angular_C_kakbk;

  // nb: this not const, cos used in construction!
  std::vector<std::vector<double>> &get_vab_kr(const DiracSpinor &phi_a,
                                               const DiracSpinor &phi_b);

  const std::vector<double> &get_vabk_r(const DiracSpinor &phi_a,
                                        const DiracSpinor &phi_b, int k);
  void form_v_abk(const std::vector<DiracSpinor> &a_orbitals,
                  const std::vector<DiracSpinor> &b_orbitals);
  void form_v_abk(const DiracSpinor &phi_a,
                  const std::vector<DiracSpinor> &b_orbitals);

  void initialise_v_abkr(const std::vector<DiracSpinor> &a_orbitals,
                         const std::vector<DiracSpinor> &b_orbitals);
  void initialise_v_abkr(const std::vector<DiracSpinor> &orbitals);

  void extend_v_abkr(const DiracSpinor &phi_a,
                     const std::vector<DiracSpinor> &b_orbitals);
  void calculate_angular(int ki);
  double get_angular_C_kiakibk(int kia, int kib, int k);
  const std::vector<double> &get_angular_C_kiakib_k(int kia, int kib);
  double get_angular_L_kiakibk(int kia, int kib, int k);
  const std::vector<double> &get_angular_L_kiakib_k(int kia, int kib);
};
