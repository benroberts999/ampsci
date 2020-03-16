#pragma once
#include <utility>
#include <vector>
// class DiracSpinor;
// class Grid;
#include "Angular/Angular_tables.hpp"
//
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <cmath>

constexpr static bool check_bounds = true;

class CoulombNew {
public:
  CoulombNew(const Grid *const in_grid,
             const std::vector<DiracSpinor> *const in_a_orbs,
             const std::vector<DiracSpinor> *const in_b_orbs = nullptr)
      : m_a_orbs(in_a_orbs),                                    //
        m_b_orbs(in_b_orbs == nullptr ? in_a_orbs : in_b_orbs), //
        m_grid(in_grid),                                        //
        m_aisb([&]() {
          return (in_b_orbs == nullptr || in_a_orbs == in_b_orbs) ? true
                                                                  : false;
        }()) //
  {
    update_y_ints();
  }
  CoulombNew &operator=(const CoulombNew &) = default;
  CoulombNew(const CoulombNew &) = default;
  ~CoulombNew() = default;

private:
  const std::vector<DiracSpinor> *const m_a_orbs;
  const std::vector<DiracSpinor> *const m_b_orbs;
  const Grid *const m_grid;
  const bool m_aisb;

private:
  std::size_t a_size = 0;
  std::size_t b_size = 0;
  std::vector<std::vector<std::vector<std::vector<double>>>> m_y_abkr = {};

public:
  Angular::Ck_ab m_Ck = Angular::Ck_ab(); //???

public:
  void update_y_ints();
  void update_y_ints(const DiracSpinor &Fn);

  const std::vector<double> &get_yk_ab(const int k, const DiracSpinor &Fa,
                                       const DiracSpinor &Fb) const;
  const std::vector<std::vector<double>> &get_y_ab(const DiracSpinor &Fa,
                                                   const DiracSpinor &Fb) const;

private:
  int max_tj() const;
  void resize_y();

public:
  //! Calculates Hartree Screening functions \f$y^k_{ab}(r)\f$
  static void calculate_y_ijk(const DiracSpinor &Fa, const DiracSpinor &Fb,
                              const int k, std::vector<double> &vabk,
                              const std::size_t maxi = 0);

  static std::pair<int, int> k_minmax(const DiracSpinor &a,
                                      const DiracSpinor &b);
};
