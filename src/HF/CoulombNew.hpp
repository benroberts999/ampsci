#pragma once
#include "Angular/Angular_tables.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <utility>
#include <vector>

constexpr static bool check_bounds = true;

class YkTable {
public:
  YkTable(const Grid *const in_grid,
          const std::vector<DiracSpinor> *const in_a_orbs,
          const std::vector<DiracSpinor> *const in_b_orbs = nullptr);
  YkTable &operator=(const YkTable &) = default;
  YkTable(const YkTable &) = default;
  ~YkTable() = default;

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
  static std::pair<int, int> k_minmax(const DiracSpinor &a,
                                      const DiracSpinor &b);
};
