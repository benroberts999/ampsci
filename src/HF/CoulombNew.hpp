#pragma once
#include <vector>
class DiracSpinor;
class Grid;

constexpr static bool check_bounds = true;

class CoulombNew {

  CoulombNew(const std::vector<DiracSpinor> *const in_a_orbs,
             const std::vector<DiracSpinor> *const in_b_orbs = nullptr);

private:
  const DiracSpinor *const m_a_orbs;
  const DiracSpinor *const m_b_orbs;
  const bool m_aisb;
  const std::vector<std::vector<std::vector<std::vector<double>>>> m_y_abkr;

private:
}
