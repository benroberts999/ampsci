#pragma once
#include <string>
#include <vector>
class Wavefunction;
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}
namespace HF {
class HartreeFock;
}

namespace MBPT {

class DiagramRPA {
public:
  DiagramRPA(const DiracOperator::TensorOperator *const h,
             const std::vector<DiracSpinor> &basis,
             const std::vector<DiracSpinor> &core);
  // Fill occ/exc, W's, and t0's; size t's

  rpa_core(const double omega, int max_its = 100, const bool print = true);

  double dV(const DiracSpinor &Fa, const DiracSpinor &Fb) const;

  //! Returns eps (convergance) of last rpa_core run
  double get_eps() const { return m_core_eps; }

  //! @brief Clears the t_am and t_ma RPA ME's
  void clear_tam() {
    for (auto i = 0ul; i < tam.size(); i++) {
      for (auto j = 0ul; j < tam.size(); j++) {
        tam[i][j] = 0.0;
        tma[i][j] = 0.0;
      }
    }
  }

private:
  const int m_k;
  const int m_imag; // need pi?
  std::vector<DiracSpinor> holes{};
  std::vector<DiracSpinor> excited{};
  double m_omega = 0.0;
  double m_core_eps = 1.0;

  // Note: W's depend ONLY on rank! can re-use!
  std::vector<std::vector<std::vector<std::vector<double>>>> Wanmb{};
  std::vector<std::vector<std::vector<std::vector<double>>>> Wabmn{};
  std::vector<std::vector<std::vector<std::vector<double>>>> Wmnab{};
  std::vector<std::vector<std::vector<std::vector<double>>>> Wmban{};

  std::vector<std::vector<double>> tam{};
  std::vector<std::vector<double>> tma{};
  std::vector<std::vector<double>> t0am{};
  std::vector<std::vector<double>> t0ma{};
};

} // namespace MBPT
