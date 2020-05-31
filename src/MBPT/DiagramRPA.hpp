#pragma once
#include "Wavefunction/DiracSpinor.hpp" //?
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

//! RPA correction to matrix elements, using Diagram technique
class DiagramRPA {
public:
  //! Normal constructor: needs core to split basis: only uses basis
  DiagramRPA(const DiracOperator::TensorOperator *const h,
             const std::vector<DiracSpinor> &basis,
             const std::vector<DiracSpinor> &core);

  //! Second constructor: copies over W matrices (depend only on k/pi)
  DiagramRPA(const DiracOperator::TensorOperator *const h,
             const DiagramRPA *const drpa);

  //! Itterates the RPA equations for core electrons
  void rpa_core(const double omega, const bool print = true);

  //! Calculates RPA correction to matrix element: <A||dV||B>
  double dV(const DiracSpinor &Fa, const DiracSpinor &Fb,
            const bool first_order = false) const;

  //! Returns eps (convergance) of last rpa_core run
  double get_eps() const { return m_core_eps; }

  //! @brief Clears the t_am and t_ma RPA ME's [RPA ME's for hole-excited]
  //! @Details If a previous run failed, can clear t_am's + re-try
  void clear_tam();

private:
  const int m_k;
  const int m_pi; // need pi?
  const int m_imag;
  std::vector<DiracSpinor> holes{};
  std::vector<DiracSpinor> excited{};
  double m_omega = 0.0;
  double m_core_eps = 1.0;
  const double eps_targ = 1.0e-10;
  const int max_its = 100;

  // t0's never change
  std::vector<std::vector<double>> t0am{};
  std::vector<std::vector<double>> t0ma{};
  // t's updated each rpa_core itteration
  std::vector<std::vector<double>> tam{};
  std::vector<std::vector<double>> tma{};

  // Note: W's depend on rank (also parity)! Can re-use!?
  std::vector<std::vector<std::vector<std::vector<double>>>> Wanmb{};
  std::vector<std::vector<std::vector<std::vector<double>>>> Wabmn{};
  std::vector<std::vector<std::vector<std::vector<double>>>> Wmnab{};
  std::vector<std::vector<std::vector<std::vector<double>>>> Wmban{};

  void fill_W_matrix(const DiracOperator::TensorOperator *const h);
  void setup_ts(const DiracOperator::TensorOperator *const h);
};

} // namespace MBPT
