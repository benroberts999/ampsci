#pragma once
#include "CorrelationPotential.hpp"
class Grid;
namespace HF {
class HartreeFock;
}

//! Many-body perturbation theory
namespace MBPT {

//******************************************************************************
/*!
@brief
@details
*/
class GoldstoneSigma2 : public CorrelationPotential {
public:
  GoldstoneSigma2(const HF::HartreeFock *const in_hf,
                  const std::vector<DiracSpinor> &basis,
                  const Sigma_params &sigp, const rgrid_params &subgridp,
                  const std::vector<double> &en_list, const std::string &atom);

  GoldstoneSigma2 &operator=(const GoldstoneSigma2 &) = delete;
  GoldstoneSigma2(const GoldstoneSigma2 &) = delete;
  ~GoldstoneSigma2() = default;

protected:
  void fill_Sigma_k(GMatrix *Gmat, const int kappa, const double en) override;
};

} // namespace MBPT
