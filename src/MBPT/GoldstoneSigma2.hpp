#pragma once
#include "CorrelationPotential.hpp"
#include <string>
#include <vector>
namespace HF {
class HartreeFock;
}

//! Many-body perturbation theory
namespace MBPT {

//******************************************************************************
/*!
@brief Correlation Potential operator (2nd order): Goldstone (diagram) technique
@details
Derives from CorrelationPotential.

Goldstone: Sigma stored in a matrix consisting of 4 'G' terms:

Four second-order diagrams:
  - Diagram (a):
\f[ G_a = \frac{|Q^k_amn><Q^k_amn|}{ [k][j] \, de_{amn} } \f]
  - Diagram (b) (exchange):
\f[ G_b = \frac{|Q^k_amn><P^k_amn|}{ [k][j] \, de_{amn} } \f]
  - Diagram (c):
\f[ G_c = \frac{|Q^k_nba><Q^k_nba|}{ [k][j] \, de_{nba} } \f]
  - Diagram (d) (exchange):
\f[ G_d = \frac{|Q^k_nba><P^k_nba|}{ [k][j] \, de_{nba} } \f]
where:
All indeces are summed over,
a & b are core states, n & m are virtual excited states,
k is multipolarity [Coloulmb expansion], and
\f$ de_{xyz} = e_v + e_x - e_y - e_z \f$

*/
class GoldstoneSigma2 final : public CorrelationPotential {
public:
  GoldstoneSigma2(const HF::HartreeFock *const in_hf,
                  const std::vector<DiracSpinor> &basis,
                  const Sigma_params &sigp, const rgrid_params &subgridp,
                  const std::vector<double> &en_list, const std::string &atom);

  GoldstoneSigma2 &operator=(const GoldstoneSigma2 &) = delete;
  GoldstoneSigma2(const GoldstoneSigma2 &) = delete;
  ~GoldstoneSigma2() = default;

protected:
  void fill_Sigma_k(GMatrix *Gmat, const int kappa, const double en) override final;
};

} // namespace MBPT
