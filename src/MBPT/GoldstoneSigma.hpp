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
class GoldstoneSigma final : public CorrelationPotential {
public:
  GoldstoneSigma(const HF::HartreeFock *const in_hf,
                 const std::vector<DiracSpinor> &basis,
                 const Sigma_params &sigp, const rgrid_params &subgridp,
                 // const std::vector<DiracSpinor> &valence,
                 const std::string &atom);

  GoldstoneSigma &operator=(const GoldstoneSigma &) = delete;
  GoldstoneSigma(const GoldstoneSigma &) = delete;
  ~GoldstoneSigma() = default;

  void formSigma(int kappa, double en, int n = 0) override final;

protected:
  // make static!? Or move to base class ? BASE!
  void Sigma2(GMatrix *Gmat_D, GMatrix *Gmat_X, int kappa, double en);
};

} // namespace MBPT
