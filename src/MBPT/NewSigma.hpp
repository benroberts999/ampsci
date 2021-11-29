#pragma once
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "MBPT/RDMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <cassert>
#include <vector>
//
#include <iostream>

namespace MBPT {

struct SigData {
  int kappa;
  double en;
  RDMatrix<double> Sig;
};

// This is a test - I'd like to re-write the Correlation potential class
// but still not sure exactly best way to do it..

//******************************************************************************

class Sigma {
public:
  std::vector<SigData> m_Sigmas{};
  std::shared_ptr<const Grid> m_rgrid; //?
  double m_r0, m_rmax;
  std::size_t m_stride;
  std::size_t m_i0, m_size;

  std::vector<double> m_lambda{}; //?

public:
  Sigma(double r0, double rmax, std::size_t stride,
        std::shared_ptr<const Grid> rgrid)
      : m_rgrid(rgrid),
        m_r0(r0),
        m_rmax(rmax),
        m_stride(stride),
        m_i0(rgrid->getIndex(r0)),
        m_size((rgrid->getIndex(rmax) - m_i0 + 1) / m_stride) {}

  //! Find Sigma with kappa, and energy en (within +/- 1%).
  const RDMatrix<double> *get(int kappa, double en = 0.0) const {
    if (en != 0.0) {
      // Find the sigma with the correct energy (to within +/- 1%)
      const auto lambda = [=](const auto &tS) {
        return (kappa == tS.kappa &&
                std::abs((en - tS.en) / (en + tS.en)) < 1.0e-2);
      };
      const auto it = std::find_if(m_Sigmas.cbegin(), m_Sigmas.cend(), lambda);
      if (it == m_Sigmas.cend())
        return nullptr;
      return &(it->Sig);
    } else {
      // just find the First Sigma that has correct energy
      const auto lambda = [=](const auto &tS) { return kappa == tS.kappa; };
      const auto it = std::find_if(m_Sigmas.cbegin(), m_Sigmas.cend(), lambda);
      if (it == m_Sigmas.cend())
        return nullptr;
      return &(it->Sig);
    }
  }

  //
};

} // namespace MBPT
