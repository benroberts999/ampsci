#pragma once
#include "Angular/Angular_tables.hpp"
#include "Coulomb/YkTable.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <vector>
class Grid;
// namespace Angular {
// class Ck_ab;
// class SixJ;
// } // namespace Angular
// namespace Coulomb {
// class YkTable;
// } // namespace Coulomb

namespace MBPT {

class CorrelationPotential {
public:
  CorrelationPotential(const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &excited);
  // XX Needs to know which are the core states!

  // No, give the class () operator!
  // One for each energy? One each kappa?
  // Probably: each kappa stored in here, each at single energy!

  DiracSpinor operator()(const DiracSpinor &Fv) const { return Sigma2Fv(Fv); }
  double operator()(const DiracSpinor &Fv, const DiracSpinor &Fw) const {
    return Sigma2vw(Fv, Fw);
  }

  // XXX Make these all const!
  DiracSpinor Sigma2Fv(const DiracSpinor &Fv) const;

  double Sigma2vw(const DiracSpinor &Fv, const DiracSpinor &Fw) const;
  // double Sigma2vw(const DiracSpinor &Fv) const; // kill this, will be
  // confusing!

private:
  std::vector<DiracSpinor> m_core;
  std::vector<DiracSpinor> m_excited;
  Coulomb::YkTable m_yec; // constains Ck
  int m_maxk;
  Angular::SixJ m_6j;
};

} // namespace MBPT
