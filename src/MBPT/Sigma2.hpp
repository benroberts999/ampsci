#pragma once
#include "Angular/Angular.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Wavefunction/DiracSpinor.hpp"

namespace MBPT {

//! Reduced two-body Sigma (2nd order correlation) operator. Sum of 6 diagrams
double S_Sigma2(int k, const DiracSpinor &v, const DiracSpinor &w,
                const DiracSpinor &x, const DiracSpinor &y,
                const Coulomb::QkTable &qk,
                const std::vector<DiracSpinor> &core,
                const std::vector<DiracSpinor> &excited,
                const Angular::SixJTable *SJ = nullptr);

//! (sum of three Goldstone diagrams)
double S_Sigma2_a(int k, const DiracSpinor &v, const DiracSpinor &w,
                  const DiracSpinor &x, const DiracSpinor &y,
                  const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const Angular::SixJTable *SJ = nullptr);

double S_Sigma2_b(int k, const DiracSpinor &v, const DiracSpinor &w,
                  const DiracSpinor &x, const DiracSpinor &y,
                  const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const Angular::SixJTable *SJ = nullptr);

double S_Sigma2_c(int k, const DiracSpinor &v, const DiracSpinor &w,
                  const DiracSpinor &x, const DiracSpinor &y,
                  const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const Angular::SixJTable *SJ = nullptr);

double S_Sigma2_d(int k, const DiracSpinor &v, const DiracSpinor &w,
                  const DiracSpinor &x, const DiracSpinor &y,
                  const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const Angular::SixJTable *SJ = nullptr);

} // namespace MBPT