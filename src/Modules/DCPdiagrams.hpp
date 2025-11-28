#pragma once
#include "Angular/CkTable.hpp"
#include "Angular/include.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/TDHF.hpp"

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Example module, designed as a "template" to help you add a new module. Note:
//! if you add a new Module, you must also update module_lists.hpp, so it will
//! be compiled into the rest of the code.
std::vector<Coulomb::meTable<double>>
compute_me(const DiracOperator::TensorOperator *const h1,
           const DiracOperator::TensorOperator *const h2,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited, const double omega,
           const ExternalField::CorePolarisation *dV1 = nullptr,
           const ExternalField::CorePolarisation *dV2 = nullptr);
double A1(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega);
double A2(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega);
double A3(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega);
double A4(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega);
double A5(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega);
double A6(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega);
double A7(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega);
double A8(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega);
double A9(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega);
double A10(const int K, const int k1, const int k2,
           std::vector<Coulomb::meTable<double>> ME_tables,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
           const DiracSpinor &v, double omega);
double A11(const int K, const int k1, const int k2,
           std::vector<Coulomb::meTable<double>> ME_tables,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
           const DiracSpinor &v, double omega);
double A12(const int K, const int k1, const int k2,
           std::vector<Coulomb::meTable<double>> ME_tables,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
           const DiracSpinor &v, double omega);
double A1(const int K, const int k1, const int k2,
          const DiracOperator::TensorOperator *const h1,
          const DiracOperator::TensorOperator *const h2,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, const double omega);
double AK_pol(const int K, const int k1, const int k2,
              std::vector<Coulomb::meTable<double>> ME_tables,
              const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
              const DiracSpinor &v, double omega);
void Polarisability(double diagrams, std::string static_operator, int K,
                    double omega);
double delta_A_KQ(const int K, const int k1, const int k2,
                  std::vector<Coulomb::meTable<double>> Table,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const DiracSpinor &Fw, const DiracSpinor &Fv, double omega);
void DCP_correction(double diagrams, std::string static_operator, int K,
                    double omega);
void DCPdiagrams(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
