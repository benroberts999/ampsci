#pragma once
#include "DiracOperator/include.hpp" //For E1 operator
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"
#include <optional>
#include <string>
#include <vector>

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

std::vector<Coulomb::meTable<double>>
compute_me_3f(const DiracOperator::TensorOperator *hpnc,
              const DiracOperator::TensorOperator *he1,
              const DiracOperator::TensorOperator *const hfs,
              const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
              const DiracSpinor &v,
              const ExternalField::CorePolarisation *dV_pnc = nullptr,
              const ExternalField::CorePolarisation *dV_e1 = nullptr,
              const ExternalField::DiagramRPA *dV_hf = nullptr);
double h1(std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
          const DiracSpinor &v, int I2, int Fv2, int Fw2, int two_k);
double h2(std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
          const DiracSpinor &v, int I2, int Fv2, int Fw2, int two_k);

void hf_pert_weak(const IO::InputBlock &input, const Wavefunction &wf);
} // namespace Module
