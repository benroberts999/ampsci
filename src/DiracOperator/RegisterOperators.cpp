#include "DiracOperator/GenerateOperator.hpp"
#include "DiracOperator/Operators/EM_multipole.hpp"
#include "DiracOperator/Operators/Ek.hpp"
#include "DiracOperator/Operators/FieldShift.hpp"
#include "DiracOperator/Operators/M1.hpp"
#include "DiracOperator/Operators/PNC.hpp"
#include "DiracOperator/Operators/QED.hpp"
#include "DiracOperator/Operators/RadialF.hpp"
#include "DiracOperator/Operators/hfs.hpp"
#include "DiracOperator/Operators/jls.hpp"
#include "DiracOperator/Operators/p.hpp"

namespace {

const DiracOperator::Register<DiracOperator::E1> r_E1{
  "E1", "Electric dipole (moment), length form: -|e|r"};

const DiracOperator::Register<DiracOperator::E1v> r_E1v{
  "E1v", "Electric dipole, v-form"};

const DiracOperator::Register<DiracOperator::E2> r_E2{
  "E2", "Electric quadrupole moment operator"};

const DiracOperator::Register<DiracOperator::Ek> r_Ek{
  "Ek", "Electric multipole moment operator, in low qr limit"};

const DiracOperator::Register<DiracOperator::ialpha> r_ialpha{
  "ialpha", "i*alpha (propto E1v)"};

const DiracOperator::Register<DiracOperator::M1> r_M1{
  "M1", "Magnetic dipole (relativistic formula)"};

const DiracOperator::Register<DiracOperator::M1nr> r_M1nr{
  "M1nr", "Non-relativistic M1"};

const DiracOperator::Register<DiracOperator::Multipole> r_Multipole{
  "Multipole",
  "Multipole transition operators (Vector,Axial,Scalar,Pseudoscalar)"};

const DiracOperator::Register<DiracOperator::hfs> r_hfs{
  "hfs", "Hyperfine structure k-pole operators"};

const DiracOperator::Register<DiracOperator::fieldshift> r_fieldshift{
  "fieldshift", "Field-shift F(r) operator"};

const DiracOperator::Register<DiracOperator::RadialF> r_r{
  "r", "radial (scalar) |r|"};

const DiracOperator::Register<DiracOperator::sigma_r> r_sigma_r{
  "sigma_r", "scalar sigma.r operator"};

const DiracOperator::Register<DiracOperator::PNCnsi> r_pnc{"pnc",
                                                           "NSI PNC operator"};

const DiracOperator::Register<DiracOperator::Vrad> r_Vrad{
  "Vrad", "QED Radiative potential"};

const DiracOperator::Register<DiracOperator::MLVP> r_MLVP{
  "MLVP", "Magnetic-Loop vacuum polarisation vertex correction to HFS"};

const DiracOperator::Register<DiracOperator::p> r_p{"p", "Momentum operator"};

const DiracOperator::Register<DiracOperator::l> r_l{"l", "Orbital L"};

const DiracOperator::Register<DiracOperator::s> r_s{"s", "Spin S (not sigma)"};

} // namespace
