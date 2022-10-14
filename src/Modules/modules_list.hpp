#pragma once
// Add new module includes here:
// (Modules don't _need_ to be in Modules/ directory
#include "DMionisation/Module_AFBindingEnergy.hpp"
#include "DMionisation/Module_AFStepFunction.hpp"
#include "DMionisation/Module_atomicKernel.hpp"
#include "Modules/HFAnomaly.hpp"
#include "Modules/basic.hpp"
#include "Modules/exampleModule.hpp"
#include "Modules/isotopeShift.hpp"
#include "Modules/ladder.hpp"
#include "Modules/lifetimes.hpp"
#include "Modules/matrixElements.hpp"
#include "Modules/pnc.hpp"
#include "Modules/polarisability.hpp"
#include "Modules/qed.hpp"
#include "Modules/runModules.hpp"
#include "Modules/screeningFactors.hpp"
#include "Modules/testFeynman.hpp"

#include <iostream>
#include <map>
#include <string>
#include <utility>
class Wavefunction;
namespace IO {
class InputBlock;
} // namespace IO

namespace Module {

//! function (function ptr) signature for modules
using module_function_t = void (*)(const IO::InputBlock &input,
                                   const Wavefunction &wf);

//! List of all available modules as pair (name, function):
//! {"ModuleName", &ModuleName}. You must add any new modules to this list
static const std::map<std::string, module_function_t> module_list{
    {"Tests", &tests},
    {"WriteOrbitals", &writeOrbitals},
    {"continuum", &continuum},
    {"AtomicKernel", &atomicKernel},
    {"BohrWeisskopf", &BohrWeisskopf},
    {"HFAnomaly", &HFAnomaly},
    {"HF_rmag", &HF_rmag},
    {"BW_eta_sp", &BW_eta_sp},
    {"screeningFactors", &screeningFactors},
    {"pnc", &calculatePNC},
    {"QED", &QED},
    {"testFeynman", &testFeynman},
    {"matrixElements", &matrixElements},
    {"lifetimes", &lifetimes},
    {"polarisability", &polarisability},
    {"dynamicPolarisability", &dynamicPolarisability},
    {"transitionPolarisability", &transitionPolarisability},
    {"structureRad", &structureRad},
    {"fieldShift", &fieldShift},
    {"ladder", &ladder},
    {"exampleModule", &exampleModule}};

} // namespace Module
