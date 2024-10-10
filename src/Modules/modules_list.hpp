#pragma once
// Add new module includes here:
// (Modules don't _need_ to be in Modules/ directory
#include "Kionisation/Module_Kionisation.hpp"
#include "Modules/Breit.hpp"
#include "Modules/HFAnomaly.hpp"
#include "Modules/VQE.hpp"
#include "Modules/basic.hpp"
#include "Modules/exampleModule.hpp"
#include "Modules/isotopeShift.hpp"
#include "Modules/ladder.hpp"
#include "Modules/lifetimes.hpp"
#include "Modules/matrixElements.hpp"
#include "Modules/muonic.hpp"
#include "Modules/pnc.hpp"
#include "Modules/polarisability.hpp"
#include "Modules/qed.hpp"
#include "Modules/runModules.hpp"
#include "Modules/screeningFactors.hpp"
#include "Modules/thirdOrder.hpp"

#include <iostream>
#include <string>
#include <tuple>
#include <vector>
class Wavefunction;
namespace IO {
class InputBlock;
} // namespace IO

namespace Module {

//! function (function ptr) signature for modules
using module_function_t = void (*)(const IO::InputBlock &input,
                                   const Wavefunction &wf);

struct ModuleInfo {
  std::string name;
  module_function_t function;
  std::string description;
};

//! List of all available modules as pair (name, function):
//! {"ModuleName", &ModuleName, "ModuleDescription"}.
//! You must add any new modules to this list.
static const std::vector<ModuleInfo> module_list{
    {"Tests", &tests, "Some basic wavefunction numerical tests"},
    {"testBasis", &testBasis, "Tests of basis and spectrum"},
    {"WriteOrbitals", &writeOrbitals, "Write orbitals to disk for plotting"},
    {"matrixElements", &matrixElements,
     "Calculates matrix elements of any operator"},
    {"CI_matrixElements", &CI_matrixElements,
     "Calculates matrix elements of any operator for CI wavefunctions"},
    {"thirdOrderME", &thirdOrderME, "Calculates Third-order matrix elements"},
    {"normalisation", &normalisation, "normalisation"},
    {"lifetimes", &lifetimes, "Calculate radiative lifetimes (E1, E2, M1)"},
    {"polarisability", &polarisability, "Calculates static polarisabilities"},
    {"dynamicPolarisability", &dynamicPolarisability,
     "Calculates dynamic polarisabilities"},
    {"transitionPolarisability", &transitionPolarisability,
     "Calculates transition polarisabilities"},
    {"structureRad", &structureRad,
     "Calculates Struct. Rad + Normalisation corrections to MEs"},
    {"fieldShift", &fieldShift,
     "Calculates field-shift constants (isotope shift)"},
    {"fieldShift2", &fieldShift2,
     "Calculates field-shift constants (isotope shift)"},
    {"QED", &QED, "QED corrections to energies/matrix elements"},
    {"Breit", &Breit, "Breit corrections to energies"},
    {"ladder", &ladder, "Calculates ladder diagrams and energy corrections"},
    {"Kionisation", &Kionisation, "Calculate atomic ionisation form-factors"},
    {"continuum", &continuum, "Compute and use continuum wavefunctions"},
    {"HFAnomaly", &HFAnomaly,
     "Calculates Bohr-Weisskopf effect and hyperfine anomaly"},
    {"screeningFactors", &screeningFactors,
     "Calculates Feynman electron screening factors"},
    {"pnc", &calculatePNC, "Calculates APV amplitudes"},
    {"muonPV", &muonPV, "For testing/playing with muonic PV"},
    {"muon", &muon,
     "Calculating muonic wavefunctions, energies, matrix elements"},
    {"VQE", &VQE, "For testing/playing with VQE method"},

    {"exampleModule", &exampleModule, "A short description of the module"}};

} // namespace Module
