#pragma once
// Add new module includes here:
// (Modules don't _need_ to be in Modules/ directory
#include "DMionisation/Module_atomicKernal.hpp"
#include "Modules/HFAnomaly.hpp"
#include "Modules/continuum.hpp"
#include "Modules/exampleModule.hpp"
#include "Modules/isotopeShift.hpp"
#include "Modules/ladder.hpp"
#include "Modules/matrixElements.hpp"
#include "Modules/pnc.hpp"
#include "Modules/polarisability.hpp"
#include "Modules/qed.hpp"
#include "Modules/runModules.hpp"
#include "Modules/screeningFactors.hpp"
#include "Modules/testFeynman.hpp"
#include "Modules/tests.hpp"
#include "Modules/writeOrbitals.hpp"
#include <iostream>
#include <string>
#include <utility>
#include <vector>
class Wavefunction;
namespace IO {
class InputBlock;
} // namespace IO

namespace Module {

// Add new modules to this list:
static const std::vector<std::pair<
    std::string, void (*)(const IO::InputBlock &input, const Wavefunction &wf)>>
    module_list{{"Tests", &Module_tests},
                {"WriteOrbitals", &writeOrbitals},
                {"AtomicKernal", &atomicKernal},
                {"BohrWeisskopf", &calculateBohrWeisskopf},
                {"HFAnomaly", &HFAnomaly},
                {"HF_rmag", &HF_rmag},
                {"screeningFactors", &screeningFactors},
                {"BW_eta_sp", &BW_eta_sp},
                {"pnc", &calculatePNC},
                {"vertexQED", &vertexQED},
                {"QED", &QED},
                {"testFeynman", &testFeynman},
                {"matrixElements", &matrixElements},
                {"lifetimes", &calculateLifetimes},
                {"polarisability", &polarisability},
                {"dynamicPolarisability", &dynamicPolarisability},
                {"transitionPolarisability", &transitionPolarisability},
                {"structureRad", &structureRad},
                {"fieldShift", &fieldShift},
                {"continuum", &continuum},
                {"ladder", &ladder},
                {"exampleModule", &exampleModule}};

inline void list_modules() {
  for (auto &[name, func] : module_list) {
    std::cout << "  " << name << '\n';
  }
}

} // namespace Module
