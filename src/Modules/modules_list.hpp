#pragma once
// Add new module includes here:
// (Modules don't _need_ to be in Modules/ directory
#include "DMionisation/Module_atomicKernal.hpp"
#include "Modules/HFAnomaly.hpp"
#include "Modules/continuum.hpp"
#include "Modules/exampleModule.hpp"
#include "Modules/fitParametric.hpp"
#include "Modules/isotopeShift.hpp"
#include "Modules/laddertests.hpp"
#include "Modules/matrixElements.hpp"
#include "Modules/pnc.hpp"
#include "Modules/polarisability.hpp"
#include "Modules/qed.hpp"
#include "Modules/runModules.hpp"
#include "Modules/screeningFactors.hpp"
#include "Modules/testFeynman.hpp"
#include "Modules/tests.hpp"
//
#include <string>
#include <utility>
#include <vector>
class Wavefunction;
namespace IO {
class InputBlock;
class InputBlock;
} // namespace IO

namespace Module {

// Add new modules to this list:
static const std::vector<std::pair<
    std::string, void (*)(const IO::InputBlock &input, const Wavefunction &wf)>>
    module_list{{"Tests", &Module_tests},
                {"WriteOrbitals", &writeOrbitals},
                {"AtomicKernal", &atomicKernal},
                {"FitParametric", &fitParametric},
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
                {"structureRad", &structureRad},
                {"fieldShift", &fieldShift},
                {"continuum", &continuum},
                {"ladder", &ladder},
                {"exampleModule", &exampleModule}};

} // namespace Module
