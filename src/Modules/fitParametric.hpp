#pragma once
#include <tuple>
#include <vector>
namespace IO {
class InputBlock;
}
class Wavefunction;
struct GridParameters;
namespace AtomData {
struct DiracSEnken;
}
namespace Nuclear {
class Nucleus;
}

namespace Module {

//! Performs fit to Greens/Teitz parametric potential. See input example
void fitParametric(const IO::InputBlock &input, const Wavefunction &wf);

// namespace FitParametric {
// std::tuple<double, double>
// performFit(const std::vector<AtomData::DiracSEnken> &states, int Z,
//            const GridParameters &gp, const Nuclear::Nucleus &nuc_params,
//            bool green, bool fit_worst);
// }

} // namespace Module
