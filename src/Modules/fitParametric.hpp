#pragma once
#include <tuple>
#include <vector>
namespace IO {
class UserInputBlock;
}
class Wavefunction;
struct GridParameters;
namespace AtomData {
struct DiracSEnken;
}
namespace Nuclear {
struct Parameters;
}

namespace Module {

//! Performs fit to Greens/Teitz parametric potential. See input example
void fitParametric(const IO::UserInputBlock &input, const Wavefunction &wf);

namespace FitParametric {
std::tuple<double, double>
performFit(const std::vector<AtomData::DiracSEnken> &states, int Z,
           const GridParameters &gp, const Nuclear::Parameters &nuc_params,
           bool green, bool fit_worst);
}

} // namespace Module
