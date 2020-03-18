#pragma once
#include <tuple>
#include <vector>
class UserInputBlock;
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
void fitParametric(const UserInputBlock &input, const Wavefunction &wf);

namespace FitParametric {
std::tuple<double, double>
performFit(const std::vector<AtomData::DiracSEnken> &states, int Z,
           const GridParameters &gp, const Nuclear::Parameters &nuc_params,
           bool green, bool fit_worst);
}

} // namespace Module
