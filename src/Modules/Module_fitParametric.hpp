#pragma once
#include "Physics/NuclearPotentials.hpp" //Nuclear::Parameters
#include <tuple>
#include <vector>
class UserInputBlock;
class Wavefunction;
// struct DiracSEnken;
#include "Physics/AtomData.hpp"
struct GridParameters;

namespace Module {
std::tuple<double, double>
fitParametric_performFit(const std::vector<AtomData::DiracSEnken> &states,
                         int Z, const GridParameters &gp,
                         const Nuclear::Parameters &nuc_params, bool green,
                         bool fit_worst);

void fitParametric(const UserInputBlock &input, const Wavefunction &wf);
} // namespace Module
