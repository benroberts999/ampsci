#pragma once
#include "Physics/Nuclear.hpp"
#include "Dirac/Wavefunction.hpp"
#include <tuple>
#include <vector>
class UserInputBlock;
struct DiracSEnken;

namespace Module {
std::tuple<double, double> fitParametric_performFit(
    const std::vector<DiracSEnken> &states, int Z, const GridParameters &gp,
    const Nuclear::Parameters &nuc_params, bool green, bool fit_worst);

void fitParametric(const UserInputBlock &input, const Wavefunction &wf);
} // namespace Module
