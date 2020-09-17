#pragma once
#include <memory>
#include <string>
#include <vector>
class Wavefunction;
namespace IO {
class UserInputBlock;
}
namespace DiracOperator {
class TensorOperator;
}
namespace MBPT {
class FeynmanSigma;
}

namespace Module {

void testFeynman(const IO::UserInputBlock &input, const Wavefunction &wf);

namespace Feyn {
void test_Q(const Wavefunction &wf, const MBPT::FeynmanSigma &Sigma);

void test_green(const Wavefunction &wf, const MBPT::FeynmanSigma &Sigma,
                double omre, const std::vector<double> &omim_v);

void test_GQ(const Wavefunction &wf, const MBPT::FeynmanSigma &Sigma,
             double omre, const std::vector<double> &omim_v);

void test_pol(const Wavefunction &wf, const MBPT::FeynmanSigma &Sigma,
              double omre, const std::vector<double> &omim_v,
              int max_l_excited);
} // namespace Feyn

} // namespace Module
