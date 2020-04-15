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

namespace Module {

//! Calculates matrix elements of any tensor operator, with RPA
void matrixElements(const IO::UserInputBlock &input, const Wavefunction &wf);

//! Calculates Bohr-Weisskopf effect for hyperfine structure
void calculateBohrWeisskopf(const IO::UserInputBlock &input,
                            const Wavefunction &wf);

//! Calculates state lifetimes (using E1 and E2 only). nb: HF energies
void calculateLifetimes(const IO::UserInputBlock &input,
                        const Wavefunction &wf);

//! Returns a ptr to the requested operator, with given properties
std::unique_ptr<DiracOperator::TensorOperator>
generateOperator(const IO::UserInputBlock &input, const Wavefunction &wf);

// ------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_E1(const IO::UserInputBlock &input, const Wavefunction &wf);
std::unique_ptr<DiracOperator::TensorOperator>
generate_Ek(const IO::UserInputBlock &input, const Wavefunction &wf);
std::unique_ptr<DiracOperator::TensorOperator>
generate_M1(const IO::UserInputBlock &input, const Wavefunction &wf);
std::unique_ptr<DiracOperator::TensorOperator>
generate_hfs(const IO::UserInputBlock &input, const Wavefunction &wf);
std::unique_ptr<DiracOperator::TensorOperator>
generate_r(const IO::UserInputBlock &input, const Wavefunction &wf);
std::unique_ptr<DiracOperator::TensorOperator>
generate_pnc(const IO::UserInputBlock &input, const Wavefunction &wf);
std::unique_ptr<DiracOperator::TensorOperator>
generate_Hrad_el(const IO::UserInputBlock &input, const Wavefunction &wf);
std::unique_ptr<DiracOperator::TensorOperator>
generate_Hrad_mag(const IO::UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
