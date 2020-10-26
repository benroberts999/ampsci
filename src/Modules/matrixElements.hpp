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

//! Calculates Structure Radiation + Normalisation of States
/*!
Note: Most input options are similar to MatrixElements module:

  Module::structureRad{ operator; options; rpa; printBoth; onlyDiagonal; omega;
n_minmax;  }

n_minmax: is input as list of ints:
 * n_minmax = min,max;
 * min: minimum n for core states kept in summations
 * max: maximum n for excited states kept in summations

For explanation of the rest, see MatrixElements module.
*/
void structureRad(const IO::UserInputBlock &input, const Wavefunction &wf);

//! Calculates state lifetimes (using E1 and E2 only). nb: HF energies
void calculateLifetimes(const IO::UserInputBlock &input,
                        const Wavefunction &wf);

//! Returns a ptr to the requested operator, with given properties
std::unique_ptr<DiracOperator::TensorOperator>
generateOperator(const IO::UserInputBlock &input, const Wavefunction &wf,
                 bool print = true);

std::unique_ptr<DiracOperator::TensorOperator>
generateOperator(std::string_view oper_name, const IO::UserInputBlock &input,
                 const Wavefunction &wf, bool print);

// ------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_E1(const IO::UserInputBlock &input, const Wavefunction &wf,
            bool print = true);
std::unique_ptr<DiracOperator::TensorOperator>
generate_Ek(const IO::UserInputBlock &input, const Wavefunction &wf,
            bool print = true);
std::unique_ptr<DiracOperator::TensorOperator>
generate_M1(const IO::UserInputBlock &input, const Wavefunction &wf,
            bool print = true);
std::unique_ptr<DiracOperator::TensorOperator>
generate_hfs(const IO::UserInputBlock &input, const Wavefunction &wf,
             bool print = true);
std::unique_ptr<DiracOperator::TensorOperator>
generate_r(const IO::UserInputBlock &input, const Wavefunction &wf,
           bool print = true);
std::unique_ptr<DiracOperator::TensorOperator>
generate_pnc(const IO::UserInputBlock &input, const Wavefunction &wf,
             bool print = true);
std::unique_ptr<DiracOperator::TensorOperator>
generate_Hrad_el(const IO::UserInputBlock &input, const Wavefunction &wf,
                 bool print = true);
std::unique_ptr<DiracOperator::TensorOperator>
generate_Hrad_mag(const IO::UserInputBlock &input, const Wavefunction &wf,
                  bool print = true);

// ------------------
const std::vector<
    std::pair<std::string, std::unique_ptr<DiracOperator::TensorOperator> (*)(
                               const IO::UserInputBlock &input,
                               const Wavefunction &wf, bool print)>>
    operator_list{{"E1", &generate_E1},
                  {"Ek", &generate_Ek},
                  {"M1", &generate_M1},
                  {"hfs", &generate_hfs},
                  {"r", &generate_r},
                  {"pnc", &generate_pnc},
                  {"Hrad_el", &generate_Hrad_el},
                  {"Hrad_mag", &generate_Hrad_mag}};

} // namespace Module
