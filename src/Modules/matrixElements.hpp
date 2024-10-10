#pragma once
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Calculates matrix elements of any tensor operator, with RPA
void matrixElements(const IO::InputBlock &input, const Wavefunction &wf);

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
void structureRad(const IO::InputBlock &input, const Wavefunction &wf);

//! Calculates matrix elements for CI wavefunctions
void CI_matrixElements(const IO::InputBlock &input, const Wavefunction &wf);

void normalisation(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
