#pragma once
#include "Maths/LinAlg_MatrixVector.hpp"
#include <string>
#include <utility>
class DiracSpinor;
class Wavefunction;
class Grid;

namespace SplineBasis {

// Forms the basis orbitals (expanded in terms of splines)
std::vector<DiracSpinor> form_basis(const std::string &states_str,
                                    const std::size_t n_spl,
                                    const std::size_t k_spl,
                                    const double r0_spl, const double rmax_spl,
                                    const Wavefunction &wf);

// Forms the underlying spline basis (which is not kept)
std::vector<DiracSpinor>
form_spline_basis(const int kappa, const std::size_t n_states,
                  const std::size_t k_spl, const double r0_spl,
                  const double rmax_spl, const Grid &rgrid, const double alpha);

std::pair<LinAlg::SqMatrix, LinAlg::SqMatrix>
fill_Hamiltonian_matrix(const std::vector<DiracSpinor> &spl_basis,
                        const Wavefunction &wf);

// Expands basis orbitals in terms of spline orbitals, by diagonalising
// Hamiltonian:
void expand_basis_orbitals(std::vector<DiracSpinor> *basis,
                           std::vector<DiracSpinor> *basis_positron,
                           const std::vector<DiracSpinor> &spl_basis,
                           const int kappa, const int max_n,
                           const LinAlg::Vector &e_values,
                           const LinAlg::SqMatrix &e_vectors,
                           const Wavefunction &wf);
} // namespace SplineBasis
