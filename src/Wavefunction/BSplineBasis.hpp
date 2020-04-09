#pragma once
#include "Maths/LinAlg_MatrixVector.hpp"
#include <string>
#include <utility>
class DiracSpinor;
class Wavefunction;
class Grid;
namespace MBPT {
class CorrelationPotential;
}
namespace IO {
class UserInputBlock;
}

/*!
@brief Constucts of spinor/orbital basis using B-splines
(DKB/Reno/Derevianko-Beloy method)

@details
Uses Maths/Bsplines to forma set of B-spline orbitals (using method from [1]).
Diagonalises B-splines over Hamiltonian to produce a set of basis orbitals.
  - [1] K. Beloy and A. Derevianko, Comput. Phys. Commun. 179, 310 (2008).

If \f$\{|i\rangle\}\f$ are the set of \f$2N\f$ DKB spline orbitals (of a given
angular symmetry), and
\f[ H_{ij} = \langle{S_i}|\hat H_{\rm HF}|{S_j}\rangle \,
, \qquad S_{ij} = \langle{S_i|S_j}\rangle. \f]
The eigenvalue problem:
\f[
H_{ij}p_i = \epsilon S_{ij}p_i,
\f]
is solved, yielding \f$2N\f$ eigenvalues \f$\epsilon\f$ with corresponding
eigenvectors \f$p\f$ (half of these are positive energy solutions, half are
negative energy solutions).

Note: form_basis() does not store the eigenvectors, instead, it expands the
basis orbitals and stores them on the regular grid (coordinate space).
i.e., for each eigenvalue, n, the corresponding basis orbital is:
\f[
|{n}\rangle = \sum_i^{2N} p_i |i\rangle\,
\f]
*/
namespace SplineBasis {

struct Parameters {
  Parameters() {}
  Parameters(std::string states, std::size_t n, std::size_t k, double r0,
             double reps, double rmax, bool positronQ);
  Parameters(IO::UserInputBlock input);
  std::string states{};
  std::size_t n{}, k{};
  double r0{}, reps{}, rmax{};
  bool positronQ{false};
};

//! @brief Forms + returns the basis orbitals (expanded in terms of splines)
/*! @details
  - states_str = which states to keep e.g., "25spd10f" (up to n=25 for
  s,p,d-states, and up to n=10 for f states)
  - n_spl: Number of splines (nb: underlying spline set is larger, see [1])
  - k_spl: k order of the B-splines. NB: must have \f$k\geq l_{\rm max}+3\f$ [1]
  - r0_spl: first internal knot (basis orbitals only non-zero after this point)
  - r0_eps: Only calculate splines for r where relative core density larger than
 r0_eps (updates r0 for each l). Typically ~1.0e-8. Set to zero to use r0_spl.
  - rmax_spl: last internal knot (basis orbitals only non-zero before this
  point)
  - wf: Wavefunction object: needed to form Hartree-Fock Hamiltonian
  - positronQ: =true will keep negative energy states (have -ve principal
  quantum number, are appended to end of the basis std::vector). If false,
  discards them.

Note: This function calls the below functions, they rarely need to be called
explicitely, unless you are trying to do something different to usual.
*/
std::vector<DiracSpinor> form_basis(const Parameters &params,
                                    const Wavefunction &wf,
                                    const bool correlationsQ = false);

//! Forms the underlying spline basis (which is not kept)
std::vector<DiracSpinor>
form_spline_basis(const int kappa, const std::size_t n_states,
                  const std::size_t k_spl, const double r0_spl,
                  const double rmax_spl, const Grid &rgrid, const double alpha);

//! Calculates  + reyurns the Hamiltonian \f$H_{ij}\f$ (and \f$S_{ij}\f$)
//! matrices
std::pair<LinAlg::SqMatrix, LinAlg::SqMatrix>
fill_Hamiltonian_matrix(const std::vector<DiracSpinor> &spl_basis,
                        const Wavefunction &wf,
                        const bool correlationsQ = false);

//! @brief Expands basis orbitals in terms of spline orbitals, by diagonalising
//! Hamiltonian
void expand_basis_orbitals(std::vector<DiracSpinor> *basis,
                           std::vector<DiracSpinor> *basis_positron,
                           const std::vector<DiracSpinor> &spl_basis,
                           const int kappa, const int max_n,
                           const LinAlg::Vector &e_values,
                           const LinAlg::SqMatrix &e_vectors,
                           const Wavefunction &wf);
} // namespace SplineBasis
