#pragma once
#include "Matrix.hpp"
#include "SparseSymmetricMatrix.hpp"
#include "Vector.hpp"
#include <tuple>

namespace LinAlg {

//! Convergence of an iterative eigensolve; see @ref LinAlg::Davidson
struct DavidsonConvergence {
  //! Number of iterations (each applies A to the new subspace vectors)
  int iterations{0};
  //! Largest residual norm, \f$ |A v_i - e_i v_i| \f$, of the solutions
  double max_residual{0.0};
  //! True if all solutions reached the requested tolerance
  bool converged{false};
};

/*!
  @brief Finds the lowest `number` eigenvalues and eigenvectors of a real
  symmetric sparse matrix, iteratively (block Davidson-Liu).

  @details
  For matrices too large to diagonalise directly, when only a few of the
  lowest solutions are required. The matrix is a
  @ref LinAlg::SparseSymmetricMatrix; the eigenvectors are returned as the
  rows of a (dense) @ref LinAlg::Matrix.

  The matrix @p A never has to be factorised or diagonalised directly:
  Instead, it is only multiplied onto vectors, \f$ A v \f$.
  Each iteration applies @p A to at most `number`
  new vectors, and convergence typically takes a few dozen iterations.
  A dense solver instead costs \f$ N^3 \f$ operations and \f$ N^2 \f$ memory,
  for an N x N matrix.

  * Initial subspace: the lowest eigenvectors of the (dense) sub-matrix of the
    basis vectors with the lowest diagonal elements.
  * Each iteration: Rayleigh-Ritz in the subspace; for each solution not yet
    converged, the residual \f$ r = A v - e v \f$, preconditioned by the
    diagonal, \f$ t_j = r_j / (e - A_{jj}) \f$, is orthogonalised and added
    to the subspace.
  * A few more solutions than requested are followed, so the subspace keeps
    the next-lowest solutions; when the subspace is full, it restarts from
    the followed solutions.
  * Converged when every requested solution has residual norm below
    @p tolerance. The eigenvalue error is then of order tolerance^2 / gap,
    and the eigenvector error of order tolerance / gap, where gap is the
    distance to the nearest other eigenvalue.

  The iteration only reaches the part of the space coupled to the initial
  subspace: if the matrix has a decoupled block containing none of the
  initial (lowest-diagonal) basis vectors, its eigenvalues are missed
  without warning.

  * [1] E. R. Davidson, J. Comput. Phys. 17, 87 (1975): the method, with the
    diagonal preconditioner.
  * [2] B. Liu, in Numerical Algorithms in Chemistry: Algebraic Methods,
    edited by C. Moler and I. Shavitt, LBL-8158 (Lawrence Berkeley
    Laboratory, 1978), p. 49: block form, for several solutions at once.
  * [3] A. Stathopoulos and C. Froese Fischer, Comput. Phys. Commun. 79, 268
    (1994): block Davidson with restarts, and DVDSON, their public Fortran
    77 implementation (CPC Program Library, catalogue identifier ACPZ). The
    same method as here, in a different implementation.

  @param A              Real symmetric sparse matrix.
  @param number         Number of lowest solutions; `0 < number < A.rows()`.
  @param tolerance      Residual norm at which a solution is converged.
  @param max_iterations Iteration limit; on reaching it the current
                        (unconverged) solutions are returned. The same
                        happens if no new direction survives
                        orthogonalisation to the subspace.
  @returns  Tuple `{e, V, convergence}`: the `number` lowest eigenvalues
            (ascending); `V(i,j)` is the j-th element of the i-th eigenvector
            (`number` rows, normalised); and the convergence details.
*/
std::tuple<Vector<double>, Matrix<double>, DavidsonConvergence>
Davidson(const SparseSymmetricMatrix<double> &A, int number,
         double tolerance = 1.0e-10, int max_iterations = 200);

} // namespace LinAlg
