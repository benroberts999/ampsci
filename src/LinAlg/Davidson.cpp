#include "Davidson.hpp"
#include "Matrix.hpp"
#include "Solvers.hpp"
#include "SparseSymmetricMatrix.hpp"
#include "Vector.hpp"
#include "qip/omp.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <gsl/gsl_cblas.h>
#include <numeric>
#include <tuple>
#include <vector>

namespace LinAlg {

//==============================================================================

std::tuple<Vector<double>, Matrix<double>, DavidsonConvergence>
Davidson(const SparseSymmetricMatrix<double> &A, int number, double tolerance,
         int max_iterations) {
  const auto dim = A.rows();
  const auto num_wanted = std::size_t(number);
  assert(number > 0 && num_wanted < dim);

  // The BLAS calls here are small, between the (OpenMP) products with A: a
  // threaded BLAS pool would only contend with those
  const qip::SingleThreadBlas single_thread_blas{};

  // Solutions followed: those requested, plus a few, so a restart keeps the
  // next-lowest solutions in the subspace
  const auto num_followed = std::min(dim, num_wanted + 5);
  // Restart once the subspace is this large: room for several iterations of
  // new vectors (at most num_wanted each) on top of the followed solutions
  const auto max_subspace = std::min(dim, num_followed + 4 * num_wanted + 20);
  // Size of the sub-matrix diagonalised for the initial subspace
  const auto num_initial =
    std::min(dim, std::max(std::size_t(400), 4 * num_followed));
  // A new vector reduced below this (relative) by orthogonalisation adds
  // nothing to the subspace
  constexpr double min_new_norm = 1.0e-8;
  // Smallest magnitude of the diagonal preconditioner denominator
  constexpr double min_denominator = 1.0e-12;

  const auto diagonal = A.diagonal();

  // Orthonormal subspace basis (rows), and A applied to it; the first
  // num_basis rows are in use, of which the last num_new are not yet in
  // A_basis
  Matrix<double> basis(max_subspace, dim);
  Matrix<double> A_basis(max_subspace, dim);
  std::size_t num_basis = 0;
  std::size_t num_new = 0;

  // Initial subspace: lowest eigenvectors of the sub-matrix of the basis
  // vectors with the lowest diagonal elements (ties: lower index first)
  {
    std::vector<std::size_t> lowest(dim);
    std::iota(lowest.begin(), lowest.end(), std::size_t(0));
    std::partial_sort(lowest.begin(), lowest.begin() + long(num_initial),
                      lowest.end(), [&diagonal](std::size_t i, std::size_t j) {
                        return diagonal[i] < diagonal[j] ||
                               (diagonal[i] == diagonal[j] && i < j);
                      });
    lowest.resize(num_initial);
    const auto sub_vectors = symmhEigensystem(A.dense_submatrix(lowest)).second;
    for (std::size_t r = 0; r < num_followed; ++r) {
      for (std::size_t a = 0; a < num_initial; ++a) {
        basis(r, lowest[a]) = sub_vectors(r, a);
      }
    }
    num_new = num_followed;
  }

  const auto dim_int = int(dim);
  Vector<double> theta;
  Matrix<double> ritz_vectors(num_followed, dim);
  Matrix<double> A_ritz_vectors(num_followed, dim);
  Matrix<double> residuals(num_wanted, dim);
  DavidsonConvergence convergence{};

  for (int iteration = 1;; ++iteration) {
    convergence.iterations = iteration;

    // A applied to the new basis vectors
    {
      Matrix<double> new_vectors(num_new, dim);
      std::copy_n(basis[num_basis], num_new * dim, new_vectors.data());
      const auto A_new = new_vectors * A;
      std::copy_n(A_new.data(), num_new * dim, A_basis[num_basis]);
      num_basis += num_new;
      num_new = 0;
    }
    const auto num_basis_int = int(num_basis);
    const auto num_followed_int = int(num_followed);

    // Rayleigh-Ritz: subspace matrix G = V (A V)^T, symmetric up to rounding
    Matrix<double> G(num_basis, num_basis);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, num_basis_int,
                num_basis_int, dim_int, 1.0, basis.data(), dim_int,
                A_basis.data(), dim_int, 0.0, G.data(), num_basis_int);
    for (std::size_t i = 0; i < num_basis; ++i) {
      for (std::size_t j = 0; j < i; ++j) {
        G(i, j) = G(j, i) = 0.5 * (G(i, j) + G(j, i));
      }
    }
    auto [G_values, G_vectors] = symmhEigensystem(std::move(G));
    theta = std::move(G_values);

    // Ritz vectors (and A applied to them) of the followed solutions
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_followed_int,
                dim_int, num_basis_int, 1.0, G_vectors.data(), num_basis_int,
                basis.data(), dim_int, 0.0, ritz_vectors.data(), dim_int);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_followed_int,
                dim_int, num_basis_int, 1.0, G_vectors.data(), num_basis_int,
                A_basis.data(), dim_int, 0.0, A_ritz_vectors.data(), dim_int);

    // Residuals of the requested solutions, r = A v - e v
    std::vector<std::size_t> unconverged;
    convergence.max_residual = 0.0;
    for (std::size_t i = 0; i < num_wanted; ++i) {
      for (std::size_t j = 0; j < dim; ++j) {
        residuals(i, j) = A_ritz_vectors(i, j) - theta[i] * ritz_vectors(i, j);
      }
      const auto norm = cblas_dnrm2(dim_int, residuals[i], 1);
      convergence.max_residual = std::max(convergence.max_residual, norm);
      if (norm >= tolerance) {
        unconverged.push_back(i);
      }
    }
    convergence.converged = unconverged.empty();
    if (convergence.converged || iteration >= max_iterations)
      break;

    // Restart from the followed solutions, if no room for the new vectors
    if (num_basis + unconverged.size() > max_subspace) {
      std::copy_n(ritz_vectors.data(), num_followed * dim, basis.data());
      std::copy_n(A_ritz_vectors.data(), num_followed * dim, A_basis.data());
      num_basis = num_followed;
    }

    // New vectors: preconditioned residuals, t_j = r_j / (e - A_jj),
    // orthogonalised (twice, for numerical stability) to the subspace
    std::vector<double> overlaps(max_subspace);
    for (const auto i : unconverged) {
      const auto row = num_basis + num_new;
      if (row >= max_subspace)
        break;
      double *t = basis[row];
      for (std::size_t j = 0; j < dim; ++j) {
        const auto denominator = theta[i] - diagonal[j];
        t[j] =
          residuals(i, j) / (std::abs(denominator) > min_denominator ?
                               denominator :
                               std::copysign(min_denominator, denominator));
      }
      const auto norm_before = cblas_dnrm2(dim_int, t, 1);
      for (int pass = 0; pass < 2; ++pass) {
        cblas_dgemv(CblasRowMajor, CblasNoTrans, int(row), dim_int, 1.0,
                    basis.data(), dim_int, t, 1, 0.0, overlaps.data(), 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, int(row), dim_int, -1.0,
                    basis.data(), dim_int, overlaps.data(), 1, 1.0, t, 1);
      }
      const auto norm = cblas_dnrm2(dim_int, t, 1);
      if (norm <= min_new_norm * norm_before)
        continue;
      cblas_dscal(dim_int, 1.0 / norm, t, 1);
      ++num_new;
    }
    // No new direction: cannot improve further
    if (num_new == 0)
      break;
  }

  Vector<double> values(num_wanted);
  std::copy_n(theta.data(), num_wanted, values.data());
  Matrix<double> vectors(num_wanted, dim);
  std::copy_n(ritz_vectors.data(), num_wanted * dim, vectors.data());
  return {std::move(values), std::move(vectors), convergence};
}

} // namespace LinAlg
