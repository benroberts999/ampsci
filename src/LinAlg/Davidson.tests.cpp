#include "LinAlg/include.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

namespace UnitTest {
// Random real symmetric sparse matrix (fixed seed). Diagonal elements
// spread * (i + u_i) / dim, with u_i uniform in [0, 1): increasing on average.
// A fraction `fill` of the off-diagonal elements are set, uniform in
// [-0.05, 0.05]
LinAlg::SparseSymmetricMatrix<double>
random_sparse_matrix(std::size_t dim, double fill, double spread);
} // namespace UnitTest

//==============================================================================
// Davidson against the dense (LAPACK) solution of the same matrix: random
// sparse matrix (~5% of elements), spread diagonal, coupling strong enough to
// need many iterations (including restarts). Lowest eigenvalues, eigenvectors
// (up to sign), residuals, and orthonormality
TEST_CASE("Davidson: random sparse matrix", "[LinAlg][Davidson]") {

  const std::size_t dim = 1500;
  const int num_wanted = 5;
  const auto A = UnitTest::random_sparse_matrix(dim, 0.05, 3.0);
  const auto A_dense = A.dense();

  const auto [e, V, convergence] = LinAlg::Davidson(A, num_wanted);
  const auto [e0, V0] = LinAlg::symmhEigensystem(A_dense);

  REQUIRE(e.size() == std::size_t(num_wanted));
  REQUIRE(V.rows() == std::size_t(num_wanted));
  REQUIRE(V.cols() == dim);
  REQUIRE(convergence.converged);

  double worst_E = 0.0, worst_overlap = 0.0, worst_residual = 0.0,
         worst_ortho = 0.0;
  for (std::size_t i = 0; i < std::size_t(num_wanted); ++i) {
    worst_E = std::max(worst_E, std::abs(e[i] - e0[i]));
    double overlap = 0.0;
    for (std::size_t j = 0; j < dim; ++j) {
      overlap += V(i, j) * V0(i, j);
    }
    worst_overlap = std::max(worst_overlap, 1.0 - std::abs(overlap));
    // Residual |A v - e v|, with the dense matrix
    double residual2 = 0.0;
    for (std::size_t j = 0; j < dim; ++j) {
      double Av_j = 0.0;
      for (std::size_t k = 0; k < dim; ++k) {
        Av_j += A_dense(j, k) * V(i, k);
      }
      residual2 += std::pow(Av_j - e[i] * V(i, j), 2);
    }
    worst_residual = std::max(worst_residual, std::sqrt(residual2));
    for (std::size_t k = 0; k <= i; ++k) {
      double vv = 0.0;
      for (std::size_t j = 0; j < dim; ++j) {
        vv += V(i, j) * V(k, j);
      }
      worst_ortho = std::max(worst_ortho, std::abs(vv - (i == k ? 1.0 : 0.0)));
    }
  }
  fmt::print("N={}: {} iterations\n", dim, convergence.iterations);
  fmt::print("E vs dense:          max diff = {:.1e}\n", worst_E);
  fmt::print("|<v|v_dense>| vs 1:  max diff = {:.1e}\n", worst_overlap);
  fmt::print("|A v - E v|:         max      = {:.1e} (reported {:.1e})\n",
             worst_residual, convergence.max_residual);
  fmt::print("<v_i|v_k> vs delta:  max diff = {:.1e}\n", worst_ortho);
  REQUIRE(worst_E < 1.0e-12);
  REQUIRE(worst_overlap < 1.0e-12);
  REQUIRE(worst_residual < 1.0e-10);
  REQUIRE(convergence.max_residual < 1.0e-10);
  REQUIRE(worst_ortho < 1.0e-12);
  // Loose bound on the iteration count (44 with the preconditioner and
  // restarts; three times that without the preconditioner)
  REQUIRE(convergence.iterations < 80);
}

//==============================================================================
// Exactly degenerate eigenvalues: two identical, decoupled blocks. The
// solutions need not equal the dense ones, but must span the same
// eigenspaces: each dense eigenvector of the lowest four lies in the span of
// the four Davidson solutions
TEST_CASE("Davidson: degenerate eigenvalues", "[LinAlg][Davidson]") {

  const std::size_t block = 400;
  const auto block_matrix =
    UnitTest::random_sparse_matrix(block, 0.1, 1.0).dense();
  LinAlg::SparseSymmetricMatrix<double> B(2 * block);
  for (std::size_t i = 0; i < block; ++i) {
    std::vector<std::uint32_t> columns, shifted_columns;
    std::vector<double> values;
    for (std::size_t j = 0; j <= i; ++j) {
      if (block_matrix(i, j) != 0.0) {
        columns.push_back(std::uint32_t(j));
        shifted_columns.push_back(std::uint32_t(j + block));
        values.push_back(block_matrix(i, j));
      }
    }
    B.set_row(i, columns, values);
    B.set_row(i + block, shifted_columns, values);
  }

  const auto [e, V, convergence] = LinAlg::Davidson(B, 4);
  const auto [e0, V0] = LinAlg::symmhEigensystem(B.dense());
  REQUIRE(convergence.converged);
  // Degenerate pairs (dense), and not degenerate between pairs
  REQUIRE(std::abs(e0[1] - e0[0]) < 1.0e-12);
  REQUIRE(std::abs(e0[3] - e0[2]) < 1.0e-12);
  REQUIRE(std::abs(e0[2] - e0[1]) > 1.0e-4);

  double worst_E = 0.0, worst_span = 0.0;
  for (std::size_t i = 0; i < 4; ++i) {
    worst_E = std::max(worst_E, std::abs(e[i] - e0[i]));
    // Norm^2 of the projection of dense vector i onto the Davidson span
    double projection2 = 0.0;
    for (std::size_t k = 0; k < 4; ++k) {
      double overlap = 0.0;
      for (std::size_t j = 0; j < 2 * block; ++j) {
        overlap += V(k, j) * V0(i, j);
      }
      projection2 += overlap * overlap;
    }
    worst_span = std::max(worst_span, std::abs(1.0 - projection2));
  }
  fmt::print("N={}: {} iterations\n", 2 * block, convergence.iterations);
  fmt::print("E vs dense:                  max diff = {:.1e}\n", worst_E);
  fmt::print("|P_Davidson v_dense|^2 vs 1: max diff = {:.1e}\n", worst_span);
  REQUIRE(worst_E < 1.0e-12);
  REQUIRE(worst_span < 1.0e-12);
}

//==============================================================================
// Matrix smaller than the initial sub-matrix: the initial subspace holds the
// exact solutions, so the first iteration is converged
TEST_CASE("Davidson: small matrix", "[LinAlg][Davidson][unit]") {

  const auto S = UnitTest::random_sparse_matrix(100, 1.0, 1.0);
  const auto [e, V, convergence] = LinAlg::Davidson(S, 3);
  const auto e0 = LinAlg::symmhEigensystem(S.dense()).first;
  REQUIRE(convergence.converged);
  REQUIRE(convergence.iterations == 1);
  for (std::size_t i = 0; i < 3; ++i) {
    REQUIRE(e[i] == Approx(e0[i]).epsilon(1.0e-13));
  }
}

//==============================================================================
// Iteration limit reached before convergence: reported as unconverged, with
// the solutions of the last iteration
TEST_CASE("Davidson: iteration limit", "[LinAlg][Davidson][unit]") {

  const auto A = UnitTest::random_sparse_matrix(1500, 0.05, 3.0);
  const auto [e, V, convergence] = LinAlg::Davidson(A, 5, 1.0e-10, 2);
  REQUIRE(!convergence.converged);
  REQUIRE(convergence.iterations == 2);
  REQUIRE(convergence.max_residual > 1.0e-10);
  REQUIRE(e.size() == 5);
  REQUIRE(V.rows() == 5);
}

//==============================================================================
LinAlg::SparseSymmetricMatrix<double>
UnitTest::random_sparse_matrix(std::size_t dim, double fill, double spread) {
  std::mt19937_64 rng{123456789};
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  LinAlg::SparseSymmetricMatrix<double> A(dim);
  for (std::size_t i = 0; i < dim; ++i) {
    std::vector<std::uint32_t> columns;
    std::vector<double> values;
    for (std::size_t j = 0; j < i; ++j) {
      if (uniform(rng) < fill) {
        columns.push_back(std::uint32_t(j));
        values.push_back(0.1 * (uniform(rng) - 0.5));
      }
    }
    columns.push_back(std::uint32_t(i));
    values.push_back(spread * (double(i) + uniform(rng)) / double(dim));
    A.set_row(i, columns, values);
  }
  return A;
}
