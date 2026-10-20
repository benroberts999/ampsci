#include "LinAlg/include.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

//==============================================================================
// Sparse symmetric matrix, against the dense matrix with the same elements:
// element storage (both triangles), diagonal, sub-matrix, and products with
// one and several vectors
TEST_CASE("LinAlg: sparse symmetric matrix", "[LinAlg][sparse][unit]") {

  const std::size_t dim = 300;
  std::mt19937_64 rng{123456789};
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);

  // ~10% of the lower triangle; some diagonal elements left out (zero)
  LinAlg::SparseSymmetricMatrix<double> A_sparse(dim);
  LinAlg::Matrix<double> A_dense(dim, dim);

  std::size_t num_set = 0;
  for (std::size_t i = 0; i < dim; ++i) {
    std::vector<std::uint32_t> columns;
    std::vector<double> values;
    for (std::size_t j = 0; j <= i; ++j) {
      const auto include = j == i ? (i % 7 != 0) : (uniform(rng) > 0.8);
      if (!include)
        continue;
      const auto value = uniform(rng);
      columns.push_back(std::uint32_t(j));
      values.push_back(value);
      A_dense(i, j) = value;
      A_dense(j, i) = value;
    }
    num_set += columns.size();
    A_sparse.set_row(i, columns, values);
  }

  REQUIRE(A_sparse.rows() == dim);
  REQUIRE(A_sparse.cols() == dim);
  REQUIRE(A_sparse.num_stored() == num_set);

  // Full dense matrix
  const auto A_full = A_sparse.dense();
  double worst_dense = 0.0;
  for (std::size_t i = 0; i < dim; ++i) {
    for (std::size_t j = 0; j < dim; ++j) {
      worst_dense =
        std::max(worst_dense, std::abs(A_full(i, j) - A_dense(i, j)));
    }
  }

  // Diagonal
  const auto diag = A_sparse.diagonal();
  double worst_diag = 0.0;
  for (std::size_t i = 0; i < dim; ++i) {
    worst_diag = std::max(worst_diag, std::abs(diag[i] - A_dense(i, i)));
  }
  REQUIRE(diag[0] == 0.0); // not stored

  // Sub-matrix, unordered indices
  const std::vector<std::size_t> indices{17, 3, 250, 0, 121, 64, 299, 8};
  const auto sub = A_sparse.dense_submatrix(indices);
  REQUIRE(sub.rows() == indices.size());
  double worst_sub = 0.0;
  for (std::size_t a = 0; a < indices.size(); ++a) {
    for (std::size_t b = 0; b < indices.size(); ++b) {
      worst_sub = std::max(
        worst_sub, std::abs(sub(a, b) - A_dense(indices[a], indices[b])));
    }
  }

  // Product with several vectors (rows of X), and with one
  const std::size_t num_vectors = 4;
  LinAlg::Matrix<double> X(num_vectors, dim);
  for (auto &x : X) {
    x = uniform(rng);
  }
  const auto Y = X * A_sparse;
  LinAlg::Matrix<double> Y_dense(num_vectors, dim);
  LinAlg::GEMM(X, A_dense, &Y_dense);
  REQUIRE(Y.rows() == num_vectors);
  REQUIRE(Y.cols() == dim);
  double worst_XA = 0.0;
  for (std::size_t r = 0; r < num_vectors; ++r) {
    for (std::size_t i = 0; i < dim; ++i) {
      worst_XA = std::max(worst_XA, std::abs(Y(r, i) - Y_dense(r, i)));
    }
  }

  LinAlg::Vector<double> b(dim);
  for (std::size_t i = 0; i < dim; ++i) {
    b[i] = uniform(rng);
  }
  const auto Ab = A_sparse * b;
  const auto Ab_dense = A_dense * b;
  REQUIRE(Ab.size() == dim);
  double worst_Ab = 0.0;
  for (std::size_t i = 0; i < dim; ++i) {
    worst_Ab = std::max(worst_Ab, std::abs(Ab[i] - Ab_dense[i]));
  }

  // fmt::print("dense() vs elements set:   max diff = {:.1e}\n", worst_dense);
  // fmt::print("diagonal() vs A_ii:        max diff = {:.1e}\n", worst_diag);
  // fmt::print("dense_submatrix() vs A_(ia, ib): max diff = {:.1e}\n", worst_sub);
  // fmt::print("X A_sparse vs dense X A_sparse:          max diff = {:.1e}\n",
  //            worst_XA);
  // fmt::print("A_sparse b vs dense A_sparse b:          max diff = {:.1e}\n",
  //            worst_Ab);
  REQUIRE(worst_dense == 0.0);
  REQUIRE(worst_diag == 0.0);
  REQUIRE(worst_sub == 0.0);
  REQUIRE(worst_XA < 1.0e-13);
  REQUIRE(worst_Ab < 1.0e-13);
}
