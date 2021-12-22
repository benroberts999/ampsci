//******************************************************************************
// Implementations:

namespace LinAlg {

//******************************************************************************
// Solves matrix equation Ax=b for x, for known square matrix A and vector b.
template <typename T> Vector<T> solve_Axeqb(Matrix<T> Am, const Vector<T> &b) {
  static_assert(std::is_same_v<T, double> ||
                    std::is_same_v<T, std::complex<double>>,
                "solve_Axeqb only available for Matrix<double> or "
                "Matrix<complex<double>>");
  assert(Am.rows() == b.size());

  Vector<T> x(Am.cols());

  auto Am_gsl = Am.as_gsl_view();
  const auto b_gsl = b.as_gsl_view();
  auto x_gsl = x.as_gsl_view();

  int sLU = 0;
  gsl_permutation *Am_perm =
      gsl_permutation_alloc(std::max(Am.rows(), Am.cols()));
  if constexpr (std::is_same_v<T, double>) {
    gsl_linalg_LU_decomp(&Am_gsl.matrix, Am_perm, &sLU);
    gsl_linalg_LU_solve(&Am_gsl.matrix, Am_perm, &b_gsl.vector, &x_gsl.vector);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    gsl_linalg_complex_LU_decomp(&Am_gsl.matrix, Am_perm, &sLU);
    gsl_linalg_complex_LU_solve(&Am_gsl.matrix, Am_perm, &b_gsl.vector,
                                &x_gsl.vector);
  }
  gsl_permutation_free(Am_perm);

  return x;
}

//*****************************************************************************
// Solves Av = ev for eigenvalues e and eigenvectors v of symmetric/Hermetian
// matrix A. Returns [e,v], where v(i,j) is the jth element of the ith
// eigenvector corresponding to ith eigenvalue, e(i). e is always real.
template <typename T>
std::pair<Vector<double>, Matrix<T>> symmhEigensystem(Matrix<T> A, bool sort) {
  assert(A.rows() == A.cols());

  // Solves Av = ev for eigenvalues e and eigenvectors v

  const auto dim = A.rows();
  // E. values of Hermetian complex matrix are _real_
  auto eigen_vv = std::make_pair(Vector<double>(dim), Matrix<T>{dim});
  auto &[e_values, e_vectors] = eigen_vv;

  // From GSL:
  // This function computes the eigenvalues and eigenvectors of the
  // symmetric (or hermetian) matrix A. The diagonal and lower triangular part
  // of A are destroyed during the computation, but the strict upper triangular
  // part is not referenced. The eigenvalues are stored in the vector eval and
  // are unordered. The corresponding eigenvectors are stored in the columns of
  // the matrix evec. For example, the eigenvector in the first column
  // corresponds to the first eigenvalue. The eigenvectors are guaranteed to be
  // mutually orthogonal and normalised to unit magnitude.
  // https://www.gnu.org/software/gsl/doc/html/eigen.html

  auto A_gsl = A.as_gsl_view();
  auto val_gsl = e_values.as_gsl_view();
  auto vec_gsl = e_vectors.as_gsl_view();
  if constexpr (std::is_same_v<T, double>) {
    gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(dim);
    gsl_eigen_symmv(&A_gsl.matrix, &val_gsl.vector, &vec_gsl.matrix, work);
    gsl_eigen_symmv_free(work);
    if (sort)
      gsl_eigen_symmv_sort(&val_gsl.vector, &vec_gsl.matrix,
                           GSL_EIGEN_SORT_VAL_ASC);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    gsl_eigen_hermv_workspace *work = gsl_eigen_hermv_alloc(dim);
    gsl_eigen_hermv(&A_gsl.matrix, &val_gsl.vector, &vec_gsl.matrix, work);
    gsl_eigen_hermv_free(work);
    if (sort)
      gsl_eigen_hermv_sort(&val_gsl.vector, &vec_gsl.matrix,
                           GSL_EIGEN_SORT_VAL_ASC);
  }

  // Eigensystems using GSL. NOTE: e-vectors are stored in COLUMNS (not rows) of
  // matrix! Therefore, we transpose the matrix (duuumb)
  auto tmp = e_vectors.transpose();
  e_vectors = std::move(tmp);

  return eigen_vv;
}

//******************************************************************************
// Solves Av = eBv for eigenvalues e and eigenvectors v of symmetric/Hermetian
// matrix pair A,B. Returns [e,v], where v(i,j) is the jth element of the ith
// eigenvector corresponding to ith eigenvalue, e(i). e is always real.
template <typename T>
std::pair<Vector<double>, Matrix<T>> symmhEigensystem(Matrix<T> A, Matrix<T> B,
                                                      bool sort) {
  assert(A.rows() == A.cols());
  assert(B.rows() == B.cols());
  assert(A.rows() == B.rows());

  // GSL: This function computes the eigenvalues and eigenvectors of the
  // generalized symmetric-definite (or hermetian) matrix pair (A, B), and
  // stores them in eval and evec respectively. The computed eigenvectors are
  // normalized to have unit magnitude. On output, B contains its Cholesky
  // decomposition and A is destroyed.
  // https://www.gnu.org/software/gsl/doc/html/eigen.html

  const auto dim = A.rows();
  // E. values of Hermetian complex matrix are _real_
  auto eigen_vv = std::make_pair(Vector<double>(dim), Matrix<T>{dim});
  auto &[e_values, e_vectors] = eigen_vv;

  auto A_gsl = A.as_gsl_view();
  auto B_gsl = B.as_gsl_view();
  auto val_gsl = e_values.as_gsl_view();
  auto vec_gsl = e_vectors.as_gsl_view();
  if constexpr (std::is_same_v<T, double>) {
    gsl_eigen_gensymmv_workspace *work = gsl_eigen_gensymmv_alloc(dim);
    gsl_eigen_gensymmv(&A_gsl.matrix, &B_gsl.matrix, &val_gsl.vector,
                       &vec_gsl.matrix, work);
    gsl_eigen_gensymmv_free(work);
    if (sort)
      gsl_eigen_gensymmv_sort(&val_gsl.vector, &vec_gsl.matrix,
                              GSL_EIGEN_SORT_VAL_ASC);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    gsl_eigen_genhermv_workspace *work = gsl_eigen_genhermv_alloc(dim);
    gsl_eigen_genhermv(&A_gsl.matrix, &B_gsl.matrix, &val_gsl.vector,
                       &vec_gsl.matrix, work);
    gsl_eigen_genhermv_free(work);
    if (sort)
      gsl_eigen_genhermv_sort(&val_gsl.vector, &vec_gsl.matrix,
                              GSL_EIGEN_SORT_VAL_ASC);
  }

  // Eigensystems using GSL. NOTE: e-vectors are stored in COLUMNS (not rows) of
  // matrix! Therefore, we transpose the matrix (duuumb)
  auto tmp = e_vectors.transpose();
  e_vectors = std::move(tmp);

  return eigen_vv;
}

//******************************************************************************
// Solves Av = ev for eigenvalues e and eigenvectors v of non-symmetric real
// matrix A. Returns [e,v], where v(i,j) is the jth element of the ith
// eigenvector corresponding to ith eigenvalue, e(i). A must be real, while e
// and v are complex.
template <typename T>
std::pair<Vector<std::complex<double>>, Matrix<std::complex<double>>>
genEigensystem(Matrix<T> A, bool sort) {
  assert(A.rows() == A.cols());
  static_assert(std::is_same_v<T, double>,
                "genEigensystem only for Real matrix");

  // Solves Av = ev for eigenvalues e and eigenvectors v
  const auto dim = A.rows();
  auto eigen_vv = std::make_pair(Vector<std::complex<double>>(dim),
                                 Matrix<std::complex<double>>(dim));
  auto &[e_values, e_vectors] = eigen_vv;

  // GSL: This function computes eigenvalues and right eigenvectors of the
  // n-by-n real nonsymmetric matrix A.
  auto A_gsl = A.as_gsl_view();
  auto val_gsl = e_values.as_gsl_view();
  auto vec_gsl = e_vectors.as_gsl_view();
  gsl_eigen_nonsymmv_workspace *work = gsl_eigen_nonsymmv_alloc(dim);
  gsl_eigen_nonsymmv(&A_gsl.matrix, &val_gsl.vector, &vec_gsl.matrix, work);
  gsl_eigen_nonsymmv_free(work);
  if (sort)
    gsl_eigen_nonsymmv_sort(&val_gsl.vector, &vec_gsl.matrix,
                            GSL_EIGEN_SORT_ABS_ASC);

  // Eigensystems using GSL. NOTE: e-vectors are stored in COLUMNS (not rows) of
  // matrix! Therefore, we transpose the matrix (duuumb)
  auto tmp = e_vectors.transpose();
  e_vectors = std::move(tmp);

  return eigen_vv;
}

} // namespace LinAlg
