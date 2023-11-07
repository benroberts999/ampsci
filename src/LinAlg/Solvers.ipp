//==============================================================================
// Implementations:

extern "C" {

using complex_double = long double;
// https://gcc.gnu.org/onlinedocs/gcc-4.5.4/gfortran/ISO_005fC_005fBINDING.html
// on apple M1, long double is just double
// nb: it actually doesn't matter, since all operations happen on fortran end

void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *,
            int *);
void zheev_(char *, char *, int *, complex_double *, int *, double *,
            complex_double *, int *, double *, int *);
void dsygv_(int *, char *, char *, int *, double *, int *, double *, int *,
            double *, double *, int *, int *);
void zhegv_(int *, char *, char *, int *, complex_double *, int *,
            complex_double *, int *, double *, complex_double *, int *,
            double *, int *);

void dsyevx_(char *, char *, char *, int *, double *, int *, double *, double *,
             int *, int *, double *, int *, double *, double *, int *, double *,
             int *, int *, int *, int *);
}

namespace LinAlg {

//==============================================================================
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

//============================================================================*
// Solves Av = ev for eigenvalues e and eigenvectors v of symmetric/Hermetian
// matrix A. Returns [e,v], where v(i,j) is the jth element of the ith
// eigenvector corresponding to ith eigenvalue, e(i). e is always real.
template <typename T>
std::pair<Vector<double>, Matrix<T>> symmhEigensystem(Matrix<T> A) {
  assert(A.rows() == A.cols());

  // E. values of Hermetian complex matrix are _real_
  auto eigen_vv = std::make_pair(Vector<double>(A.rows()), std::move(A));
  auto &[e_values, e_vectors] = eigen_vv;

  int dim = (int)e_vectors.rows();
  char jobz{'V'};
  char uplo{'U'};
  int info = 0;

  std::vector<T> work(3 * e_vectors.rows());
  int workspace_size = (int)work.size();

  if constexpr (std::is_same_v<T, double>) {

    // For double (symmetric real) matrix
    dsyev_(&jobz, &uplo, &dim, e_vectors.data(), &dim, e_values.data(),
           work.data(), &workspace_size, &info);

  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    // static_assert(sizeof(complex_double) == 16); // LAPACK assumption
    // nb: this isn't actually required

    // For complex double (Hermetian) matrix
    std::vector<double> Rwork(3 * e_vectors.rows() - 2ul);
    complex_double *evec_ptr =
        reinterpret_cast<complex_double *>(e_vectors.data());
    complex_double *work_ptr = reinterpret_cast<complex_double *>(work.data());
    zheev_(&jobz, &uplo, &dim, evec_ptr, &dim, e_values.data(), work_ptr,
           &workspace_size, Rwork.data(), &info);
  }

  if (info != 0) {
    std::cerr << "\nError 91: symmhEigensystem " << info << " " << dim << " "
              << std::endl;
    if (info < 0) {
      std::cerr << "The " << -info << "-th argument had an illegal value\n";
    } else {
      std::cerr << "The algorithm failed to converge; " << info
                << " off-diagonal elements of an intermediate "
                   "tridiagonal form did not converge to zero.\n";
    }
    std::cerr << info << " " << A.size() << "\n";
  }

  return eigen_vv;
}

//============================================================================*
template <typename T>
std::pair<Vector<double>, Matrix<T>> symmhEigensystem(Matrix<T> A, int number) {
  assert(A.rows() == A.cols());
  assert(number < (int)A.rows());

  // E. values of Hermetian complex matrix are _real_
  auto eigen_vv = std::make_pair(Vector<double>(A.rows()),
                                 Matrix<T>(std::size_t(A.rows())));
  auto &[e_values, e_vectors] = eigen_vv;

  int dim = (int)A.rows();
  char jobz{'V'};
  char uplo{'U'};
  char range{'I'};
  double vl{0.0};
  double vu{0.0};
  int il{1};
  int iu{number};
  double abstol{1.0e-12}; //?
  int m{0};
  int info{0};

  std::vector<T> work(8 * A.rows());
  std::vector<int> iwork(5 * A.rows());
  std::vector<int> ifail(A.rows());
  int workspace_size = (int)work.size();

  // For double (symmetric real) matrix
  dsyevx_(&jobz, &range, &uplo, &dim, A.data(), &dim, &vl, &vu, &il, &iu,
          &abstol, &m, e_values.data(), e_vectors.data(), &dim, work.data(),
          &workspace_size, iwork.data(), ifail.data(), &info);

  if (info != 0) {
    std::cerr << "\nError 135: symmhEigensystem " << info << " " << dim << " "
              << std::endl;
  }

  return eigen_vv;
}

//==============================================================================
// Solves Av = eBv for eigenvalues e and eigenvectors v of symmetric/Hermetian
// matrix pair A,B. Returns [e,v], where v(i,j) is the jth element of the ith
// eigenvector corresponding to ith eigenvalue, e(i). e is always real.
template <typename T>
std::pair<Vector<double>, Matrix<T>> symmhEigensystem(Matrix<T> A,
                                                      Matrix<T> B) {
  assert(A.rows() == A.cols());
  assert(B.rows() == B.cols());
  assert(A.rows() == B.rows());

  // E. values of Hermetian complex matrix are _real_
  auto eigen_vv = std::make_pair(Vector<double>(A.rows()), std::move(A));
  auto &[e_values, e_vectors] = eigen_vv;

  int itype = 1;
  int dim = (int)e_vectors.rows();
  char jobz{'V'};
  char uplo{'U'};
  int info = 0;

  std::vector<T> work(3 * e_vectors.rows());
  int workspace_size = (int)work.size();

  if constexpr (std::is_same_v<T, double>) {

    // For double (symmetric real) matrix
    dsygv_(&itype, &jobz, &uplo, &dim, e_vectors.data(), &dim, B.data(), &dim,
           e_values.data(), work.data(), &workspace_size, &info);

  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    // static_assert(sizeof(complex_double) == 16); // LAPACK assumption
    // nb: this isn't actually required

    // For complex double (Hermetian) matrix
    std::vector<double> Rwork(3 * e_vectors.rows() - 2ul);
    complex_double *evec_ptr =
        reinterpret_cast<complex_double *>(e_vectors.data());
    complex_double *Bmat_ptr = reinterpret_cast<complex_double *>(B.data());
    complex_double *work_ptr = reinterpret_cast<complex_double *>(work.data());
    zhegv_(&itype, &jobz, &uplo, &dim, evec_ptr, &dim, Bmat_ptr, &dim,
           e_values.data(), work_ptr, &workspace_size, Rwork.data(), &info);
    e_vectors.conj_in_place(); // ?? why?
  }

  if (info != 0) {
    std::cerr << "\nError 198: symmhEigensystem " << info << " " << dim << " "
              << B.size() << std::endl;
    if (info < 0) {
      std::cerr << "The " << -info << "-th argument had an illegal value\n";
    } else if (info <= dim) {
      std::cerr << "DSYEV failed to converge; " << info
                << " off-diagonal elements of an intermediate tridiagonal form "
                   "did not converge to zero;";
    } else {
      std::cerr << "The leading minor of order " << info
                << " of B is not positive definite. "
                   "The factorization of B could not be completed and no "
                   "eigenvalues or eigenvectors were computed.";
    }
  }

  return eigen_vv;
}

//==============================================================================
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
