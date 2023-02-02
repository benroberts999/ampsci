#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <tuple>
#include <type_traits>

//! Contains classes and functions which use general N-step Adams Moulton method
//! to solve systems of 2x2 ODEs, up to N=12.
namespace AdamsMoulton {

//==============================================================================

//! Pure-virtual struct, holds the derivative matrix for 2x2 system of ODEs.
//! Derive from this, and implement a(t),b(t),c(t),d(t) to define the 2x2
//! ODE.
/*! @details
This struct is used by ODESolver2D, and defined the Derivative matrix, and the
(optional) inhomogenous term.

The system of ODEs is defined such that:

\f[ \frac{dF(t)}{dt} = D(t) * F(t) + S(t) \f]

Where F is a 2D set of functions:

\f[
  F(t) = \begin{pmatrix}
    f(t)\\
    g(t)
  \end{pmatrix},
\f]

D is the 2x2 "derivative matrix":

\f[
  D(t) = \begin{pmatrix}
    a(t) & b(t)\\
    c(t) & d(t)
  \end{pmatrix},
\f]

and S(t) is the (optional) 2D inhomogenous term:

\f[
  S(t) = \begin{pmatrix}
    s_f(t)\\
    s_g(t)
  \end{pmatrix}.
\f]

D (and, optionally S) must be provided by the user by implementing this
DerivativeMatrix. \n
D is defined by four functions, a,b,c,d. These four functions
must be overriden with definitions to define the ODE system.

Template parameters, T and Y

Each of the four functions is a function of t, which has type T.
Usually, T=double or complex<double>, but it may also be an index type (e.g.,
std::size_t) if the derivative matrix is only known numerically at certain grid
points/stored in an array.

The return type, Y, is also a template parameter.
It is usually double, but may be any type (e.g., float or complex<double>).
DerivativeMatrix will work with any type, though only arithmetic types will
work in the ODE solver.

See documentation of ODESolver2D for example.

\par

Note: in tests, deriving from a struct was significantly faster than using
std::function, slightly faster than using function pointers, and performed
about equally to directly implementing the DerivativeMatrix.
*/
template <typename T = double, typename Y = double> struct DerivativeMatrix {
  //! a,b,c,d are derivative matrix functions; all must be user implemented
  virtual Y a(T t) const = 0;
  virtual Y b(T t) const = 0;
  virtual Y c(T t) const = 0;
  virtual Y d(T t) const = 0;
  //! Sf and Sg are optional inhomogenous terms
  virtual Y Sf(T) const { return Y(0); };
  virtual Y Sg(T) const { return Y(0); };
  virtual ~DerivativeMatrix() = default;
};

//==============================================================================

// User-defined type-trait: Checks whether T is a std::complex type
template <typename T> struct is_complex : std::false_type {};
// User-defined type-trait: Checks whether T is a std::complex type
template <typename T> struct is_complex<std::complex<T>> : std::true_type {};
//! User-defined type-trait: Checks whether T is a std::complex type
/*! @details
  For example:
  static_assert(!is_complex_v<double>);
  static_assert(!is_complex_v<float>);
  static_assert(is_complex_v<std::complex<double>>);
  static_assert(is_complex_v<std::complex<float>>);
  static_assert(is_complex<std::complex<float>>::value);
*/
template <typename T> constexpr bool is_complex_v = is_complex<T>::value;

//==============================================================================

//! Inner product of two std::arrays.
/*! @details
\f[ inner_product(a, b) = \sum_{i=0}^{N-1} a_i * b_i \f]
Where `N = min(a.size(), b.size())`.
The types of the arrays may be different (T and U).
However, U must be convertable to T; the return-type is T (same as first
array).
*/
template <typename T, typename U, std::size_t N, std::size_t M>
constexpr T inner_product(const std::array<T, N> &a,
                          const std::array<U, M> &b) {
  static_assert(std::is_convertible_v<U, T>,
                "In inner_product, type of second array (U) must be "
                "convertable to that of dirst (T)");
  constexpr std::size_t Size = std::min(N, M);

  if constexpr (Size == 0) {
    return T{0};
  } else if constexpr (!std::is_same_v<T, U> && is_complex_v<T>) {
    // This is to avoid float conversion warning in case that U=double,
    // T=complex<float>; want to case b to float, then to complex<float>
    T res = a[0] * static_cast<typename T::value_type>(b[0]);
    for (std::size_t i = 1; i < Size; ++i) {
      res += a[i] * static_cast<typename T::value_type>(b[i]);
    }
    return res;
  } else {
    T res = a[0] * static_cast<T>(b[0]);
    for (std::size_t i = 1; i < Size; ++i) {
      res += a[i] * static_cast<T>(b[i]);
    }
    return res;
  }
}

//==============================================================================

namespace helper {

// Simple struct for storing "Raw" Adams ("B") coefficients.
/*
Stored as integers, with 'denominator' factored out. \n
Converted to double  ("A" coefs) once, at compile time (see below). \n
Adams coefficients, a_k, defined such that: \n
 \f[ F_{n+K} = F_{n+K-1} + dx * \sum_{k=0}^K a_k y_{n+k} \f]
where:
 \f[ y = d(F)/dr \f]
Note: the 'order' of the coefs is reversed compared to some sources.
The final coefficient is separated, such that: \n
 \f[ a_k = b_k / denom \f]
for k = {0,1,...,K-1} \n
and
 \f[ a_K = b_K / denom \f]
*/
template <std::size_t K> struct AdamsB {
  long denom;
  std::array<long, K> bk;
  long bK;
};

// List of Adams coefficients data
/*
Note: there is (of course) no 0-step Adams method.
The entry at [0] is invalid, and will produce 0.
It is included so that the index of this list matches the order of the method.
Program will not compile (static_asser) is [0] is requested.
Note: assumes that the kth element corresponds to k-order AM method.
*/
static constexpr auto ADAMS_data = std::tuple{
    AdamsB<0>{1, {}, 0}, // invalid entry, but want index to match order
    AdamsB<1>{2, {1}, 1},
    AdamsB<2>{12, {-1, 8}, 5},
    AdamsB<3>{24, {1, -5, 19}, 9},
    AdamsB<4>{720, {-19, 106, -264, 646}, 251},
    AdamsB<5>{1440, {27, -173, 482, -798, 1427}, 475},
    AdamsB<6>{60480, {-863, 6312, -20211, 37504, -46461, 65112}, 19087},
    AdamsB<7>{
        120960, {1375, -11351, 41499, -88547, 123133, -121797, 139849}, 36799},
    AdamsB<8>{3628800,
              {-33953, 312874, -1291214, 3146338, -5033120, 5595358, -4604594,
               4467094},
              1070017},
    AdamsB<9>{7257600,
              {57281, -583435, 2687864, -7394032, 13510082, -17283646, 16002320,
               -11271304, 9449717},
              2082753},
    AdamsB<10>{479001600,
               {-3250433, 36284876, -184776195, 567450984, -1170597042,
                1710774528, -1823311566, 1446205080, -890175549, 656185652},
               134211265},
    AdamsB<11>{958003200,
               {5675265, -68928781, 384709327, -1305971115, 3007739418,
                -4963166514, 6043521486, -5519460582, 3828828885, -2092490673,
                1374799219},
               262747265},
    AdamsB<12>{2615348736000,
               {-13695779093, 179842822566, -1092096992268, 4063327863170,
                -10344711794985, 19058185652796, -26204344465152,
                27345870698436, -21847538039895, 13465774256510, -6616420957428,
                3917551216986},
               703604254357}};

} // namespace helper

//==============================================================================

//! Stores maximum K (order of AM method) for which we have coefficients
//! implemented.
static constexpr std::size_t K_max =
    std::tuple_size_v<decltype(helper::ADAMS_data)> - 1;

//==============================================================================

//! Holds the K+1 Adams-Moulton ak coefficients for the K-step AM method. Final
//! one, aK, is stored separately.
/*! @details
The Adams coefficients, a_k, are defined such that: \n
 \f[ F_{n+K} = F_{n+K-1} + dx * \sum_{k=0}^K a_k y_{n+k} \f]
where:
 \f[ y = d(F)/dr \f]
Note: the 'order' of the coefs is reversed compared to some sources. \n
The final coefficient is separated, such that: \n
 \f[ a_k = b_k / denom \f]
for k = {0,1,...,K-1} \n
and
 \f[ a_K = b_K / denom \f]
They are stored as doubles regardless of other template parameters.
*/
template <std::size_t K, typename = std::enable_if_t<(K > 0)>,
          typename = std::enable_if_t<(K <= K_max)>>
struct AM_Coefs {

private:
  // Forms the (double) ak coefficients, from the raw (int) bk ones
  static constexpr std::array<double, K> make_ak() {
    const auto &am = std::get<K>(helper::ADAMS_data);
    static_assert(
        am.bk.size() == K,
        "Kth Entry in ADAMS_data must correspond to K-order AM method");
    std::array<double, K> tak{};
    for (std::size_t i = 0; i < K; ++i) {
      tak.at(i) = double(am.bk.at(i)) / double(am.denom);
    }
    return tak;
  }

  // Forms the final (double) aK coefficient, from the raw (int) bK one
  static constexpr double make_aK() {
    const auto &am = std::get<K>(helper::ADAMS_data);
    return (double(am.bK) / double(am.denom));
  }

public:
  //! First K coefficients: ak for k={0,1,...,K-1}
  static constexpr std::array<double, K> ak{make_ak()};
  //! Final aK coefficients: ak for k=K
  static constexpr double aK{make_aK()};
};

//==============================================================================
//! Solves a 2D system of ODEs using a K-step Adams-Moulton method
/*! @details
The system of ODEs is defined such that:

\f[ \frac{dF(t)}{dt} = D(t) * F(t) + S(t) \f]

Where F is a 2D set of functions:

\f[
  F(t) = \begin{pmatrix}
    f(t)\\
    g(t)
  \end{pmatrix},
\f]

D is the 2x2 "derivative matrix":

\f[
  D(t) = \begin{pmatrix}
    a(t) & b(t)\\
    c(t) & d(t)
  \end{pmatrix},
\f]

and S(t) is the (optional) 2D inhomogenous term:

\f[
  S(t) = \begin{pmatrix}
    s_f(t)\\
    s_g(t)
  \end{pmatrix}.
\f]

See struct `DerivativeMatrix` - which is a pure virtual struct that must be
implmented by the user in order to define the ODE.

The step-size, dt, must remain constant (since it must remain consistant
between the K+1 and previous K points). It may be positive or negative,
however (or complex). It's perfectly possible to have a non-uniform grid - this
just introduces a Jacobian into the Derivative matrix; dt must still be
constant.

The template parameter, T, is the type of the argument of the Derivative
Matrix (i.e., type of `t`). This is often `double` or `complex<double>`, but may
also be an index type (e.g., std::size_t) if the derivative matrix is only known
numerically at certain grid points/stored in an array.

The template parameter, Y, is the type of the function value F(t), and the
type of dt, and the return value of the Derivative Matrix. This is often
`double`, but may also be another floating-point type, or std::complex.

The first K points of the function F, and derivative dF/dt, must be known.
You may directly access the f,g (function) and df,dg (derivative) arrays, to
set these points.

Alternatively, you may use the provided function
  ```cpp
  void solve_initial_K(T t0, Y f0, Y g0);
  ```
which automatically sets the first K values for F (and dF), given a single
initial value for F, f0=f(t0), fg=g(t0), by using successive N-step AM
methods, for N={1,2,...,K-1}.

For now, just a 2x2 system. In theory, simple to generalise to an N*N system,
though requires a matrix inversion.

\par

**Example:** Bessel's differential equation

\f[
  x^2 \frac{d^2y}{dx^2} + x \frac{dy}{dx} + (x^2-n^2)y = 0
\f]

With y(0) = 1.0, the solutions are the Bessel functions, y(x) = J_n(x)

This can be re-written as a pair of coupled first-order ODEs:

\f[
  \begin{aligned}
  \frac{dy}{dx} &= p \\
  \frac{dp}{dx} &= \left[\left(\frac{n}{x}\right)^2 - 1\right]y - \frac{1}{x}p
  \end{aligned} 
\f]

Putting this into the notation/form required for the solver we have:

\f[
  F(x) = \begin{pmatrix}
    y(x)\\
    p(x)
  \end{pmatrix}
\f]

with the "derivative matrix":

\f[
  D(x) = \begin{pmatrix}
    0 & 1 \\
    \left(\frac{n}{x}\right)^2 - 1 & \frac{-1}{x}
  \end{pmatrix}
\f]

i.e.,

for n=1:

```cpp
  struct BesselDerivative : AdamsMoulton::DerivativeMatrix<double, double> {
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double t) const final { return 1.0/t/t - 1.0; }
    double d(double t) const final { return -1.0 / t; }
  };
```

Or, more generally (for example):

```cpp
  struct BesselDerivative : AdamsMoulton::DerivativeMatrix<double, double> {
    int n;
    BesselDerivative(int tn) : n(tn) {}
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double t) const final { return std::pow(n/t,2) - 1.0; }
    double d(double t) const final { return -1.0 / t; }
  };
```

Minimal example: -- see full examples included elsewhere

```
  // Construct the Derivative matrix (BesselDerivative defined above) with n=0
  int n = 0;
  BesselDerivative D{n};

  // Set the step size:
  double dt = 0.01;

  // Construct the Solver, using K=6-step method:
  AdamsMoulton::ODESolver2D<6> ode{dt, &D};

  // Since 1/t appears in D, we cannot start at zero. Instead, begin at small t
  double t0 = 1.0e-6;

  // Set initial points:
  // Note: these are *approximate*, since they are technically f(0.0)
  double f0 = 1.0;
  double g0 = 0.0;

  // Use automatic solver for first K points:
  ode.solve_initial_K(t0, f0, g0);

  // Drive forwards another 100 steps
  for (int i = 0; i < 100; ++i) {
    ode.drive();
    std::cout << ode.last_t() << " " << ode.last_f() << '\n';
  }
```

*/
template <std::size_t K, typename T = double, typename Y = double>
class ODESolver2D {
  static_assert(K > 0, "Order (K) for Adams method must be K>0");
  static_assert(K <= K_max,
                "Order (K) requested for Adams method too "
                "large. Adams coefficients are implemented up to K_max-1 only");
  static_assert(
      is_complex_v<Y> || std::is_floating_point_v<Y>,
      "Template parameter Y (function values and dt) must be floating point "
      "or complex");
  static_assert(is_complex_v<Y> || std::is_floating_point_v<Y> ||
                    std::is_integral_v<T>,
                "Template parameter T (derivative matrix argument) must be "
                "floating point, complex, or integral");

private:
  // Stores the AM coefficients
  static constexpr AM_Coefs<K> am{};
  // step size:
  Y m_dt;
  // previous 't' value
  // T m_t{T(0)}; // Should be invalid (nan), but no nan for int
  // Pointer to the derivative matrix
  const DerivativeMatrix<T, Y> *m_D;

public:
  //! Arrays to store the previous K values of f and g.
  /*!
  Note: Stored in the same order, regardless of sign of dt (i.e. whether
  driving forwards or backwards). So f[0] is the 'oldest' value, f[K-1] is
  newest.
  */
  std::array<Y, K> f{}, g{};

  //! Arrays to store the previous K values of derivatives, df and dg
  std::array<Y, K> df{}, dg{};

  //! Array to store the previous K values of t: f.at(i) = f(t.at(i))
  std::array<T, K> t{};

  Y S_scale{1.0};

public:
  //! Construct the solver. dt is the (constant) step size
  //! @param dt -- (constant) step size
  //! @param D -- Pointer to derivative matrix
  /*! @details
   The step-size, dt, may be positive (to drive forwards) or negative (to
   drive backwards); it may also be complex. \n
   Note: a pointer to the DerivativeMatrix, tD, is stored. This may not be null,
   and must outlive the ODESolver2D.
  */
  ODESolver2D(Y dt, const DerivativeMatrix<T, Y> *D) : m_dt(dt), m_D(D) {
    assert(dt != Y{0.0} && "Cannot have zero step-size in ODESolver2D");
    assert(D != nullptr && "Cannot have null Derivative Matrix in ODESolver2D");
  }

  //! Returns the AM order (number of steps), K
  constexpr std::size_t K_steps() const { return K; }

  //! Returns most recent f value. Can also access f array directly
  Y last_f() { return f.back(); }
  //! Returns most recent g value. Can also access g array directly
  Y last_g() { return g.back(); }
  //! Returns most recent t value; last_f() := f(last_t())
  T last_t() { return t.back(); }
  //! Returns the step size
  Y dt() { return m_dt; }

  //! Returns derivative, df/dt(t), given f(t),g(t),t
  Y dfdt(Y ft, Y gt, T tt) const {
    return m_D->a(tt) * ft + m_D->b(tt) * gt + S_scale * m_D->Sf(tt);
  }
  //! Returns derivative, dg/dt(t), given f(t),g(t),t
  Y dgdt(Y ft, Y gt, T tt) const {
    return m_D->c(tt) * ft + m_D->d(tt) * gt + S_scale * m_D->Sg(tt);
  }

  //! Drives the DE system to next value, F(t), assuming system has already
  //! been solved for the K previous values {t-K*dt, ..., t-dt}.
  /*! @details
  Note: t must align with previous values: i.e., t_next = last_t + dt \n
  We re-send t to the function to avoid large build-up of small errors. \n
  There's an overload that avoids this, but care should be taken. \n
  For example, the 10,000th point along t grid may not exactly line up with t0
  + 10000*dt, particularly for complicated non-linear grids. \n

  The type of t_next (T) must match type required by DerivativeMatrix. \n
  This is often T=double if analytic formula for derivative is known. \n
  If   we have a numerical approximation to the derivative stored, e.g., on a
  grid, then we may instead have T=std::size_t (or whatever). \n
  */
  void drive(T t_next) {

    // assert (t_next = Approx[next_t(last_t())])

    // Use AM method to determine new values, given previous K values:
    const auto sf = f.back() + m_dt * (inner_product(df, am.ak) +
                                       am.aK * S_scale * m_D->Sf(t_next));
    const auto sg = g.back() + m_dt * (inner_product(dg, am.ak) +
                                       am.aK * S_scale * m_D->Sg(t_next));
    const auto [a, b, c, d] = std::tuple{m_D->a(t_next), m_D->b(t_next),
                                         m_D->c(t_next), m_D->d(t_next)};
    const auto a0 = m_dt * static_cast<Y>(am.aK);
    const auto a02 = a0 * a0;
    const auto det_inv =
        Y{1.0} / (Y{1.0} - (a02 * (b * c - a * d) + a0 * (a + d)));
    const auto fi = (sf - a0 * (d * sf - b * sg)) * det_inv;
    const auto gi = (sg - a0 * (-c * sf + a * sg)) * det_inv;

    // Shift values along. nb: rotate({1,2,3,4}) = {2,3,4,1}
    // We keep track of previous K values in order to determine next value
    std::rotate(f.begin(), f.begin() + 1, f.end());
    std::rotate(g.begin(), g.begin() + 1, g.end());
    std::rotate(df.begin(), df.begin() + 1, df.end());
    std::rotate(dg.begin(), dg.begin() + 1, dg.end());
    std::rotate(t.begin(), t.begin() + 1, t.end());

    // Sets new values:
    t.back() = t_next;
    f.back() = fi;
    g.back() = gi;
    df.back() = dfdt(fi, gi, t_next);
    dg.back() = dgdt(fi, gi, t_next);
  }

  //! Overload of drive(T t) for 'default' case, where next t is defined as
  //! last_t + dt (for arithmetic/complex types), or last_t++/last_t-- for
  //! integral types (grid index).
  void drive() { drive(next_t(last_t())); }

  //! Automatically sets the first K values for F (and dF), given a single
  //! initial value for F, f0=f(t0), fg=g(t0).
  /*! @details
  It does this by using successive N-step AM methods, for N={1,2,...,K-1}. \n
  i.e., \n
  F[0] must be given \n
  F[1] determined using F[0] and a 1-step AM method \n
  F[2] determined using F[0],F[1] and a 2-step AM method \n
  ... \n
  F[K-1] determined using F[0],F[1],...,F[K-2] and a (K-1)-step AM method \n
  */
  void solve_initial_K(T t0, Y f0, Y g0) {
    t.at(0) = t0;
    f.at(0) = f0;
    g.at(0) = g0;
    df.at(0) = dfdt(f0, g0, t0);
    dg.at(0) = dgdt(f0, g0, t0);
    first_k_i<1>(next_t(t0));
  }

private:
  // only used in solve_initial_K(). Use this, because it works for double t,
  // where next value is t+dt, and for integral t, where next value is t++ is
  // driving forward (dt>0), or t-- if driving backwards (dt<0)
  T next_t(T last_t) {
    if constexpr (std::is_integral_v<T>) {
      return (m_dt > 0.0) ? last_t + 1 : last_t - 1;
    } else {
      return last_t + m_dt;
    }
  }

  // Used recursively by solve_initial_K() to find first K points
  template <std::size_t ik> void first_k_i(T t_next) {
    if constexpr (ik >= K) {
      (void)t_next; // suppress unused variable warning on old g++ versions
      return;
    } else {
      constexpr AM_Coefs<ik> ai{};
      // nb: ai.ak is smaller than df; inner_product still works
      const auto sf = f.at(ik - 1) + m_dt * (inner_product(df, ai.ak) +
                                             am.aK * S_scale * m_D->Sf(t_next));
      const auto sg = g.at(ik - 1) + m_dt * (inner_product(dg, ai.ak) +
                                             am.aK * S_scale * m_D->Sg(t_next));
      const auto a0 = m_dt * static_cast<Y>(ai.aK);
      const auto a02 = a0 * a0;
      const auto [a, b, c, d] = std::tuple{m_D->a(t_next), m_D->b(t_next),
                                           m_D->c(t_next), m_D->d(t_next)};
      const auto det_inv =
          Y{1.0} / (Y{1.0} - (a02 * (b * c - a * d) + a0 * (a + d)));
      const auto fi = (sf - a0 * (d * sf - b * sg)) * det_inv;
      const auto gi = (sg - a0 * (-c * sf + a * sg)) * det_inv;
      // Sets new values:
      t.at(ik) = t_next;
      f.at(ik) = fi;
      g.at(ik) = gi;
      df.at(ik) = dfdt(fi, gi, t_next);
      dg.at(ik) = dgdt(fi, gi, t_next);
      // call recursively
      first_k_i<ik + 1>(next_t(t_next));
    }
  }
};
} // namespace AdamsMoulton