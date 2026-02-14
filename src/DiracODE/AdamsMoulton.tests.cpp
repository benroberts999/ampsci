#include "AdamsMoulton.hpp"
#include "catch2/catch.hpp"
#include <cmath>

//------------------------------------------------------------------------------
TEST_CASE("Simple: integrate forwards", "[AdamsMoulton][unit]") {

  // y''(x) = -y(x)
  // y(0)=0, y'(0)=1
  // => y(x) = sin(x)
  struct DM : AdamsMoulton::DerivativeMatrix<double> {
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double) const final { return -1.0; }
    double d(double) const final { return 0.0; }
  };
  DM D;

  double dt = 0.001;
  AdamsMoulton::ODESolver2D<5> ode{dt, &D};

  REQUIRE(ode.dt() == dt);

  double t0 = 0.0;
  double f0 = 0.0;
  double g0 = 1.0;
  ode.solve_initial_K(t0, f0, g0);

  REQUIRE(ode.dfdt(ode.last_f(), ode.last_g(), ode.last_t()) ==
          Approx(D.a(ode.last_t()) * ode.last_f() +
                 D.b(ode.last_t()) * ode.last_g() + D.Sf(ode.last_t())));

  REQUIRE(ode.dgdt(ode.last_f(), ode.last_g(), ode.last_t()) ==
          Approx(D.c(ode.last_t()) * ode.last_f() +
                 D.d(ode.last_t()) * ode.last_g() + D.Sg(ode.last_t())));

  // test the initial points:
  double t = 0.0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    REQUIRE(ode.f[i] == Approx(std::sin(t)).epsilon(1.0e-5));
    REQUIRE(ode.g[i] == Approx(std::cos(t)).epsilon(1.0e-5));
    REQUIRE(ode.t[i] == Approx(t));
    t += ode.dt();
  }
  REQUIRE(ode.last_t() == ode.t.back());
  REQUIRE(ode.last_f() == ode.f.back());
  REQUIRE(ode.last_g() == ode.g.back());

  int num_steps = 1000;
  for (int i = 0; i < num_steps; ++i) {
    ode.drive();
  }
  REQUIRE(ode.last_f() == Approx(std::sin(ode.last_t())).epsilon(1.0e-5));
  REQUIRE(ode.last_g() == Approx(std::cos(ode.last_t())).epsilon(1.0e-5));

  REQUIRE(ode.last_t() ==
          Approx(t0 + (num_steps + (int)ode.K_steps() - 1) * ode.dt()));
}

//------------------------------------------------------------------------------
TEST_CASE("Simple: integrate backwards", "[AdamsMoulton][unit]") {

  // y''(x) = -y(x)
  // y(0)=0, y'(0)=1
  // => y(x) = sin(x)
  struct DM : AdamsMoulton::DerivativeMatrix<double> {
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double) const final { return -1.0; }
    double d(double) const final { return 0.0; }
  };
  DM D;

  double dt = -0.001;
  AdamsMoulton::ODESolver2D<5> ode{dt, &D};

  double t0 = 0.0;
  double f0 = 0.0;
  double g0 = 1.0;
  ode.solve_initial_K(t0, f0, g0);

  // test the initial points:
  double t = 0.0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    REQUIRE(ode.f[i] == Approx(std::sin(t)).epsilon(1.0e-5));
    REQUIRE(ode.g[i] == Approx(std::cos(t)).epsilon(1.0e-5));
    REQUIRE(ode.t[i] == Approx(t));
    t += ode.dt();
  }

  int num_steps = 1000;
  for (int i = 0; i < num_steps; ++i) {
    ode.drive();
  }
  REQUIRE(ode.last_f() == Approx(std::sin(ode.last_t())).epsilon(1.0e-5));
  REQUIRE(ode.last_g() == Approx(std::cos(ode.last_t())).epsilon(1.0e-5));
  REQUIRE(ode.last_t() < t0);

  REQUIRE(ode.last_t() ==
          Approx(t0 + (num_steps + (int)ode.K_steps() - 1) * ode.dt()));
}

//------------------------------------------------------------------------------
TEST_CASE("Simple: use drive(t)", "[AdamsMoulton][unit]") {

  // y''(x) = -y(x)
  // y(0)=0, y'(0)=1
  // => y(x) = sin(x)
  struct DM : AdamsMoulton::DerivativeMatrix<double> {
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double) const final { return -1.0; }
    double d(double) const final { return 0.0; }
  };
  DM D;

  double dt = 0.001;
  AdamsMoulton::ODESolver2D<5> ode{dt, &D};

  REQUIRE(ode.dt() == dt);

  double t0 = 0.0;
  double f0 = 0.0;
  double g0 = 1.0;
  ode.solve_initial_K(t0, f0, g0);

  // test the initial points:
  double t = 0.0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    REQUIRE(ode.f[i] == Approx(std::sin(t)).epsilon(1.0e-5));
    REQUIRE(ode.g[i] == Approx(std::cos(t)).epsilon(1.0e-5));
    REQUIRE(ode.t[i] == Approx(t));
    t += ode.dt();
  }
  REQUIRE(ode.last_f() == ode.f.back());
  REQUIRE(ode.last_g() == ode.g.back());
  REQUIRE(ode.last_t() == ode.t.back());

  int num_steps = 1000;
  for (int i = 0; i < num_steps; ++i) {
    ode.drive(t);
    t += ode.dt();
  }
  REQUIRE(ode.last_f() == Approx(std::sin(ode.last_t())).epsilon(1.0e-5));
  REQUIRE(ode.last_g() == Approx(std::cos(ode.last_t())).epsilon(1.0e-5));

  REQUIRE(ode.last_t() ==
          Approx(t0 + (num_steps + (int)ode.K_steps() - 1) * ode.dt()));
  REQUIRE(ode.last_t() == t - ode.dt());
}

//------------------------------------------------------------------------------
TEST_CASE("Simple: array index argument", "[AdamsMoulton][unit]") {

  // y''(x) = -y(x)
  // y(0)=0, y'(0)=1
  // => y(x) = sin(x)
  struct DM : AdamsMoulton::DerivativeMatrix<int, double> {
    double a(int) const final { return 0.0; }
    double b(int) const final { return 1.0; }
    double c(int) const final { return -1.0; }
    double d(int) const final { return 0.0; }
  };
  DM D;

  double dt = 0.001;
  AdamsMoulton::ODESolver2D<5, int, double> ode{dt, &D};

  REQUIRE(ode.dt() == dt);

  int t0 = 0;
  double f0 = 0.0;
  double g0 = 1.0;
  ode.solve_initial_K(t0, f0, g0);

  // test the initial points:
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    double t = double(i) * ode.dt();
    REQUIRE(ode.f[i] == Approx(std::sin(t)).epsilon(1.0e-5));
    REQUIRE(ode.g[i] == Approx(std::cos(t)).epsilon(1.0e-5));
    REQUIRE(ode.t[i] == int(i));
  }
  REQUIRE(ode.last_f() == ode.f.back());
  REQUIRE(ode.last_g() == ode.g.back());
  REQUIRE(ode.last_t() == ode.t.back());
  REQUIRE(ode.last_t() == ode.K_steps() - 1);

  int num_steps = 1000;
  for (int i = 0; i < num_steps; ++i) {
    ode.drive();
  }
  REQUIRE(ode.last_f() ==
          Approx(std::sin(ode.last_t() * ode.dt())).epsilon(1.0e-5));
  REQUIRE(ode.last_g() ==
          Approx(std::cos(ode.last_t() * ode.dt())).epsilon(1.0e-5));

  REQUIRE(ode.last_t() == (num_steps + (int)ode.K_steps() - 1));
}

//------------------------------------------------------------------------------
template <std::size_t K>
void helper(const AdamsMoulton::DerivativeMatrix<double> &D, double dt) {

  AdamsMoulton::ODESolver2D<K> ode{dt, &D};
  REQUIRE(ode.K_steps() == K);
  REQUIRE(ode.f.size() == K);
  REQUIRE(ode.g.size() == K);
  REQUIRE(ode.t.size() == K);
  REQUIRE(ode.df.size() == K);
  REQUIRE(ode.dg.size() == K);
  double t0 = 0.0;
  double f0 = 0.0;
  double g0 = 1.0;
  ode.solve_initial_K(t0, f0, g0);
  // test the initial points:
  double t = 0.0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    REQUIRE(ode.f[i] == Approx(std::sin(t)).epsilon(1.0e-5));
    REQUIRE(ode.g[i] == Approx(std::cos(t)).epsilon(1.0e-5));
    REQUIRE(ode.t[i] == Approx(t).epsilon(1.0e-5));
    t += ode.dt();
  }

  int num_steps = 1000;
  for (int i = 0; i < num_steps; ++i) {
    ode.drive();
  }
  REQUIRE(ode.last_f() == Approx(std::sin(ode.last_t())).epsilon(1.0e-5));
}

//------------------------------------------------------------------------------
TEST_CASE("Simple: K=[1,12]", "[AdamsMoulton][unit]") {
  // y''(x) = -y(x)
  // y(0)=0, y'(0)=1
  // => y(x) = sin(x)
  struct DM : AdamsMoulton::DerivativeMatrix<double> {
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double) const final { return -1.0; }
    double d(double) const final { return 0.0; }
  };
  DM D;

  double dt = 0.001;

  helper<1>(D, dt);
  helper<2>(D, dt);
  helper<3>(D, dt);
  helper<4>(D, dt);
  helper<5>(D, dt);
  helper<6>(D, dt);
  helper<7>(D, dt);
  helper<8>(D, dt);
  helper<9>(D, dt);
  helper<10>(D, dt);
  helper<11>(D, dt);
  helper<12>(D, dt);
}

//------------------------------------------------------------------------------
TEST_CASE("Complex", "[AdamsMoulton][unit]") {

  /*
  y''(z) = - i*y(z)
  The exact solution, with y(0)=0, y'(0)=1, is:
  y(z) = - (-1)^(3/4) * sin[(-1)^(1/4)*z]
  */

  using namespace std::complex_literals;

  auto y_expected = [](auto x) {
    const auto r2 = std::sqrt(2.0);
    const auto pow1 = -1.0 / r2 + 1.0i / r2; // (-1)^(3/4)
    const auto pow2 = 1.0 / r2 + 1.0i / r2;  // (-1)^(1/4)
    return -pow1 * std::sin(pow2 * x);
  };

  // Define the DerivativeMatrix for the Bessel equation
  struct ComplexDerivative
    : AdamsMoulton::DerivativeMatrix<std::complex<double>,
                                     std::complex<double>> {
    std::complex<double> a(std::complex<double>) const final { return 0.0; }
    std::complex<double> b(std::complex<double>) const final { return 1.0; }
    std::complex<double> c(std::complex<double>) const final { return -1.0i; }
    std::complex<double> d(std::complex<double>) const final { return 0.0; }
  };

  ComplexDerivative D{};

  constexpr std::size_t N_step = 5; // use 5-step method

  using ComplexAdams = AdamsMoulton::ODESolver2D<N_step, std::complex<double>,
                                                 std::complex<double>>;

  std::complex<double> z0 = 0.0;
  std::complex<double> y0 = 0.0;
  std::complex<double> dy0 = 1.0;

  // First, solve along the positive real axis, starting at 0:
  std::complex<double> dz = {1.0e-3, 0.0};
  auto ode = ComplexAdams{dz, &D};
  ode.solve_initial_K(z0, y0, dy0);

  // Then, solve along the positive imaginary axis, starting at 0:
  std::complex<double> idz = {0.0, 1.0e-3};
  auto ode_iz = ComplexAdams{idz, &D};
  ode_iz.solve_initial_K(z0, y0, dy0);

  // Integrate outwards to |z|=10
  double abs_z_max = 10.0;
  while (std::abs(ode.last_t()) < abs_z_max) {
    ode.drive();
    ode_iz.drive();
  }

  REQUIRE(std::abs(ode.last_t()) ==
          Approx(abs_z_max).margin(1.01 * std::abs(dz)));

  REQUIRE(ode.last_f().real() ==
          Approx(y_expected(ode.last_t()).real()).epsilon(1.0e-5));
  REQUIRE(ode.last_f().imag() ==
          Approx(y_expected(ode.last_t()).imag()).epsilon(1.0e-5));

  REQUIRE(ode_iz.last_f().real() ==
          Approx(y_expected(ode_iz.last_t()).real()).epsilon(1.0e-5));
  REQUIRE(ode_iz.last_f().imag() ==
          Approx(y_expected(ode_iz.last_t()).imag()).epsilon(1.0e-5));
}

//------------------------------------------------------------------------------
TEST_CASE("Inhomogenous", "[AdamsMoulton][unit]") {

  // y''(x) = -y(x) + sin(x)
  // The exact solution, with y(0)=0 and y'(0)=1, is:
  // y(x) = 0.5 * [ 3*sin(x) - x*cos(x) ]

  struct DM : AdamsMoulton::DerivativeMatrix<double> {
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double) const final { return -1.0; }
    double d(double) const final { return 0.0; }
    double Sf(double) const final { return 0.0; }
    double Sg(double t) const final { return std::sin(t); }
  };
  DM D;

  double dt = 0.001;
  AdamsMoulton::ODESolver2D<5> ode{dt, &D};

  REQUIRE(ode.dt() == dt);

  double t0 = 0.0;
  double f0 = 0.0;
  double g0 = 1.0;
  ode.solve_initial_K(t0, f0, g0);

  auto y_exact = [](double x) {
    return 0.5 * (3.0 * std::sin(x) - x * std::cos(x));
  };

  REQUIRE(ode.dfdt(ode.last_f(), ode.last_g(), ode.last_t()) ==
          Approx(D.a(ode.last_t()) * ode.last_f() +
                 D.b(ode.last_t()) * ode.last_g() + D.Sf(ode.last_t())));

  REQUIRE(ode.dgdt(ode.last_f(), ode.last_g(), ode.last_t()) ==
          Approx(D.c(ode.last_t()) * ode.last_f() +
                 D.d(ode.last_t()) * ode.last_g() + D.Sg(ode.last_t())));

  // test the initial points:
  double t = 0.0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    REQUIRE(ode.f[i] == Approx(y_exact(t)).epsilon(1.0e-5));
    t += ode.dt();
  }
  REQUIRE(ode.last_f() == ode.f.back());
  REQUIRE(ode.last_g() == ode.g.back());

  int num_steps = 1000;
  for (int i = 0; i < num_steps; ++i) {
    ode.drive();
  }
  REQUIRE(ode.last_f() == Approx(y_exact(ode.last_t())).epsilon(1.0e-5));

  REQUIRE(ode.last_t() ==
          Approx(t0 + (num_steps + (int)ode.K_steps() - 1) * ode.dt()));
}

//------------------------------------------------------------------------------
TEST_CASE("Inhomogenous with scale", "[AdamsMoulton][unit]") {

  // y''(x) = -y(x) + sin(x)
  // The exact solution, with y(0)=0 and y'(0)=1, is:
  // y(x) = 0.5 * [ 3*sin(x) - x*cos(x) ]
  // The inhomog term is scaled by 2x in D, we re-scale back to one in solution!

  struct DM : AdamsMoulton::DerivativeMatrix<double> {
    double a(double) const final { return 0.0; }
    double b(double) const final { return 1.0; }
    double c(double) const final { return -1.0; }
    double d(double) const final { return 0.0; }
    double Sf(double) const final { return 0.0; }
    double Sg(double t) const final { return 2.0 * std::sin(t); }
  };
  DM D;

  double dt = 0.001;
  AdamsMoulton::ODESolver2D<5> ode{dt, &D};
  ode.S_scale = 0.5;

  REQUIRE(ode.dt() == dt);

  double t0 = 0.0;
  double f0 = 0.0;
  double g0 = 1.0;
  ode.solve_initial_K(t0, f0, g0);

  auto y_exact = [](double x) {
    return 0.5 * (3.0 * std::sin(x) - x * std::cos(x));
  };

  REQUIRE(ode.dfdt(ode.last_f(), ode.last_g(), ode.last_t()) ==
          Approx(D.a(ode.last_t()) * ode.last_f() +
                 D.b(ode.last_t()) * ode.last_g() +
                 ode.S_scale * D.Sf(ode.last_t())));

  REQUIRE(ode.dgdt(ode.last_f(), ode.last_g(), ode.last_t()) ==
          Approx(D.c(ode.last_t()) * ode.last_f() +
                 D.d(ode.last_t()) * ode.last_g() +
                 ode.S_scale * D.Sg(ode.last_t())));

  // test the initial points:
  double t = 0.0;
  for (std::size_t i = 0; i < ode.f.size(); ++i) {
    REQUIRE(ode.f[i] == Approx(y_exact(t)).epsilon(1.0e-5));
    t += ode.dt();
  }
  REQUIRE(ode.last_f() == ode.f.back());
  REQUIRE(ode.last_g() == ode.g.back());

  int num_steps = 1000;
  for (int i = 0; i < num_steps; ++i) {
    ode.drive();
  }
  REQUIRE(ode.last_f() == Approx(y_exact(ode.last_t())).epsilon(1.0e-5));

  REQUIRE(ode.last_t() ==
          Approx(t0 + (num_steps + (int)ode.K_steps() - 1) * ode.dt()));
}