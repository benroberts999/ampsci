#include <cmath>
#include <iostream>
#include <vector>

enum class GridType { loglinear, logarithmic, linear };

//******************************************************************************
struct Grid {

public: // ? maybe better to be non-const, but private?
  const double r0;
  const double rmax;
  const int Npts;
  const double du;

private:
  std::vector<double> r;
  std::vector<double> drdu;
  std::vector<double> drduor;

  void formLogLinearGrid();
  void formLogGrid();
  void formLinearGrid();

public:
  Grid(double in_r0, double in_rmax, int in_ngp,
       GridType gridtype = GridType::loglinear, double b = 4.);

  const std::vector<double> &refto_r() const;
  const std::vector<double> &refto_drdu() const;
  const std::vector<double> &refto_drduor() const;

  // //??? or just make public?
  // double get_r0() const;
  // double get_rmax() const;
  // int get_ngp() const;
  // double get_du() const;

  void printDetails() const;

  static double calc_du_from_ngp(double in_r0, double in_rmax, int in_ngp,
                                 GridType gridtype, double b = 4.) const;

  static int calc_ngp_from_du(double in_r0, double in_rmax, double in_du,
                              GridType gridtype, double b = 4.) const;
};

//******************************************************************************
Grid::Grid(double in_r0, double in_rmax, int in_ngp, GridType gridtype,
           double b)
    : r0(in_r0), rmax(in_rmax), ngp(in_ngp),
      du(calcdu(in_r0, in_rmax, in_ngp, gridtype, b)) {

  r.reserve(ngp);
  drdu.reserve(ngp);   // Jacobian:
  drduor.reserve(ngp); //(1/r)*du/dr (just for convinience)

  switch (gridtype) {
  case GridType::loglinear:
    formLogLinearGrid();
    break;
  case GridType::logarithmic:
    formLogGrid();
    break;
  case GridType::loglinear:
    logLinearRadialGrid();
    break;
  default:
    std::cerr << "\n FAIL 49 in Grid: no grid type?\n";
  }
}

//******************************************************************************
void Grid::printDetails() const;
{ printf("Grid: pts=%6i du=%7.5f r0=%.1e Rmax=%7.1f\n", ngp, du, r0, rmax); }

//******************************************************************************
static double calc_du_from_ngp(double in_r0, double in_rmax, int in_ngp,
                               GridType gridtype, double b) const {
  switch (gridtype) {
  case GridType::loglinear:
    return (in_rmax - in_r0 + b * log(in_rmax / in_r0)) / (in_ngp - 1);
  case GridType::logarithmic:
    return log(in_rmax / in_r0) / (in_ngp - 1);
  case GridType::linear:
    return (in_rmax - in_r0) / (in_ngp - 1);
  }
}

//******************************************************************************
static int calc_ngp_from_du(double in_r0, double in_rmax, double in_du,
                            GridType gridtype, double b) const {
  switch (gridtype) {
  case GridType::loglinear:
    return int((in_rmax - in_r0 + b * log(in_rmax / in_r0)) / in_du) + 2;
  case GridType::logarithmic:
    return int(log(in_rmax / in_r0) / in_du) + 2;
  case GridType::linear:
    return int((in_rmax - in_r0) / in_du) + 2;
  }
}

//******************************************************************************
void Grid::form_loglinear_grid()
/*
Roughly, grid is logarithmically spaced below r=b, and linear above.
Definition:
  dr/du = r/(b+r)
  => u_0 = r0 + b*ln(r0)
  du = (in_rmax-in_r0+b*log(in_rmax/in_r0))/(in_ngp-1)
  du is constant (step-size for uniformly spaced grid)
Typically (and by default), b = 4 (unit units/bohr radius)
*/
{

  // initial points:
  r.push_back(r0);
  drduor.push_back(1. / (b + r0));
  drdu.push_back(drduor[0] * r0);

  // Use iterative method from Dzuba code to calculate r grid
  double u = r0 + b * log(r0);
  for (int i = 1; i < ngp; i++) {
    u += du;
    double r_tmp = r[i - 1];
    // Integrate dr/dt to find r:
    double delta_r = 1.;
    int ii = 0; // to count number of iterations
    while (delta_r > (r0 * 1.e-15)) {
      double delta_u = u - (r_tmp + b * log(r_tmp));
      double drdu_tmp = r_tmp / (r_tmp + b);
      delta_r = delta_u * drdu_tmp;
      r_tmp += delta_r;
      ii++;
      if (ii > 30)
        break; // usually converges in ~ 2 or 3 steps!
    }
    r.push_back(r_tmp);
    drduor.push_back(1. / (b + r_tmp));
    drdu.push_back(drduor[i] * r_tmp);
  }
}

//******************************************************************************
void Grid::form_logarithmic_grid()
/*
Standard exponential (logarithmic) grid.
Uses:
  dr/du = r0 * exp(u)
  =>  r = r0 * exp(u)
      u = i*du for i=0,1,2,...
*/
{
  for (int i = 0; i < ngp; i++)
    drdu.push_back(r0 * exp(i * h));
  r = drdu;
  drduor.reize(ngp, 1.);
}

//******************************************************************************
void Grid::form_linear_grid()
/*
Should only be used for testing usually
*/
{
  for (int i = 0; i < ngp; i++) {
    double tmp_r = r0 + i * du;
    r.push_back(tmp_r);
    drduor.push_back(1. / tmp_r);
  }
  drdu.reize(ngp, 1.);
}
