#include <cmath>
#include <iostream>
#include <vector>

enum class GridType { loglinear, logarithmic, linear };

/*
To do:
-add warnings for -ve grids (exp)
-don't allow 1 pt, or 0 (for exp) etc.

Make a template
*/

//******************************************************************************
struct Grid {

public:
  const double r0;   // Minimum grid value
  const double rmax; // Maximum grid value
  const int ngp;     // Number of Grid Points
  const double du;   // Uniform grid step size

  std::vector<double> r;      // Grid values r[i]
  std::vector<double> drdu;   // Jacobian (dr/du)[i]
  std::vector<double> drduor; // Convinient: (1/r)*(dr/du)[i]

private:
  void form_loglinear_grid(double b);
  void form_logarithmic_grid();
  void form_linear_grid();

public:
  Grid(double in_r0, double in_rmax, int in_ngp,
       GridType gridtype = GridType::loglinear, double b = 4.);

  int findNextIndex(double x) const;
  int findNearestIndex(double x) const;

  void printDetails() const;

  static double calc_du_from_ngp(double in_r0, double in_rmax, int in_ngp,
                                 GridType gridtype, double b = 4.);

  static int calc_ngp_from_du(double in_r0, double in_rmax, double in_du,
                              GridType gridtype, double b = 4.);
};

//******************************************************************************
inline double Grid::calc_du_from_ngp(double in_r0, double in_rmax, int in_ngp,
                                     GridType gridtype, double b) {
  if (in_ngp == 1)
    return 0;
  switch (gridtype) {
  case GridType::loglinear:
    return (in_rmax - in_r0 + b * log(in_rmax / in_r0)) / (in_ngp - 1);
  case GridType::logarithmic:
    return log(in_rmax / in_r0) / (in_ngp - 1);
  case GridType::linear:
    return (in_rmax - in_r0) / (in_ngp - 1);
  }
  std::cerr << "\nFAIL 63 in Grid: wrong type?\n";
  return 1.;
}

//******************************************************************************
inline int Grid::calc_ngp_from_du(double in_r0, double in_rmax, double in_du,
                                  GridType gridtype, double b) {
  switch (gridtype) {
  case GridType::loglinear:
    return int((in_rmax - in_r0 + b * log(in_rmax / in_r0)) / in_du) + 2;
  case GridType::logarithmic:
    return int(log(in_rmax / in_r0) / in_du) + 2;
  case GridType::linear:
    return int((in_rmax - in_r0) / in_du) + 2;
  }
  std::cerr << "\nFAIL 84 in Grid: wrong type?\n";
  return 1;
}

//******************************************************************************
inline Grid::Grid(double in_r0, double in_rmax, int in_ngp, GridType gridtype,
                  double b)
    : r0(in_r0), rmax(in_rmax), ngp(in_ngp),
      du(calc_du_from_ngp(in_r0, in_rmax, in_ngp, gridtype, b)) {

  r.reserve(ngp);
  drdu.reserve(ngp);   // Jacobian:
  drduor.reserve(ngp); //(1/r)*du/dr (just for convinience)

  switch (gridtype) {
  case GridType::loglinear:
    form_loglinear_grid(b);
    break;
  case GridType::logarithmic:
    form_logarithmic_grid();
    break;
  case GridType::linear:
    form_linear_grid();
    break;
  default:
    std::cerr << "\n FAIL 49 in Grid: no grid type?\n";
  }
}

//******************************************************************************
inline int Grid::findNextIndex(double x) const
// Returns index correspoding to given value
// Note: finds NEXT LARGEST grid point (greater then or equal to.)
// Note: this is slow - don't rely on it inside big loops!
// For linear or exponential, faster to use formula.
// But for log-linear, can't
{
  for (int i = 0; i < ngp; i++)
    if (x < r[i])
      return i;
  return ngp - 1;
}

//******************************************************************************
inline int Grid::findNearestIndex(double x) const
// Note: this is slow - don't rely on it inside big loops!
// Returns index correspoding to given value
{
  if (x < r0)
    return 0;
  if (x > rmax)
    return ngp - 1;
  int index = 0;
  for (int i = 0; i < ngp; i++)
    if (x < r[i]) {
      index = i;
      break;
    }
  // we have (in order): r[i-1], x, r[i]
  if (fabs(x - r[index - 1]) < fabs(r[index] - x))
    return index - 1;
  else
    return index;
}

//******************************************************************************
inline void Grid::printDetails() const {
  printf("Grid: pts=%6i du=%7.5f r0=%.1e Rmax=%7.1f\n", ngp, du, r0, rmax);
}

//******************************************************************************
inline void Grid::form_loglinear_grid(double b)
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
inline void Grid::form_logarithmic_grid()
/*
Standard exponential (logarithmic) grid.
Uses:
  dr/du = r0 * exp(u)
  =>  r = r0 * exp(u)
      u = i*du for i=0,1,2,...
*/
{
  for (int i = 0; i < ngp; i++)
    drdu.push_back(r0 * exp(i * du));
  r = drdu;
  drduor.resize(ngp, 1.);
}

//******************************************************************************
inline void Grid::form_linear_grid()
/*
Should only be used for testing usually
*/
{
  for (int i = 0; i < ngp; i++) {
    double tmp_r = r0 + i * du;
    r.push_back(tmp_r);
    drduor.push_back(1. / tmp_r);
  }
  drdu.resize(ngp, 1.);
}
