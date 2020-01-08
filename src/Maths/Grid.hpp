#pragma once
#include <string>
#include <vector>

// Fix the "extendedGrid to extend using a specific grid-type?"

//******************************************************************************
enum class GridType { loglinear, logarithmic, linear };

struct GridParameters {
  std::size_t num_points;
  double r0, rmax, b;
  GridType type;
  GridParameters(std::size_t innum_points, double inr0, double inrmax,
                 double inb = 4.0, GridType intype = GridType::loglinear,
                 double indu = 0);

  // indu is optional. Only used if innum_points = 0 [then, uses du to find N]
  GridParameters(std::size_t innum_points, double inr0, double inrmax,
                 double inb = 4.0, const std::string &str_type = "loglinear",
                 double indu = 0);

  static GridType parseType(const std::string &str_type);
  static std::string parseType(GridType type);
};

//******************************************************************************
//******************************************************************************
class Grid {

public:
  const double r0;              // Minimum grid value
  const double rmax;            // Maximum grid value
  const std::size_t num_points; // Number of Grid Points
  const double du;              // Uniform grid step size

public:
  const GridType gridtype;
  const double b;
  const std::vector<double> r;      // Grid values r[i]
  const std::vector<double> drduor; // Convinient: (1/r)*(dr/du)[i]
  const std::vector<double> drdu;   // Jacobian (dr/du)[i]

public:
  Grid(double in_r0, double in_rmax, std::size_t in_num_points,
       GridType in_gridtype, double in_b = 0);
  Grid(const GridParameters &in);

  std::size_t getIndex(double x, bool require_nearest = false) const;

  std::string gridParameters() const;

  // Static functions: can be called outside of instantialised object
  static double calc_du_from_num_points(double in_r0, double in_rmax,
                                        std::size_t in_num_points,
                                        GridType in_gridtype, double in_b = 0);
  static std::size_t calc_num_points_from_du(double in_r0, double in_rmax,
                                             double in_du, GridType in_gridtype,
                                             double in_b = 0);

  std::vector<double> inverse_r() const;

protected:
  // secondary constructor: used for derivated class
  Grid(const Grid &t_gr, const double new_rmax);

private: // helper fns
  static std::vector<double> form_r(const GridType type, const double r0,
                                    const std::size_t num_points,
                                    const double du, const double b);
  static std::vector<double> form_drduor(const GridType type,
                                         const std::vector<double> &in_r,
                                         const double b);
  static std::vector<double> form_drdu(const GridType type,
                                       const std::vector<double> &in_r,
                                       const std::vector<double> &in_drduor);
};

//******************************************************************************
class ExtendedGrid : public Grid {
  // For now, can only extend using linear method. Could update later (but no
  // real need)
public:
  ExtendedGrid(const Grid &t_gr, const double new_rmax)
      : Grid(t_gr, new_rmax) {}
};
