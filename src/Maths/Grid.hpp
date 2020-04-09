#pragma once
#include <string>
#include <vector>

// Fix the "extendedGrid to extend using a specific grid-type?"

//******************************************************************************
enum class GridType { loglinear, logarithmic, linear };

//! Parmaters used to create Grid
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
//! Holds grid, including type + Jacobian (dr/du)
class Grid {

public:
  //! Minimum grid value
  const double r0;
  //! Maximum grid value
  const double rmax;
  //! Number of Grid Points
  const std::size_t num_points;
  //! Uniform (u) grid step size
  const double du;

public:
  //! GridType: loglinear, logarithmic, linear [enum]
  const GridType gridtype;
  //! log-linear grid 'turning point'
  const double b;
  //! Grid values r[i] = r(u)
  const std::vector<double> r;
  // Convinient: (1/r)*(dr/du)[i]
  const std::vector<double> drduor;
  //! Jacobian (dr/du)[i]
  const std::vector<double> drdu;

private:
  Grid(double in_r0, double in_rmax, std::size_t in_num_points,
       GridType in_gridtype, double in_b = 0);

public:
  Grid(const GridParameters &in);

  //! Given value, returns grid index. if not require_nearest, returns lower
  std::size_t getIndex(double x, bool require_nearest = false) const;

  //! String human-readable grid parameters
  std::string gridParameters() const;

  GridParameters params() const {
    return GridParameters{num_points, r0, rmax, b, gridtype, du};
    // num_points, r0, rmax, b, grid_type, du_tmp
  }

  // Static functions: can be called outside of instantialised object
  //! Given r0/rmax + num_points, calculates du
  static double calc_du_from_num_points(double in_r0, double in_rmax,
                                        std::size_t in_num_points,
                                        GridType in_gridtype, double in_b = 0);
  //! Given r0/rmax + du, calculates num_points
  static std::size_t calc_num_points_from_du(double in_r0, double in_rmax,
                                             double in_du, GridType in_gridtype,
                                             double in_b = 0);

  //! Calculates+returns vector of 1/r
  std::vector<double> rpow(double k) const;

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
//! @brief Created an extended grid (extending out to new_rmax), based on
//! esisting  Grid. Grid-poins match up to end of smaller.
//! @details Used mainly for continuum states. For now, can only extend using
//! linear method. Could update later (but no real need)
class ExtendedGrid : public Grid {
public:
  ExtendedGrid(const Grid &t_gr, const double new_rmax)
      : Grid(t_gr, new_rmax) {}
};
