#pragma once
#include <string>
#include <vector>

//==============================================================================
enum class GridType { loglinear, logarithmic, linear };

//! Parmaters used to create Grid
struct GridParameters {
  std::size_t num_points;
  double r0, rmax, b;
  GridType type;

  //! indu is optional. Only used if innum_points = 0 [then, uses du to find N]
  GridParameters(std::size_t innum_points, double inr0, double inrmax,
                 double inb = 4.0, GridType intype = GridType::loglinear,
                 double indu = 0);

  //! indu is optional. Only used if innum_points = 0 [then, uses du to find N]
  GridParameters(std::size_t innum_points, double inr0, double inrmax,
                 double inb = 4.0, const std::string &str_type = "loglinear",
                 double indu = 0);

  static GridType parseType(const std::string &str_type);
  static std::string parseType(GridType type);
};

//==============================================================================
//==============================================================================
//! Holds grid, including type + Jacobian (dr/du)
class Grid {

private:
  // Minimum grid value
  double m_r0;
  // Uniform (u) grid step size
  double m_du;
  // GridType: loglinear, logarithmic, linear [enum]
  GridType gridtype;
  // log-linear grid 'turning point'
  double m_b;
  // Grid values r[i] = r(u)
  std::vector<double> m_r;
  // Convinient: (1/r)*(dr/du)[i]
  std::vector<double> m_drduor;
  // Jacobian (dr/du)[i]
  std::vector<double> m_drdu;

public:
  //! Manual constructor
  Grid(double in_r0, double in_rmax, std::size_t in_num_points,
       GridType in_gridtype, double in_b = 0);

  //! Constructor using GridParameters class
  Grid(const GridParameters &in);

  //! Minium (first) grid point
  auto r0() const { return m_r0; }
  //! Maximum (final) grid point
  auto rmax() const { return m_r.back(); }
  //! Number of grid points
  auto num_points() const { return m_r.size(); }
  //! Linear step size dr = (dr/dr)*du
  auto du() const { return m_du; }
  //! Grid-type (linear, logarithmic, loglinear)
  auto type() const { return gridtype; }
  //! log-linear grid 'turning point' (~ roughly log for r<b, lin for r>b)
  auto loglin_b() const { return m_b; }

  //! Grid points, r
  const std::vector<double> &r() const { return m_r; };
  auto r(std::size_t i) const { return m_r.at(i); };
  //! Jacobian (dr/du)[i]
  const std::vector<double> &drdu() const { return m_drdu; };
  auto drdu(std::size_t i) const { return m_drdu.at(i); };
  //! Convinient: (1/r)*(dr/du)
  const std::vector<double> &drduor() const { return m_drduor; };
  auto drduor(std::size_t i) const { return m_drduor.at(i); };

  //! Given value, returns grid index. if not require_nearest, returns lower
  std::size_t getIndex(double x, bool require_nearest = false) const;

  //! String human-readable grid parameters
  std::string gridParameters() const;

  //! Calculates+returns vector of 1/r
  std::vector<double> rpow(double k) const;

  //! Extends grid to new_r_max. Note: This is the only non-const function; use
  //! with caution
  /*! @details Note: New grid not guarenteed to end exactly at new_rmax; usually
      will go a bit further (might not line up exactly on existing grid)
  */
  void extend_to(double new_rmax);

  //! Returns set of paramets; may be used to construct grid
  GridParameters params() const {
    return GridParameters{num_points(), m_r0, rmax(), m_b, gridtype, m_du};
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

private: // helper fns
  static double next_r_loglin(double b, double u, double r_guess);
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
