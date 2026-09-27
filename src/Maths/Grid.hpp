#pragma once
#include <string>
#include <vector>

//==============================================================================

//! Grid type: determines how r is distributed along the uniform u grid
enum class GridType {
  //! Log-linear: logarithmic for small r, linear for large r (default)
  loglinear,
  //! Purely logarithmic: r = r0 * exp(u)
  logarithmic,
  //! Purely linear: r = r0 + u * du
  linear
};

//==============================================================================

/*!
  @brief Parameters used to construct a Grid.
  @details
  Bundles all grid construction parameters. Either specify the number of points
  directly, or set `innum_points = 0` and provide `indu` to have the number of
  points calculated automatically.
*/
struct GridParameters {
  std::size_t num_points;
  double r0;   //!< Minimum grid point
  double rmax; //!< Maximum grid point
  double b; //!< Log-linear turning point (~logarithmic for r<b, linear for r>b)
  GridType type;

  /*!
    @brief Construct GridParameters.
    @param innum_points Number of grid points. If 0, calculated from `indu`.
    @param inr0    Minimum grid point.
    @param inrmax  Maximum grid point.
    @param inb     Log-linear turning point.
    @param intype  Grid type (loglinear, logarithmic, linear).
    @param indu    Uniform step size; only used if `innum_points == 0`.
  */
  GridParameters(std::size_t innum_points = 1, double inr0 = 1.0,
                 double inrmax = 1.0, double inb = 4.0,
                 GridType intype = GridType::loglinear, double indu = 0);

  /*!
    @brief Construct GridParameters with grid type given as a string.
    @param innum_points Number of grid points. If 0, calculated from `indu`.
    @param inr0      Minimum grid point.
    @param inrmax    Maximum grid point.
    @param inb       Log-linear turning point.
    @param str_type  Grid type as string: "loglinear", "logarithmic", or "linear".
    @param indu      Uniform step size; only used if `innum_points == 0`.
  */
  GridParameters(std::size_t innum_points, double inr0, double inrmax,
                 double inb, const std::string &str_type = "loglinear",
                 double indu = 0);

  //! Converts a string ("loglinear", "logarithmic", "linear") to GridType
  static GridType parseType(const std::string &str_type);
  //! Converts a GridType to its string representation
  static std::string parseType(GridType type);
};

//==============================================================================
//==============================================================================

/*!
  @brief Non-uniform radial grid with Jacobian, suitable for atomic structure
  calculations.
  @details
  Defines a mapping from a uniform grid u = {0, du, 2du, ...} to a
  non-uniform grid r(u), along with the Jacobian dr/du and the convenience
  quantity (1/r)(dr/du).

  Three grid types are supported (see `GridType`):
  - **loglinear** (default): logarithmic at small r, linear at large r;
    controlled by the turning point `b`.
  - **logarithmic**: r = r0 * exp(u); good for tightly bound states.
  - **linear**: uniform spacing in r.

  Grid points and Jacobians are stored as `std::vector<double>` and can be
  iterated over directly (const iterators only).
*/
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
  //! Construct grid directly from parameters
  Grid(double in_r0, double in_rmax, std::size_t in_num_points,
       GridType in_gridtype, double in_b = 0);

  //! Construct grid from a GridParameters struct
  Grid(const GridParameters &in);

  //! Minimum (first) grid point r[0]
  auto r0() const { return m_r.front(); }
  //! Minimum (first) grid point r[0]; alias for r0()
  auto front() const { return m_r.front(); }
  //! Maximum (last) grid point r[N-1]
  auto rmax() const { return m_r.back(); }
  //! Maximum (last) grid point r[N-1]; alias for rmax()
  auto back() const { return m_r.back(); }
  //! Number of grid points
  auto num_points() const { return m_r.size(); }
  //! Number of grid points; alias for num_points()
  auto size() const { return num_points(); }
  //! Uniform step size du
  auto du() const { return m_du; }
  //! Grid type (loglinear, logarithmic, or linear)
  auto type() const { return gridtype; }
  //! Log-linear turning point b: roughly logarithmic for r<b, linear for r>b
  auto loglin_b() const { return m_b; }

  //! Full grid vector r
  const std::vector<double> &r() const { return m_r; }
  //! Grid point r[i], with range checking
  auto r(std::size_t i) const { return m_r.at(i); };
  //! Grid point r[i], with range checking; alias for r(i)
  auto at(std::size_t i) const { return m_r.at(i); };
  //! Grid point r[i], with range checking; alias for r(i)
  auto operator()(std::size_t i) const { return r(i); }

  //! Full Jacobian vector dr/du
  const std::vector<double> &drdu() const { return m_drdu; };
  //! Jacobian (dr/du)[i], with range checking
  auto drdu(std::size_t i) const { return m_drdu.at(i); };

  //! Full vector of (1/r)(dr/du)
  const std::vector<double> &drduor() const { return m_drduor; };
  //! (1/r)(dr/du)[i], with range checking
  auto drduor(std::size_t i) const { return m_drduor.at(i); };

  //! Returns the grid index closest to x. If `require_nearest` is false,
  //! returns the index of the largest grid point  x.
  std::size_t getIndex(double x, bool require_nearest = false) const;

  //! Returns a human-readable string of grid parameters
  std::string gridParameters() const;

  //! Returns a vector of r^k for each grid point
  std::vector<double> rpow(double k) const;

  //! Const iterator to first grid point
  auto begin() const { return m_r.cbegin(); }
  //! Const iterator past last grid point
  auto end() const { return m_r.cend(); }
  //! Const iterator to first grid point
  auto cbegin() const { return m_r.cbegin(); }
  //! Const iterator past last grid point
  auto cend() const { return m_r.cend(); }
  //! Const reverse iterator to last grid point
  auto rbegin() const { return m_r.crbegin(); }
  //! Const reverse iterator before first grid point
  auto rend() const { return m_r.crend(); }
  //! Const reverse iterator to last grid point
  auto crbegin() const { return m_r.crbegin(); }
  //! Const reverse iterator before first grid point
  auto crend() const { return m_r.crend(); }

  /*!
    @brief Extends the grid to at least `new_rmax`.
    @details
    This is the only mutating function; use with caution. The extended grid
    is not guaranteed to end exactly at `new_rmax`  it will typically extend
    slightly beyond, as new points must land on the existing uniform u grid.
  */
  void extend_to(double new_rmax);

  //! Returns a GridParameters struct that can be used to reconstruct this grid
  GridParameters params() const {
    return GridParameters{num_points(), m_r0, rmax(), m_b, gridtype, m_du};
  }

  //! Given r0, rmax, and num_points, calculates the uniform step size du
  static double calc_du_from_num_points(double in_r0, double in_rmax,
                                        std::size_t in_num_points,
                                        GridType in_gridtype, double in_b = 0);
  //! Given r0, rmax, and du, calculates the number of grid points
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
