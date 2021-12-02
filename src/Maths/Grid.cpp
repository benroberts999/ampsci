#include "Maths/Grid.hpp"
#include "qip/Methods.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

//******************************************************************************
GridParameters::GridParameters(std::size_t innum_points, double inr0,
                               double inrmax, double inb,
                               const std::string &str_type, double indu)
    : num_points(innum_points == 0 ?
                     Grid::calc_num_points_from_du(inr0, inrmax, indu,
                                                   parseType(str_type), inb) :
                     innum_points),
      r0(inr0),
      rmax(inrmax),
      b(inb),
      type(parseType(str_type)) {}
// indu is optional. Only used if innum_points = 0
GridParameters::GridParameters(std::size_t innum_points, double inr0,
                               double inrmax, double inb, GridType intype,
                               double indu)
    : num_points(innum_points == 0 ? Grid::calc_num_points_from_du(
                                         inr0, inrmax, indu, intype, inb) :
                                     innum_points),
      r0(inr0),
      rmax(inrmax),
      b(inb),
      type(intype) {}
//------------------------------------------------------------------------------
GridType GridParameters::parseType(const std::string &str_type) {
  if (str_type == "loglinear")
    return GridType::loglinear;
  if (str_type == "logarithmic")
    return GridType::logarithmic;
  if (str_type == "linear")
    return GridType::linear;
  return GridType::loglinear;
}

std::string GridParameters::parseType(GridType type) {
  switch (type) {
  case GridType::linear:
    return "linear";
  case GridType::logarithmic:
    return "logarithmic";
  case GridType::loglinear:
    return "log-linear";
  }
  return "?GridType?";
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
Grid::Grid(const GridParameters &in)
    : Grid(in.r0, in.rmax, in.num_points, in.type, in.b) {}
//------------------------------------------------------------------------------
Grid::Grid(double in_r0, double in_rmax, std::size_t in_num_points,
           GridType in_gridtype, double in_b)
    : m_r0(in_r0),
      m_du(calc_du_from_num_points(in_r0, in_rmax, in_num_points, in_gridtype,
                                   in_b)),
      gridtype(in_gridtype),
      b(gridtype == GridType::loglinear ? in_b : 0.0),
      m_r(form_r(gridtype, m_r0, in_num_points, m_du, b)),
      m_drduor(form_drduor(gridtype, m_r, b)),
      m_drdu(form_drdu(gridtype, m_r, m_drduor)) {}

//******************************************************************************
std::size_t Grid::getIndex(double x, bool require_nearest) const
// Returns index correspoding to given value
// Note: finds NEXT LARGEST grid point (greater then or equal to.),
// unluess require_nearest=true, when will give closest point.
// For linear or exponential, faster to use formula.
// But for log-linear, can't.
// I don't think this works for "backwards" grids, maybe not negative
// grids either
{
  auto low = std::lower_bound(m_r.begin(), m_r.end(), x);
  auto index = std::size_t(std::distance(m_r.begin(), low));

  if (index >= m_r.size())
    index--;

  if (!require_nearest || index == 0)
    return index;

  // Must resturn /nearest/ index (we have (in order): r[i-1], x, r[i])
  if (std::fabs(x - m_r[index - 1]) < std::fabs(m_r[index] - x))
    return index - 1;
  else
    return index;
}

//******************************************************************************
std::string Grid::gridParameters() const {
  std::ostringstream out;
  switch (gridtype) {
  case GridType::linear:
    out << "Linear ";
    break;
  case GridType::logarithmic:
    out << "Logarithmic ";
    break;
  case GridType::loglinear:
    out << "Log-linear (b=" << b << ") ";
  }
  out << "grid: " << m_r0 << "->" << rmax() << ", N=" << num_points()
      << ", du=" << m_du;
  return out.str();
}

//******************************************************************************
std::vector<double> Grid::rpow(double k) const {
  std::vector<double> rk;
  rk.reserve(m_r.size());
  for (const auto ir : m_r) {
    rk.push_back(std::pow(ir, k));
  }
  return rk;
}

//******************************************************************************
// Extends grid to new_r_max. Note: This is the only non-const function; use
// with caution
void Grid::extend_to(double new_rmax) {

  const auto old_max_r = rmax();
  if (old_max_r >= new_rmax)
    return;

  // Number of points total grid should have, and number of points needed for
  // 'extra' part of grid:
  const auto total_points =
      calc_num_points_from_du(r0(), new_rmax, du(), gridtype, b);
  const auto new_points = total_points - num_points() + 1;

  // Form the 'extra' part of the grid (from old max to new max)
  const auto x_r = form_r(gridtype, old_max_r, new_points, m_du, b);
  const auto x_drduor = form_drduor(gridtype, x_r, b);
  const auto x_drdu = form_drdu(gridtype, x_r, x_drduor);

  // The first point of new grid matches last of old, so drop these
  m_r.pop_back();
  m_drduor.pop_back();
  m_drdu.pop_back();

  // Insert new grids onto end of existing ones
  m_r.insert(m_r.end(), x_r.begin(), x_r.end());
  m_drduor.insert(m_drduor.end(), x_drduor.begin(), x_drduor.end());
  m_drdu.insert(m_drdu.end(), x_drdu.begin(), x_drdu.end());
}

//******************************************************************************
double Grid::next_r_loglin(double b, double u, double r_guess) {
  // Solve eq. u = r + b ln(r) to find r
  // Use Newtons method
  // => f(r) = r + b ln(r) - u = 0
  // dfdr = b(1/r + 1/(b+r))
  const auto f_u = [b, u](double tr) { return tr + b * std::log(tr) - u; };
  const auto dr = 0.1 * r_guess;
  const auto delta_targ = r_guess * 1.0e-18;
  const auto [ri, delta_r] = qip::Newtons(f_u, r_guess, dr, delta_targ, 30);

  if (std::abs(delta_r / ri) > 1.0e-10) {
    std::cerr << "\nWARNING Grid:194: Converge? " << ri << " " << delta_r / ri
              << "\n";
  }
  return ri;
}
//******************************************************************************
std::vector<double> Grid::form_r(const GridType type, const double r0,
                                 const std::size_t num_points, const double du,
                                 const double b) {
  std::vector<double> r;
  r.reserve(num_points);

  if (type == GridType::loglinear) {
    // u = r + b ln(r), du/dr = r/(b+r)
    auto u = r0 + b * std::log(r0);
    r.push_back(r0);
    for (auto i = 1ul; i < num_points; ++i) {
      u += du;
      r.push_back(next_r_loglin(b, u, r.back()));
    }
  } else if (type == GridType::logarithmic) {
    double u = 0.0;
    r.push_back(r0);
    for (auto i = 1ul; i < num_points; ++i) {
      u += du;
      r.push_back(r0 * std::exp(u));
    }
  } else if (type == GridType::linear) {
    double u = 0.0;
    double tr = r0;
    for (auto i = 0ul; i < num_points; ++i) {
      r.push_back(tr);
      u += du;
      tr += du;
    }
  } else {
    std::cerr << "\nFAIL252 in Grid: grid type?\n";
    std::abort();
  }
  return r;
}

//******************************************************************************
std::vector<double> Grid::form_drduor(const GridType type,
                                      const std::vector<double> &in_r,
                                      const double b) {
  std::vector<double> drduor;
  drduor.reserve(in_r.size());
  if (type == GridType::loglinear) {
    for (const auto &r : in_r) {
      drduor.push_back(1.0 / (b + r));
    }
  } else if (type == GridType::logarithmic) {
    drduor.resize(in_r.size(), 1.0);
  } else if (type == GridType::linear) {
    for (const auto &r : in_r) {
      drduor.push_back(1.0 / r);
    }
  } else {
    std::cerr << "\nFAIL274 in Grid: grid type?\n";
    std::abort();
  }
  return drduor;
}
//******************************************************************************
std::vector<double> Grid::form_drdu(const GridType type,
                                    const std::vector<double> &in_r,
                                    const std::vector<double> &in_drduor) {
  std::vector<double> drdu;
  auto size = in_r.size();
  drdu.reserve(size);
  if (type == GridType::loglinear) {
    for (auto i = 0ul; i < size; ++i) {
      drdu.push_back(in_drduor[i] * in_r[i]);
    }
  } else if (type == GridType::logarithmic) {
    drdu = in_r;
  } else if (type == GridType::linear) {
    drdu.resize(size, 1.0);
  } else {
    std::cerr << "\nFAIL293 in Grid: grid type?\n";
    std::abort();
  }
  return drdu;
}

//******************************************************************************
double Grid::calc_du_from_num_points(double r0, double rmax,
                                     std::size_t num_points, GridType gridtype,
                                     double b) {
  if (num_points == 1)
    return 0;
  switch (gridtype) {
  case GridType::loglinear:
    if (b <= 0)
      std::cerr << "\nFAIL57 in Grid: cant have b=0 for log-linear grid!\n";
    return (rmax - r0 + b * std::log(rmax / r0)) / (double(num_points - 1));
  case GridType::logarithmic:
    return std::log(rmax / r0) / double(num_points - 1);
  case GridType::linear:
    return (rmax - r0) / double(num_points - 1);
  }
  std::cerr << "\nFAIL 63 in Grid: wrong type?\n";
  return 1.0;
}

//******************************************************************************
std::size_t Grid::calc_num_points_from_du(double r0, double rmax, double du,
                                          GridType gridtype, double b) {
  switch (gridtype) {
  case GridType::loglinear:
    if (b <= 0)
      std::cerr << "\nFAIL57 in Grid: cant have b=0 for log-linear grid!\n";
    return std::size_t((rmax - r0 + b * std::log(rmax / r0)) / du) + 2;
  case GridType::logarithmic:
    return std::size_t(std::log(rmax / r0) / du) + 2;
  case GridType::linear:
    return std::size_t((rmax - r0) / du) + 2;
  }
  std::cerr << "\nFAIL 84 in Grid: wrong type?\n";
  return 1;
}
