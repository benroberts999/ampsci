#include "Maths/Grid.hpp"
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
    : num_points(innum_points == 0
                     ? Grid::calc_num_points_from_du(inr0, inrmax, indu,
                                                     parseType(str_type), inb)
                     : innum_points),
      r0(inr0),
      rmax(inrmax),
      b(inb),
      type(parseType(str_type)) {}
// indu is optional. Only used if innum_points = 0
GridParameters::GridParameters(std::size_t innum_points, double inr0,
                               double inrmax, double inb, GridType intype,
                               double indu)
    : num_points(innum_points == 0 ? Grid::calc_num_points_from_du(
                                         inr0, inrmax, indu, intype, inb)
                                   : innum_points),
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
Grid::Grid(const GridParameters &in)
    : Grid(in.r0, in.rmax, in.num_points, in.type, in.b) {}
//------------------------------------------------------------------------------
Grid::Grid(double in_r0, double in_rmax, std::size_t in_num_points,
           GridType in_gridtype, double in_b)
    : r0(in_r0),
      rmax(in_rmax),
      num_points(in_num_points),
      du(calc_du_from_num_points(in_r0, in_rmax, in_num_points, in_gridtype,
                                 in_b)),
      gridtype(in_gridtype),                           //
      b(gridtype == GridType::loglinear ? in_b : 0.0), //
      r(form_r(gridtype, r0, num_points, du, b)),      //
      drduor(form_drduor(gridtype, r, b)),             //
      drdu(form_drdu(gridtype, r, drduor))             //
{}
//------------------------------------------------------------------------------

//******************************************************************************
Grid::Grid(const Grid &t_gr, const double new_rmax)
    // This initialiser is only used by the derived class "ExtendedGrid"
    // For now, grid is just extended linearly (or truncated)
    // XXX Make it use helper funnction?
    : r0(t_gr.r0), // thing
      rmax([&]() {
        return (new_rmax >= t_gr.rmax) ? new_rmax
                                       : t_gr.r[t_gr.getIndex(new_rmax)];
      }()), // rmax
      num_points([&]() {
        return (new_rmax >= t_gr.rmax)
                   ? t_gr.num_points +
                         std::size_t((new_rmax - t_gr.rmax) / t_gr.du)
                   : t_gr.getIndex(new_rmax);
      }()),                    // num_points
      du(t_gr.du),             // du
      gridtype(t_gr.gridtype), // gridtype
      b(t_gr.b),               // b
      r([&]() {
        if (new_rmax >= t_gr.rmax) {
          return [&]() {
            // extend radial vector:
            auto tmp_r = t_gr.r;
            auto r_i = t_gr.rmax;
            for (auto i = t_gr.num_points; i < this->num_points; i++) {
              r_i += du;
              tmp_r.push_back(r_i);
            }
            return tmp_r;
          }();
        }
        return std::vector<double>(t_gr.r.begin(),
                                   t_gr.r.begin() + long(this->num_points));
      }()), // r
      drduor([&]() {
        std::vector<double> temp_drduor = t_gr.drduor;
        for (auto i = t_gr.num_points; i < this->num_points; i++) {
          temp_drduor.push_back(1.0 / this->r[i]);
        }
        return temp_drduor;
      }()), // drduor
      drdu([&]() {
        if (new_rmax > t_gr.rmax) {
          auto temp_drdu = t_gr.drdu;
          temp_drdu.resize(this->num_points, 1.0);
          return temp_drdu;
        }
        return std::vector<double>(t_gr.drdu.begin(),
                                   t_gr.drdu.begin() + long(this->num_points));
      }()) // drdu
{}

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
  auto low = std::lower_bound(r.begin(), r.end(), x);
  auto index = std::size_t(std::distance(r.begin(), low));

  if (index >= r.size())
    index--;

  if (!require_nearest || index == 0)
    return index;

  // Must resturn /nearest/ index (we have (in order): r[i-1], x, r[i])
  if (std::fabs(x - r[index - 1]) < std::fabs(r[index] - x))
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
  out << "grid: " << r0 << "->" << rmax << ", N=" << num_points
      << ", du=" << du;
  return out.str();
}

//******************************************************************************
std::vector<double> Grid::rpow(double k) const {
  std::vector<double> rk;
  rk.reserve(num_points);
  for (const auto ir : r) {
    rk.push_back(std::pow(ir, k));
  }
  return rk;
}

//******************************************************************************
std::vector<double> Grid::form_r(const GridType type, const double r0,
                                 const std::size_t num_points, const double du,
                                 const double b) {
  std::vector<double> r;
  r.reserve(num_points);

  if (type == GridType::loglinear) {
    // u = r + b ln(r), du/dr = r/(b+r)
    r.push_back(r0);
    auto u = r0 + b * std::log(r0);
    auto r_prev = r0;
    for (auto i = 1ul; i < num_points; i++) {
      u += du;
      double r_tmp = r_prev;
      // Solve eq. u = r + b ln(r) to find r
      double delta_r = 1.0;
      int ii = 0; // to count number of iterations
      while (std::fabs(delta_r) > (r_tmp * 1.0e-17)) {
        double delta_u = u - (r_tmp + b * std::log(r_tmp));
        double drdu_tmp = r_tmp / (r_tmp + b);
        delta_r = delta_u * drdu_tmp;
        r_tmp += delta_r;
        if (++ii > 20) {
          if (std::fabs(delta_r / r_tmp) > 1.0e-6) {
            std::cerr << "WARNING Grid:194: Converge? " << i << " " << r_tmp
                      << " " << delta_r / r_tmp << "\n";
          }
          break; // usually converges in ~ 2 or 3 steps!
        }
      }
      r.push_back(r_tmp);
      r_prev = r_tmp;
    }
  } else if (type == GridType::logarithmic) {
    for (std::size_t i = 0; i < num_points; i++) {
      r.push_back(r0 * std::exp(double(i) * du));
    }
  } else if (type == GridType::linear) {
    for (std::size_t i = 0; i < num_points; i++) {
      double tmp_r = r0 + double(i) * du;
      r.push_back(tmp_r);
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
  return 1.;
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
