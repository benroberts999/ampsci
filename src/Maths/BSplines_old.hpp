#pragma once
#include "Maths/Grid.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <fstream>
#include <gsl/gsl_bspline.h>
#include <iostream>
#include <numeric>
#include <utility>

//! @brief Uses GSL to generate set of B-splines and their derivatives
/*! @details Splines are generates upon construction; derivates are not.
Uses GSL: https://www.gnu.org/software/gsl/doc/html/bspline.html
*/
class BSplines {
public:
  //! n is number of spines, k is order. r0/rmax are first/last internal knots
  BSplines(std::size_t in_n, std::size_t in_k, const Grid &in_grid,
           double in_r0, double in_rmax)
      : m_number_n(in_n),
        m_order_k(in_k),
        m_rgrid_ptr(&in_grid),
        // min_index must be at least 1
        m_rmin_index(in_grid.getIndex(in_r0) + 1),
        // max_index must be at least 2 below num_points
        m_rmax_index(in_rmax <= 0.0 ? in_grid.num_points() - 2 :
                                      in_grid.getIndex(in_rmax) - 1) {
    if (in_n != 0 && (in_k > in_n || in_k < 2)) {
      std::cerr << "Fail 24 in BSplines: k>n ? " << m_order_k << " "
                << m_number_n << "\n";
      std::abort();
    }

    construct_splines_gsl();

    if (verbose) {
      std::cout << "B-splines: " << m_number_n << " of order " << m_order_k
                << ". " << GridParameters::parseType(m_rgrid_ptr->type())
                << " subgrid ";
      std::cout << m_rgrid_ptr->r()[m_rmin_index] << ","
                << m_rgrid_ptr->r()[m_rmax_index] << "=[";
      std::cout << m_rmin_index << "," << m_rmax_index << "]\n";
      print_knots();
      write_splines();
    }
  }

  BSplines &operator=(const BSplines &) = delete; // copy assignment
  BSplines(const BSplines &) = delete;            // copy constructor

  ~BSplines() {
    gsl_bspline_free(gsl_bspl_work);
    gsl_vector_free(gsl_breakpts_vec);
    gsl_vector_free(gsl_bspl_vec);
    gsl_matrix_free(gsl_bspl_deriv_mat);
  }

private: // data
  const std::size_t m_number_n;
  const std::size_t m_order_k;
  const Grid *const m_rgrid_ptr;
  // m_rmin_index: must be at least 1
  const std::size_t m_rmin_index;
  // Must be at least 2 below grid.num_points() (inclusive + must have 1 past
  // end)
  const std::size_t m_rmax_index; // inclusive

  gsl_bspline_workspace *gsl_bspl_work = nullptr;
  gsl_vector *gsl_breakpts_vec = nullptr;
  gsl_vector *gsl_bspl_vec = nullptr;
  gsl_matrix *gsl_bspl_deriv_mat = nullptr;

  std::vector<double> m_breaks = {};
  std::vector<std::vector<double>> m_Bk = {};
  std::vector<std::vector<double>> m_dBkdr1 = {};
  std::vector<std::vector<double>> m_dBkdr2 = {};
  std::vector<std::pair<std::size_t, std::size_t>> m_ends = {};

  bool verbose = false;

public:
  std::size_t get_n() const { return m_number_n; }
  std::size_t get_k() const { return m_order_k; }
  const Grid &get_grid() const { return *m_rgrid_ptr; }
  //****************************************************************************
  //! returns nth spline
  const std::vector<double> &get_spline(std::size_t n) const {
    return m_Bk[n]; // add bounds-check?
  }
  //! returns first-order derivative of nth spline
  const std::vector<double> &get_spline_deriv(std::size_t n) const {
    return m_dBkdr1[n]; // add bounds-check?
  }
  //! returns second-order derivative of nth spline
  const std::vector<double> &get_spline_deriv2(std::size_t n) const {
    return m_dBkdr2[n]; // add bounds-check?
  }
  //! returns [first,last] non-zero grid point of nth spline ("p0,pinf")
  const std::pair<std::size_t, std::size_t> &get_ends(std::size_t n) const {
    return m_ends[n]; // add bounds-check?
  }

  //****************************************************************************
  //! Writes splines to text file
  void write_splines(const std::string &ofname = "Bspl.txt",
                     bool deriv = false) const {
    std::ofstream of(ofname);

    const auto &splines = (deriv) ? m_dBkdr1 : m_Bk;

    for (std::size_t ir = 0; ir < m_rgrid_ptr->num_points(); ++ir) {
      of << m_rgrid_ptr->r()[ir] << " ";
      auto sum = 0.0;
      for (const auto &Bi : splines) {
        sum += Bi[ir];
      }
      for (const auto &Bi : splines) {
        of << Bi[ir] << " ";
      }
      of << sum << "\n";
    }
  }

  //****************************************************************************
  void print_knots() const {

    std::cout << "Break points:\n ";
    int n_rows_to_print = 9;
    int count = 0;
    for (const auto &t : m_breaks) {
      printf("%.1e ", t);
      ++count;
      if (count % n_rows_to_print == 0)
        std::cout << "\n ";
    }
    if (count % n_rows_to_print != 0)
      std::cout << "\n";

    std::cout << "Knots:\n ";
    count = 0;
    while (count < int(m_number_n + m_order_k)) {
      auto t = gsl_bspl_work->knots->data[count];
      printf("%.1e ", t);
      ++count;
      if (count % n_rows_to_print == 0)
        std::cout << "\n ";
    }
    if (count % n_rows_to_print != 0)
      std::cout << "\n";
  }

  //****************************************************************************
  //****************************************************************************
private:
  //****************************************************************************
  void construct_splines_gsl() {
    setup_splines();
    calculate_splines();
    set_end_points();
  }

  //****************************************************************************
  void setup_splines() {
    const auto nbreak = m_number_n + 2 - m_order_k;

    gsl_bspl_work = gsl_bspline_alloc(m_order_k, nbreak);
    gsl_bspl_vec = gsl_vector_alloc(m_number_n);
    // brakpoints are "internal" knots
    m_breaks = break_points(nbreak);

    gsl_breakpts_vec = gsl_vector_alloc(nbreak);
    for (std::size_t i = 0; i < nbreak; i++) {
      gsl_vector_set(gsl_breakpts_vec, i, m_breaks[i]);
    }
    gsl_bspline_knots(gsl_breakpts_vec, gsl_bspl_work);
  }

  //----------------------------------------------------------------------------
  std::vector<double> break_points(std::size_t nbreaks) const {
    // Breakpoints: (0, r0, ..., rmax)
    // knots are the same, but have k-fold "multiplicity" at the ends
    std::vector<double> breaks;

    // Set breakpoints according to grid (align exactly w/ grid points)
    // Go 1 past beg/end, to ensure f_spl(r) is in spline interp. region
    auto points =
        qip::uniform_range(m_rmin_index - 1, m_rmax_index + 1, nbreaks - 1);
    breaks.reserve(points.size());
    breaks.push_back(0.0);
    for (auto p : points) {
      breaks.push_back(m_rgrid_ptr->r()[p]);
    }

    // Ensure there are at least k grid points between knots:
    // First point may be less than k points from r0
    std::adjacent_difference(points.begin(), points.end(), points.begin());
    const auto min = *std::min_element(points.begin() + 2, points.end());
    if (min <= m_order_k) {
      std::cerr
          << "\nWARNING 233 in BSplines: not enough points between knots for k="
          << m_order_k << ", " << min << "\n";
    }

    return breaks;
  }

  //----------------------------------------------------------------------------
  void calculate_splines() {
    m_Bk.clear();
    m_Bk.resize(m_number_n);
    for (auto &Bk : m_Bk) {
      Bk.resize(m_rgrid_ptr->num_points());
    }

    for (std::size_t ir = 0; ir < m_rgrid_ptr->num_points(); ++ir) {
      if (ir < m_rmin_index || ir > m_rmax_index)
        // if (ir > m_rmax_index)
        continue;
      auto r = m_rgrid_ptr->r()[ir];
      gsl_bspline_eval(r, gsl_bspl_vec, gsl_bspl_work);
      for (std::size_t j = 0; j < m_number_n; ++j) {
        m_Bk[j][ir] = gsl_vector_get(gsl_bspl_vec, j);
      }
    }
  }

public:
  //----------------------------------------------------------------------------
  //! Calculates + stores 1st and 2nd order derivatives
  void derivitate() {

    const auto n_max_deriv = 2ul;

    m_dBkdr1.clear();
    m_dBkdr1.resize(m_number_n);
    for (auto &dBkdr : m_dBkdr1) {
      dBkdr.resize(m_rgrid_ptr->num_points());
    }
    m_dBkdr2 = m_dBkdr1;

    gsl_bspl_deriv_mat = gsl_matrix_alloc(m_number_n, n_max_deriv + 1);
    for (std::size_t ir = 0; ir < m_rgrid_ptr->num_points(); ++ir) {
      if (ir < m_rmin_index || ir > m_rmax_index)
        // if (ir > m_rmax_index)
        continue;
      auto r = m_rgrid_ptr->r()[ir];

      gsl_bspline_deriv_eval(r, n_max_deriv, gsl_bspl_deriv_mat, gsl_bspl_work);

      for (std::size_t j = 0; j < m_number_n; ++j) {
        m_dBkdr1[j][ir] = gsl_matrix_get(gsl_bspl_deriv_mat, j, 1);
      }
      for (std::size_t j = 0; j < m_number_n; ++j) {
        m_dBkdr2[j][ir] = gsl_matrix_get(gsl_bspl_deriv_mat, j, 2);
      }
    }
  }

private:
  //----------------------------------------------------------------------------
  void set_end_points() {
    // stores first/last non-zero points
    // p0 is first non-zero point
    // pinf is first zero point (i.e., last non-zero point + 1)
    m_ends.clear();
    for (const auto &Bk : m_Bk) {
      auto is_nonzero = [](auto a) { return std::abs(a) > 0.0; };
      auto first_nonzero = std::find_if(Bk.begin(), Bk.end(), is_nonzero);
      auto first_zero = std::find_if_not(first_nonzero, Bk.end(), is_nonzero);
      auto p0 =
          static_cast<std::size_t>(std::distance(Bk.begin(), first_nonzero));
      auto pinf =
          static_cast<std::size_t>(std::distance(Bk.begin(), first_zero));
      m_ends.emplace_back(p0, pinf);
    }
  }
};
