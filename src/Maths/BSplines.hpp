#pragma once
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include <algorithm>
#include <fstream>
#include <gsl/gsl_bspline.h>
#include <iostream>
#include <utility>
// class Grid;
// #include "Dirac/DiracSpinor.hpp"
// #include "Maths/NumCalc_quadIntegrate.hpp"

class BSplines {
public:
  BSplines(std::size_t in_n, std::size_t in_k, const Grid &in_grid,
           double in_r0,
           double in_rmax)
      : m_number_n(in_n),                      //
        m_order_k(in_k),                       //
        m_rgrid_ptr(&in_grid),                 //
        m_rmin_index(in_grid.getIndex(in_r0)), //
        m_rmax_index(in_rmax <= 0.0 ? in_grid.num_points - 1
                                    : in_grid.getIndex(in_rmax)) //
  {
    if (in_k > in_n) {
      std::cerr << "Fail 24 in BSplines: k>n ? " << m_order_k << " "
                << m_number_n << "\n";
      // std::abort();
    }
    std::cout << "Constructing " << m_number_n << " B-splines of order "
              << m_order_k << " on: \n "
              << GridParameters::parseType(m_rgrid_ptr->gridtype)
              << " sub-grid: ";
    std::cout << m_rgrid_ptr->r[m_rmin_index] << " -> "
              << m_rgrid_ptr->r[m_rmax_index] << "  (";
    std::cout << m_rmin_index << " -> " << m_rmax_index << ")\n";
    construct_splines_gsl();
    // print_knots();
  }

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
  const std::size_t m_rmin_index;
  const std::size_t m_rmax_index; // inclusive

  gsl_bspline_workspace *gsl_bspl_work;
  gsl_vector *gsl_breakpts_vec;
  gsl_vector *gsl_bspl_vec;
  gsl_matrix *gsl_bspl_deriv_mat;

  std::vector<double> m_knots;
  std::vector<std::vector<double>> m_Bk;
  std::vector<std::vector<double>> m_dBkdr1;
  std::vector<std::vector<double>> m_dBkdr2;
  std::vector<std::pair<std::size_t, std::size_t>> m_ends;

public:
  std::size_t get_n() const { return m_number_n; }
  std::size_t get_k() const { return m_order_k; }
  const Grid &get_grid() const { return *m_rgrid_ptr; }
  //****************************************************************************
  const std::vector<double> &get_spline(std::size_t n) const {
    return m_Bk[n]; // add bounds-check?
  }
  const std::vector<double> &get_spline_deriv(std::size_t n) const {
    return m_dBkdr1[n]; // add bounds-check?
  }
  const std::vector<double> &get_spline_deriv2(std::size_t n) const {
    return m_dBkdr2[n]; // add bounds-check?
  }
  const std::pair<std::size_t, std::size_t> &get_ends(std::size_t n) const {
    return m_ends[n]; // add bounds-check?
  }

  void derivitate() {

    // m_dBkdr.clear();
    // for (auto &Bk : m_Bk) {
    //   m_dBkdr.push_back(
    //       NumCalc::derivative(Bk, m_rgrid_ptr->drdu, m_rgrid_ptr->du,
    //       n_deriv));
    // }
    // Almost correct /\, but small deviations at small r, large r?

    auto n_max_deriv = 2;

    m_dBkdr1.clear();
    m_dBkdr1.resize(m_number_n);
    for (auto &dBkdr : m_dBkdr1) {
      dBkdr.resize(m_rgrid_ptr->num_points);
    }
    m_dBkdr2 = m_dBkdr1;

    gsl_bspl_deriv_mat = gsl_matrix_alloc(m_number_n, n_max_deriv + 1);
    for (std::size_t ir = 0; ir < m_rgrid_ptr->num_points; ++ir) {
      if (ir < m_rmin_index || ir > m_rmax_index)
        continue;
      auto r = m_rgrid_ptr->r[ir];

      gsl_bspline_deriv_eval(r, n_max_deriv, gsl_bspl_deriv_mat, gsl_bspl_work);

      for (std::size_t j = 0; j < m_number_n; ++j) {
        double Bj = gsl_matrix_get(gsl_bspl_deriv_mat, j, 1);
        m_dBkdr1[j][ir] = Bj;
      }

      for (std::size_t j = 0; j < m_number_n; ++j) {
        double Bj = gsl_matrix_get(gsl_bspl_deriv_mat, j, 2);
        m_dBkdr2[j][ir] = Bj;
      }
    }
  }

  //****************************************************************************
  void write_splines(const std::string ofname = "Bspl.txt",
                     bool deriv = false) const {
    std::ofstream of(ofname);

    const auto &splines = (deriv) ? m_dBkdr1 : m_Bk;

    for (std::size_t ir = 0; ir < m_rgrid_ptr->num_points; ++ir) {
      of << m_rgrid_ptr->r[ir] << " ";
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
    int n_rows_to_print = 7;
    int count = 0;
    for (const auto &t : m_knots) {
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
    m_knots = break_points(nbreak);

    gsl_breakpts_vec = gsl_vector_alloc(nbreak);
    for (std::size_t i = 0; i < nbreak; i++) {
      gsl_vector_set(gsl_breakpts_vec, i, m_knots[i]);
    }
    gsl_bspline_knots(gsl_breakpts_vec, gsl_bspl_work);
  }

  //----------------------------------------------------------------------------
  std::vector<double> break_points(std::size_t nbreaks) const {
    std::vector<double> breaks;
    for (std::size_t i = 0; i < nbreaks; i++) {
      auto dindex = (double(i) / double(nbreaks - 1)) *
                    double(m_rmax_index - m_rmin_index);
      auto r = m_rgrid_ptr->r[m_rmin_index + std::size_t(dindex)];
      breaks.push_back(r);
    }
    return breaks;
  }

  //----------------------------------------------------------------------------
  void calculate_splines() {
    m_Bk.clear();
    m_Bk.resize(m_number_n);
    for (auto &Bk : m_Bk) {
      Bk.resize(m_rgrid_ptr->num_points);
    }

    for (std::size_t ir = 0; ir < m_rgrid_ptr->num_points; ++ir) {
      if (ir < m_rmin_index || ir > m_rmax_index)
        continue;
      auto r = m_rgrid_ptr->r[ir];
      gsl_bspline_eval(r, gsl_bspl_vec, gsl_bspl_work);
      for (std::size_t j = 0; j < m_number_n; ++j) {
        double Bj = gsl_vector_get(gsl_bspl_vec, j);
        m_Bk[j][ir] = Bj;
      }
    }
  }

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
