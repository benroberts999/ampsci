#include <iostream>
class Grid;
//
#include "Dirac/DiracSpinor.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include <algorithm>
#include <fstream>
#include <gsl/gsl_bspline.h>

class BSplines {
public:
  BSplines(int in_n, int in_k, const Grid &in_grid, double in_r0,
           double in_rmax)
      : m_number_n(in_n),                      //
        m_order_k(in_k),                       //
        m_rgrid_ptr(&in_grid),                 //
        m_rmin_index(in_grid.getIndex(in_r0)), //
        m_rmax_index(in_rmax <= 0.0 ? in_grid.num_points - 1
                                    : in_grid.getIndex(in_rmax)) //
  {
    // constuct_knots();
    // construct_splines();
    construct_splines_gsl();
    write_splines();
  }

private:
  const std::size_t m_number_n;
  const std::size_t m_order_k;
  const Grid *const m_rgrid_ptr;
  const std::size_t m_rmin_index;
  const std::size_t m_rmax_index; // inclusive
  std::vector<std::size_t> m_knots_t;

public:
  std::vector<std::vector<double>> m_Bk;

public:
  std::vector<std::size_t> constuct_knots(int k) {
    std::vector<std::size_t> knots_t;
    knots_t.reserve(m_number_n + k);
    for (int i = 0; i < k; ++i) {
      knots_t.push_back(m_rmin_index); // right? or  t = 0.0?
    }
    for (int i = 0; i < m_number_n - k - 1; ++i) {
      auto ti = static_cast<std::size_t>(
          (double(i + m_rmin_index + 1) / double(m_number_n - k)) *
          double(m_rmax_index - 1));
      knots_t.push_back(ti);
    }
    for (int i = 0; i < k; ++i) {
      knots_t.push_back(m_rmax_index); // right? or  t = 0.0?
    }
    print_knots(knots_t);
    return knots_t;
  }

  void print_knots(std::vector<std::size_t> &knots_t) {
    std::cout << "knots: ";
    auto ti_cout = [](auto x) { std::cout << x << " "; };
    auto rti_cout = [&](auto x) { printf("%.1e ", m_rgrid_ptr->r[x]); };
    std::for_each(knots_t.begin(), knots_t.end(), ti_cout);
    std::cout << "\n = ";
    std::for_each(knots_t.begin(), knots_t.end(), rti_cout);
    std::cout << "\n";
  }

  void construct_splines_gsl() {

    const auto nbreak = m_number_n + 2 - m_order_k;

    gsl_bspline_workspace *bw;
    gsl_vector *B;
    bw = gsl_bspline_alloc(m_order_k, nbreak);
    B = gsl_vector_alloc(m_number_n);

    // gsl_bspline_knots_uniform(0.01, 10.0, bw);

    gsl_vector *breakpts;
    breakpts = gsl_vector_alloc(nbreak);
    for (int i = 0; i < nbreak; i++) {
      auto dindex = (double(i + m_rmin_index) / double(nbreak - 1)) *
                    double(m_rmax_index);
      auto r = m_rgrid_ptr->r[int(dindex)];
      std::cout << r << "\n";
      gsl_vector_set(breakpts, i, r);
    }
    gsl_bspline_knots(breakpts, bw);
    std::cout << "\n\n";

    m_Bk.resize(m_number_n);
    for (int j = 0; j < m_number_n; ++j) {
      m_Bk[j].resize(m_rgrid_ptr->num_points);
    }

    for (std::size_t ir = 0; ir < m_rgrid_ptr->num_points; ++ir) {
      auto r = m_rgrid_ptr->r[ir];
      gsl_bspline_eval(r, B, bw);
      for (int j = 0; j < m_number_n; ++j) {
        double Bj = gsl_vector_get(B, j);
        m_Bk[j][ir] = Bj;
      }
    }
  }

  void write_splines() {
    std::ofstream of("out.txt");

    for (std::size_t ir = 0; ir < m_rgrid_ptr->num_points; ++ir) {
      of << m_rgrid_ptr->r[ir] << " ";
      auto sum = 0.0;
      for (const auto &Bi : m_Bk) {
        sum += Bi[ir];
      }
      for (const auto &Bi : m_Bk) {
        of << Bi[ir] << " ";
      }
      of << sum << "\n";
    }
  }
};
