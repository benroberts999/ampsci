#include <iostream>
class Grid;
//
#include "Dirac/DiracSpinor.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include <algorithm>
#include <fstream>

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
    constuct_knots();
    construct_splines();
    write_splines();
  }

private:
  const int m_number_n;
  const int m_order_k;
  const Grid *const m_rgrid_ptr;
  const std::size_t m_rmin_index;
  const std::size_t m_rmax_index; // inclusive
  std::vector<std::size_t> m_knots_t;

public:
  std::vector<std::vector<double>> m_Bk;

public:
  void constuct_knots() {
    m_knots_t.reserve(m_number_n + m_order_k);
    for (int i = 0; i < m_order_k; ++i) {
      m_knots_t.push_back(m_rmin_index); // right? or  t = 0.0?
    }
    for (int i = 0; i < m_number_n; ++i) {
      auto ti = static_cast<std::size_t>(
          (double(i + 1) / double(m_number_n + 1)) * double(m_rmax_index - 1));
      m_knots_t.push_back(ti);
    }
    for (int i = 0; i < m_order_k; ++i) {
      m_knots_t.push_back(m_rmax_index); // right? or  t = 0.0?
    }
    print_knots();
  }

  void print_knots() {
    std::cout << "knots: ";
    auto ti_cout = [](auto x) { std::cout << x << " "; };
    auto rti_cout = [&](auto x) { printf("%.1e ", m_rgrid_ptr->r[x]); };
    std::for_each(m_knots_t.begin(), m_knots_t.end(), ti_cout);
    std::cout << "\n = ";
    std::for_each(m_knots_t.begin(), m_knots_t.end(), rti_cout);
    std::cout << "\n";
  }

  void construct_splines() {
    // std::vector<std::vector<double>> B_prev(m_number_n); //
    m_Bk.resize(m_number_n);
    for (auto i = 0; i < m_number_n; ++i) {
      auto &Bi = m_Bk[i];
      for (std::size_t ir = 0; ir < m_rgrid_ptr->num_points; ++ir) {
        auto Br = (m_knots_t[i] <= ir && ir < m_knots_t[i + 1]) ? 1.0 : 0;
        Bi.push_back(Br);
      }
    }
    // auto B_prev = m_Bk;

    for (int k = 2; k <= m_order_k; k++) {
      auto B_prev = m_Bk;
      for (std::size_t i = 0; i < (std::size_t)m_number_n; ++i) {
        for (std::size_t ir = 0; ir < m_rgrid_ptr->num_points; ++ir) {
          auto ri = m_rgrid_ptr->r[ir];
          auto rti = m_rgrid_ptr->r[m_knots_t[i]];
          auto rtipk = m_rgrid_ptr->r[m_knots_t[i + k]];
          auto rtip1 = m_rgrid_ptr->r[m_knots_t[i + 1]];
          auto rtipkm1 = m_rgrid_ptr->r[m_knots_t[i + k - 1]];
          auto ratio_1 = (ri - rti) / (rtipkm1 - rti);
          auto ratio_2 = (rtipk - ri) / (rtipk - rtip1);
          auto Bki = ratio_1 * B_prev[i][ir];
          // std::cout << i << " ";
          // printf("%.1e, %.1e, %.1e, %.1e, %.1e | %.1e, %.1e : ", ri, rti,
          // rtipk,
          //        rtip1, rtipkm1, ratio_1, ratio_2);
          // std::cout << Bki << " ";
          if (i + 1 < (std::size_t)m_number_n)
            Bki += ratio_2 * B_prev[i + 1][ir];
          // std::cout << Bki << "\n";
          m_Bk[i][ir] = Bki;
        }
      }
    }

    // m_Bk = B_prev;
  }

  void write_splines() {
    std::ofstream of("out.txt");

    for (std::size_t ir = 0; ir < m_rgrid_ptr->num_points; ++ir) {
      of << m_rgrid_ptr->r[ir] << " ";
      for (const auto &Bi : m_Bk) {
        of << Bi[ir] << " ";
      }
      of << "\n";
    }
  }
};
