#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "Maths/Grid.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Wavefunction/Wavefunction.hpp"
// #include "IO/InputBlock.hpp"

namespace DiracOperator {

// XXX Check the 4*pi!???

//! Matrix element of tensor operator: J_L(qr)*C^L
/*!@details
Common matrix element for scattering problems.
Defined as matrix element of J_L(qr)*C^L.
Note: Does NOT include the Sqrt(2L+1) term (or 4pi).
JL is spherical bessel function.
Uses a look-up table for jL, faster than calculating on-the-fly.
Note: different to other operators: *rank* (rank=L) and q are variable.
Must use set_L_q (which sets L and q value) prior to use.
*/
class jL : public TensorOperator {
public:
  //! Contruction takes radial grid, a q grid, and a maximum L. Fills lookup table
  jL(const Grid &r_grid, const Grid &q_grid, std::size_t max_l)
      : TensorOperator(-1, Parity::blank),
        m_max_l(max_l),
        m_j_lq_r(m_max_l + 1, q_grid.num_points()),
        m_r_grid(&r_grid),
        m_q_grid(&q_grid) {
    fill_table();
  }

protected:
  std::size_t m_max_l;
  LinAlg::Matrix<std::vector<double>> m_j_lq_r;
  const Grid *m_r_grid;
  const Grid *m_q_grid;

private:
  // fills lookup J_l_q_r table
  void fill_table() {
    std::cout << "Filling jL lookup table:\n";
    for (std::size_t l = 0; l <= m_max_l; ++l) {
#pragma omp parallel for
      for (std::size_t iq = 0; iq < m_q_grid->num_points(); ++iq) {
        m_j_lq_r.at(l, iq) = SphericalBessel::fillBesselVec_kr(
            int(l), m_q_grid->r(iq), m_r_grid->r());
      }
    }
  }

public:
  //! Current value of L (should = rank)
  std::size_t L() const { return std::size_t(m_rank); }
  //! Maximum L value in table.
  std::size_t max_L() const { return m_max_l; }

  //! Sets the current L and q values for use. Note: NOT thread safe!
  virtual void set_L_q(std::size_t L, double q) {
    assert(L <= m_max_l && "L must be <= max L");
    m_rank = int(L);
    // opposite parity for 'pseudo' cases?
    m_parity = Angular::evenQ(m_rank) ? Parity::even : Parity::odd;
    const auto iq = m_q_grid->getIndex(q);
    m_vec = m_j_lq_r.at(L, iq);
  }

public:
  virtual double angularF(const int ka, const int kb) const override {
    // 4 pi?
    return Angular::Ck_kk(m_rank, ka, kb);
  }
  virtual double angularCff(int, int) const override { return 1.0; }
  virtual double angularCgg(int, int) const override { return 1.0; }
  virtual double angularCfg(int, int) const override { return 0.0; }
  virtual double angularCgf(int, int) const override { return 0.0; }

  virtual std::string name() const override {
    return std::string("j") + std::to_string(m_rank);
  }
  std::string units() const override final { return "au"; }
};

//------------------------------------------------------------------------------
//! Matrix element of tensor operator: gamma^0 J_L(qr) C^L
class g0jL : public jL {
public:
  g0jL(const Grid &r_grid, const Grid &q_grid, std::size_t max_l)
      : jL(r_grid, q_grid, max_l) {}

public:
  double angularCff(int, int) const override final { return 1.0; }
  double angularCgg(int, int) const override final { return -1.0; }
  double angularCfg(int, int) const override final { return 0.0; }
  double angularCgf(int, int) const override final { return 0.0; }

  std::string name() const override final {
    return std::string("g0j") + std::to_string(m_rank);
  }
};

//------------------------------------------------------------------------------
//! Matrix element of tensor operator: i gamma^5 J_L(qr) C^L
class ig5jL : public jL {
public:
  ig5jL(const Grid &r_grid, const Grid &q_grid, std::size_t max_l)
      : jL(r_grid, q_grid, max_l) {}

public:
  void set_L_q(std::size_t L, double q) override final {
    assert(L <= m_max_l && "L must be <= max L");
    m_rank = int(L);
    // Note: opposite parity for 'pseudo' cases
    m_parity = Angular::evenQ(m_rank) ? Parity::odd : Parity::even;
    const auto iq = m_q_grid->getIndex(q);
    m_vec = m_j_lq_r.at(L, iq);
  }

  double angularF(const int ka, const int kb) const override final {
    // 4 pi?
    return -1.0 * Angular::Ck_kk(m_rank, ka, -kb);
  }

  double angularCff(int, int) const override final { return 0.0; }
  double angularCgg(int, int) const override final { return 0.0; }
  double angularCfg(int, int) const override final { return 1.0; }
  double angularCgf(int, int) const override final { return -1.0; }

  std::string name() const override final {
    return std::string("ig5j") + std::to_string(m_rank);
  }
};

//------------------------------------------------------------------------------
//! Matrix element of tensor operator: i gamma^0gamma^5 J_L(qr) C^L
class ig0g5jL : public jL {
public:
  ig0g5jL(const Grid &r_grid, const Grid &q_grid, std::size_t max_l)
      : jL(r_grid, q_grid, max_l) {}

public:
  virtual void set_L_q(std::size_t L, double q) {
    assert(L <= m_max_l && "L must be <= max L");
    m_rank = int(L);
    // Note: opposite parity for 'pseudo' cases
    m_parity = Angular::evenQ(m_rank) ? Parity::odd : Parity::even;
    const auto iq = m_q_grid->getIndex(q);
    m_vec = m_j_lq_r.at(L, iq);
  }

  virtual double angularF(const int ka, const int kb) const override {
    // 4 pi?
    return -1.0 * Angular::Ck_kk(m_rank, ka, -kb);
  }

  virtual double angularCff(int, int) const override { return 0.0; }
  virtual double angularCgg(int, int) const override { return 0.0; }
  virtual double angularCfg(int, int) const override { return 1.0; }
  virtual double angularCgf(int, int) const override { return 1.0; }

  std::string name() const override final {
    return std::string("ig5j") + std::to_string(m_rank);
  }
};

} // namespace DiracOperator
