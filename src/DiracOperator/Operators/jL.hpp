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
Must use set_L_q (which sets L and q value) prior to use. Note: this is not thread-safe, so cannot parallelise over L or q.
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
  //! Constructing from existing operator: copies JL table - faster. Can copy from a jL of different type (g0,g5 etc)
  jL(const jL &other)
      : TensorOperator(-1, Parity::blank),
        m_max_l(other.m_max_l),
        m_j_lq_r(other.m_j_lq_r),
        m_r_grid(other.m_r_grid),
        m_q_grid(other.m_q_grid) {}

  // jL(const jL &) = default;
  jL &operator=(const jL &) = delete;

protected:
  std::size_t m_max_l;
  LinAlg::Matrix<std::vector<double>> m_j_lq_r;
  const Grid *m_r_grid;
  const Grid *m_q_grid;

private:
  // fills lookup J_l_q_r table
  void fill_table() {
    std::cout << "Filling jL lookup table: ";
    for (std::size_t l = 0; l <= m_max_l; ++l) {
#pragma omp parallel for
      for (std::size_t iq = 0; iq < m_q_grid->num_points(); ++iq) {
        m_j_lq_r.at(l, iq) = SphericalBessel::fillBesselVec_kr(
            int(l), m_q_grid->r(iq), m_r_grid->r());
      }
    }
    std::cout << "done.\n";
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
    m_parity = Angular::evenQ(m_rank) ? Parity::even : Parity::odd;
    const auto iq = m_q_grid->getIndex(q);
    m_vec = m_j_lq_r.at(L, iq);
  }

  //! Directly calculate reduced matrix element without needing to call set_L_q - this *is* thread safe
  double rme(const DiracSpinor &a, const DiracSpinor &b, std::size_t L,
             double q) const {
    const auto iq = m_q_grid->getIndex(q);
    const auto &jlqr = m_j_lq_r.at(L, iq);
    const auto &gr = *m_r_grid;
    const auto cff = angularCff(a.kappa(), b.kappa());
    const auto cgg = angularCgg(a.kappa(), b.kappa());
    const auto cfg = angularCfg(a.kappa(), b.kappa());
    const auto cgf = angularCgf(a.kappa(), b.kappa());
    auto max = std::min(a.max_pt(), b.max_pt());
    const auto Rf = cff == 0.0 ? 0.0 :
                                 NumCalc::integrate(1.0, 0, max, a.f(), b.f(),
                                                    jlqr, gr.drdu());
    const auto Rg = cgg == 0.0 ? 0.0 :
                                 NumCalc::integrate(1.0, 0, max, a.g(), b.g(),
                                                    jlqr, gr.drdu());
    const auto Rfg = cfg == 0.0 ? 0.0 :
                                  NumCalc::integrate(1.0, 0, max, a.f(), b.g(),
                                                     jlqr, gr.drdu());
    const auto Rgf = cgf == 0.0 ? 0.0 :
                                  NumCalc::integrate(1.0, 0, max, a.g(), b.f(),
                                                     jlqr, gr.drdu());

    const auto kb = cfg == 0.0 ? b.kappa() : -b.kappa();
    const auto s = cfg == 0.0 ? 1.0 : -1.0;
    return s * Angular::Ck_kk(int(L), a.kappa(), kb) *
           (cff * Rf + cgg * Rg + cfg * Rfg + cgf * Rgf) * gr.du();
  }
  //! Checks if specific ME is zero (when not useing set_L_q)
  bool is_zero(const DiracSpinor &a, const DiracSpinor &b,
               std::size_t L) const {
    const auto kb =
        angularCfg(a.kappa(), b.kappa()) == 0.0 ? b.kappa() : -b.kappa();
    return !Angular::Ck_kk_SR(int(L), a.kappa(), kb);
  }

public:
  virtual double angularF(const int ka, const int kb) const override {
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
  g0jL(const jL &other) : jL(other) {}

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
//! Matrix element of tensor operator: i gamma^5 J_L(qr) C^L. nb: i makes ME real
class ig5jL : public jL {
public:
  ig5jL(const Grid &r_grid, const Grid &q_grid, std::size_t max_l)
      : jL(r_grid, q_grid, max_l) {}
  ig5jL(const jL &other) : jL(other) {}

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
//! Matrix element of tensor operator: i gamma^0gamma^5 J_L(qr) C^L. nb: i makes ME real
class ig0g5jL : public jL {
public:
  ig0g5jL(const Grid &r_grid, const Grid &q_grid, std::size_t max_l)
      : jL(r_grid, q_grid, max_l) {}
  ig0g5jL(const jL &other) : jL(other) {}

public:
  void set_L_q(std::size_t L, double q) override final {
    assert(L <= m_max_l && "L must be <= max L");
    m_rank = int(L);
    // Note: opposite parity for 'pseudo' cases
    m_parity = Angular::evenQ(m_rank) ? Parity::odd : Parity::even;
    const auto iq = m_q_grid->getIndex(q);
    m_vec = m_j_lq_r.at(L, iq);
  }

  virtual double angularF(const int ka, const int kb) const override {
    return -1.0 * Angular::Ck_kk(m_rank, ka, -kb);
  }

  virtual double angularCff(int, int) const override { return 0.0; }
  virtual double angularCgg(int, int) const override { return 0.0; }
  virtual double angularCfg(int, int) const override { return 1.0; }
  virtual double angularCgf(int, int) const override { return 1.0; }

  std::string name() const override final {
    return std::string("ig0g5j") + std::to_string(m_rank);
  }
};

} // namespace DiracOperator
