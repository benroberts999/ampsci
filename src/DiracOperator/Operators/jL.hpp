#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "Maths/Grid.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace DiracOperator {

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
  jL(const Grid &r_grid, const Grid &q_grid, std::size_t max_l,
     bool subtract_one = false)
      : TensorOperator(-1, Parity::blank),
        m_max_l(max_l),
        m_j_lq_r(m_max_l + 1, q_grid.num_points()),
        m_r_grid(&r_grid),
        m_q_grid(&q_grid),
        m_subtract_one(subtract_one) {
    fill_table();
  }
  //! Constructing from existing operator: copies JL table - faster. Can copy from a jL of different type (g0,g5 etc)
  jL(const jL &other)
      : TensorOperator(-1, Parity::blank),
        m_max_l(other.m_max_l),
        m_j_lq_r(other.m_j_lq_r),
        m_r_grid(other.m_r_grid),
        m_q_grid(other.m_q_grid),
        m_subtract_one(other.m_subtract_one) {}

  // jL(const jL &) = default;
  jL &operator=(const jL &) = delete;

protected:
  std::size_t m_max_l;
  LinAlg::Matrix<std::vector<double>> m_j_lq_r;
  const Grid *m_r_grid;
  const Grid *m_q_grid;
  bool m_subtract_one;

private:
  // fills lookup J_l_q_r table
  void fill_table() {
    std::cout << "Filling jL lookup table: " << std::flush;
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

  const auto &q_grid() const { return *m_q_grid; }
  const auto &r_grid() const { return *m_r_grid; }

  //! Sets the current L and q values for use. Note: NOT thread safe!
  virtual void set_L_q(std::size_t L, double q) {
    assert(L <= m_max_l && "L must be <= max L");
    m_rank = int(L);
    m_parity = Angular::evenQ(m_rank) ? Parity::even : Parity::odd;
    const auto iq = m_q_grid->getIndex(q);
    m_vec = m_j_lq_r.at(L, iq);
  }

  // This is _not_ thread safe
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override {
    // Fa * radial_rhs(Fa.kappa(),Fb) = h.radialIntegral(Fa, Fb)

    // const auto &gr = Fb.grid();
    DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
    dF.min_pt() = Fb.min_pt();
    dF.max_pt() = Fb.max_pt();

    if (isZero(kappa_a, Fb.kappa())) {
      dF.min_pt() = Fb.min_pt();
      dF.max_pt() = Fb.min_pt();
      return dF;
    }

    const auto &df = Fb.f();
    const auto &dg = Fb.g();

    const auto cff = angularCff(kappa_a, Fb.kappa());
    const auto cgg = angularCgg(kappa_a, Fb.kappa());
    const auto cfg = angularCfg(kappa_a, Fb.kappa());
    const auto cgf = angularCgf(kappa_a, Fb.kappa());

    if (m_subtract_one && m_rank == 0 && Fb.kappa() == kappa_a) {
      // jL -> (jL-1)
      // Vector
      //   (ff + gg)*jL - (ff+gg)
      // = (ff + gg)*(jL-1)

      // Scalar:
      // (ff - gg)*jL - (ff+gg)
      // = ff*(jL-1) - gg*(jL+1)
      if (cgg > 0.0) {
        //vector
        for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
          dF.f(i) = m_constant * (m_vec[i] - 1.0) * cff * df[i];
          dF.g(i) = m_constant * (m_vec[i] - 1.0) * cgg * dg[i];
        }
      } else if (cgg < 0.0) {
        //scalar
        for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
          dF.f(i) = m_constant * (m_vec[i] - 1.0) * cff * df[i];
          dF.g(i) = m_constant * (m_vec[i] + 1.0) * cgg * dg[i];
        }
      } else {
        assert(false && "Error 111 - subtract_one for pseudo case?");
      }
    } else {
      for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
        dF.f(i) = m_constant * m_vec[i] * (cff * df[i] + cfg * dg[i]);
        dF.g(i) = m_constant * m_vec[i] * (cgf * df[i] + cgg * dg[i]);
      }
    }

    return dF;
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

    double Rf{0.0}, Rg{0.0};
    if (m_subtract_one && a.kappa() == b.kappa() && L == 0) {
      auto jl_neg1 = jlqr;
      for (auto &el : jl_neg1) {
        el -= 1.0;
      }
      // jL -> (jL-1)
      // Vector
      //   (ff + gg)*jL - (ff+gg)
      // = (ff + gg)*(jL-1)

      // Scalar:
      // (ff - gg)*jL - (ff+gg)
      // (ff - gg)*jL - (ff-gg)
      // = ff*(jL-1) - gg*(jL+1)

      Rf = NumCalc::integrate(1.0, 0, max, a.f(), b.f(), jl_neg1, gr.drdu());

      if (cgg > 0.01) {
        //vector
        Rg = NumCalc::integrate(1.0, 0, max, a.g(), b.g(), jl_neg1, gr.drdu());
      } else if (cgg < -0.01) {
        //scalar (-ve sign is inside cgg)
        auto jl_p1 = jlqr;
        for (auto &el : jl_p1) {
          el += 1.0;
        }
        Rg = NumCalc::integrate(1.0, 0, max, a.g(), b.g(), jl_p1, gr.drdu());
      } else {
        assert(false && "Error: subtract 1 in pseudo-case?");
      }
    } else {
      Rf = cff == 0.0 ?
               0.0 :
               NumCalc::integrate(1.0, 0, max, a.f(), b.f(), jlqr, gr.drdu());
      Rg = cgg == 0.0 ?
               0.0 :
               NumCalc::integrate(1.0, 0, max, a.g(), b.g(), jlqr, gr.drdu());
    }
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

  virtual std::string name() const override { return std::string("jL"); }
  std::string units() const override final { return "au"; }
};

//------------------------------------------------------------------------------
//! Matrix element of tensor operator: gamma^0 J_L(qr) C^L
class g0jL : public jL {
public:
  g0jL(const Grid &r_grid, const Grid &q_grid, std::size_t max_l,
       bool subtract_one = false)
      : jL(r_grid, q_grid, max_l, subtract_one) {}
  g0jL(const jL &other) : jL(other) {}

public:
  double angularCff(int, int) const override final { return 1.0; }
  double angularCgg(int, int) const override final { return -1.0; }
  double angularCfg(int, int) const override final { return 0.0; }
  double angularCgf(int, int) const override final { return 0.0; }

  std::string name() const override final { return std::string("g0jL"); }
};

//------------------------------------------------------------------------------
//! Matrix element of tensor operator: i gamma^5 J_L(qr) C^L. nb: i makes ME real
class ig5jL : public jL {
public:
  ig5jL(const Grid &r_grid, const Grid &q_grid, std::size_t max_l)
      : jL(r_grid, q_grid, max_l, false) {}
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

  std::string name() const override final { return std::string("ig5jL"); }
};

//------------------------------------------------------------------------------
//! Matrix element of tensor operator: i gamma^0gamma^5 J_L(qr) C^L. nb: i makes ME real
class ig0g5jL : public jL {
public:
  ig0g5jL(const Grid &r_grid, const Grid &q_grid, std::size_t max_l)
      : jL(r_grid, q_grid, max_l, false) {}
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

  std::string name() const override final { return std::string("ig0g5jL"); }
};

} // namespace DiracOperator
