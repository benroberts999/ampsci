#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "Angular/Wigner369j.hpp"

namespace DiracOperator {

//==============================================================================
//! @brief Nuclear-spin independent PNC operator (Qw)
/*! @details
  \f[
    h_{PNCnsi} = -\frac{G_F \, Q_W}{2\sqrt{2}} \, \rho(r) \, \gamma^5.
  \f]
  \f[
    h_{PNCnsi} = \frac{G_F \, Q_W}{2\sqrt{2}} \, \rho(r) \, \gamma_5.
  \f]

  - Ouput is in units of (Qw * 1.e-11.) by default. To get (Qw/-N), multiply
    by (-N) [can go into optional 'factor']

  Scalar, ME = Radial integral (F'|h|F)
  \f[
  (F'|h|F) = -\frac{G_F \, Q_W}{2\sqrt{2}}\int (f'g - g'f) rho(r) dr
  \f]

  Generates rho(r) according to fermi distribution, given c and t [c and t in
  FERMI / femptometers].
*/
class PNCnsi final : public ScalarOperator {
public:
  PNCnsi(double c, double t, const Grid &rgrid, double factor = 1.0,
         const std::string &in_units = "iQw*e-11")
    : ScalarOperator(Parity::odd, factor * PhysConst::GFe11 / std::sqrt(8.0),
                     Nuclear::fermiNuclearDensity_tcN(t, c, 1.0, rgrid),
                     {0, -1, +1, 0}, Realness::imaginary),
      m_unit(in_units) {}

  std::string name() const override final { return "pnc-nsi"; }
  std::string units() const override final { return m_unit; }

  std::unique_ptr<TensorOperator> clone() const override final {
    return std::make_unique<PNCnsi>(*this);
  }

  static std::unique_ptr<TensorOperator> generate(const IO::InputBlock &input,
                                                  const Wavefunction &wf) {
    input.check(
      {{{"c",
         "Half-density radius for Fermi rho(r). [defaut: from wavefunction]"},
        {"t", "skin thickness [2.3]"},
        {"N", "Neutron number, for units [default: from wavefunction]"},
        {"print", "Write details to screen [true]"}}});
    if (input.has_option("help"))
      return nullptr;
    const auto r_rms = wf.get_rrms();
    const auto c = input.get("c", Nuclear::c_hdr_formula_rrms_t(r_rms));
    const auto t = input.get("t", Nuclear::default_t);
    const auto N = input.get("N", wf.Anuc() - wf.Znuc());
    std::string units =
      N == 1 ? "i(-Qw)e-11" : "i(-Qw/" + std::to_string(N) + ")e-11";
    if (input.get("print", true))
      std::cout << "pnc: with c=" << c << ", t=" << t << " [" << units << "]\n";
    return std::make_unique<PNCnsi>(c, t, wf.grid(), -1.0 * N, "units");
  }

private:
  const std::string m_unit{"iQw*e-11"};
};


// anapole moment/spin dependent pnc
class PNCsd : public TensorOperator {
  public:
    PNCsd(double c, double t, double I, const Grid &rgrid)
      : TensorOperator(1,  Parity::odd, I * PhysConst::GFe11 / std::sqrt(2.0),
        Nuclear::fermiNuclearDensity_tcN(t, c, 1.0, rgrid), Realness::imaginary)
       {}
  
    std::unique_ptr<TensorOperator> clone() const override {
      return std::make_unique<PNCsd>(*this);
    }
  
    std::string name() const override {
      return std::string("PNCsd");
    }
    std::string units() const override {
      return std::string("iκ e-11");
    }
  
    double angularF(const int, const int) const override final {
      return 1.0;
    }

    double angularCff(int kappa_a, int kappa_b) const override final {
      (void)kappa_a, (void)kappa_b;
      return 0.0;
    }
  
    double angularCgg(int, int) const override final { return 0.0; }
    
    double angularCfg(int ka, int kb) const override final { return 2.0*Angular::S_kk(ka,-kb); }
    
    double angularCgf(int ka, int kb) const override final { return -2.0*Angular::S_kk(-ka,kb); }
  
    static std::unique_ptr<TensorOperator> generate(const IO::InputBlock &input,
      const Wavefunction &wf) {
        input.check(
          {{{"c",
              "Half-density radius for Fermi rho(r). [defaut: from wavefunction]"},
            {"t", "skin thickness [2.3]"},
            {"print", "Write details to screen [true]"}}});
        if (input.has_option("help"))
          return nullptr;
        const auto r_rms = wf.get_rrms();
        const auto c = input.get("c", Nuclear::c_hdr_formula_rrms_t(r_rms));
        const auto t = input.get("t", Nuclear::default_t);
        if (input.get("print", true))
          std::cout << "spin-dep pnc: with c=" << c << ", t=" << t << "\n";
        return std::make_unique<PNCsd>(c, t, wf.grid());
      }

  };
} // namespace DiracOperator
