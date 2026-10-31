#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
//#include "Potentials/NuclearPotentials.hpp" //not necessary bc
#include "jls.hpp"

namespace DiracOperator {

//! anomalous moment in the electric field operator: <a||h2||b>
/*! @details
\f[ <a||h1||b> = -(i*a_e*alpha/2) R \f]
\f[R = \int ( 2<k_a||S||k_b> f_af_b - 2<-k_a||S||-k_b> g_ag_b ) \, dr\f]
Units of |mu_B|.
*/

class h2 final : public TensorOperator {
public:
  h2(const Wavefunction &wf)
    : TensorOperator(1, Parity::odd, +1.0, {}, Realness::imaginary, true),
      m_r(wf.grid().r()) {}
  h2 &operator=(const h2 &) = delete;
  h2(const h2 &) = default;
  ~h2() = default;
  std::string name() const override final { return std::string("h2"); }
  std::string units() const override final { return std::string("?"); }

  double angularF(const int ka, const int kb) const override final {
    if (std::abs(ka) == std::abs(kb)) {
      double twoj_plus_1 =
        Angular::jindex_to_twoj(Angular::kappa_to_jindex(ka));
      return (-1.0) * std::sqrt(twoj_plus_1);
    } else {
      return 0.0;
    }
  }
  double angularCff(int, int) const override final { return 0.0; }
  double angularCgg(int, int) const override final { return 0.0; }
  double angularCfg(int, int) const override final { return 1.0; }
  double angularCgf(int, int) const override final { return 1.0; }

  static std::vector<double> nuclear_density(const Grid &gr,
                                             const Nuclear::Nucleus &nucleus) {

    const double r_rms = nucleus.r_rms() / PhysConst::aB_fm; // convert to au
    const auto nucleus_type =
      (r_rms > 0.0 || nucleus.type() == Nuclear::ChargeDistro::custom) ?
        nucleus.type() :
        Nuclear::ChargeDistro::point;

    switch (nucleus_type) {

    case Nuclear::ChargeDistro::Fermi: {
      return Nuclear::fermiNuclearDensity_tcN(nucleus.t(), nucleus.c(), 1.0,
                                              gr);
    }
    //spherical, point, Gaussian
    case Nuclear::ChargeDistro::spherical: {
      const double rN = sqrt(5.0 / 3.0) * r_rms;
      const double rho_inside = 3.0 / (4.0 * M_PI * pow(rN, 3.0));

      std::vector<double> rho(gr.num_points(), 0.0);
      for (std::size_t i = 0; i < gr.num_points(); ++i) {
        if (gr.r(i) < rN) {
          rho[i] = rho_inside;
        }
      }
      return rho;
    }
    case Nuclear::ChargeDistro::Gaussian: {
      std::vector<double> rho(gr.num_points(), 0.0);

      const double factor = std::pow(3.0 / (2.0 * M_PI * r_rms * r_rms), 1.5);

      for (std::size_t i = 0; i < gr.num_points(); ++i) {
        const double r = gr.r(i);
        rho[i] = factor * std::exp((-3.0) * r * r / (2.0 * r_rms * r_rms));
      }
      return rho;
    }
    case Nuclear::ChargeDistro::point: {
      std::cout << "\nCannot use this operator for a point-like nucleus.\n";
      std::abort();
    }
    case Nuclear::ChargeDistro::custom: {
      std::cout
        << "\nCannot use this operator for a custom charge distribution.\n";
      std::abort();
    }
    case Nuclear::ChargeDistro::Error: {
      std::cout << "\nCannot use this operator for case "
                   "Nuclear::ChargeDistro::Error.\n";
      std::abort();
    }
    }
  }

  static std::vector<double> form_H_AMM(const Wavefunction &wf) {

    // now calculate derivate analytically?
    //need to get Z, \rho_nuc(r), and core f(r)/g(r)/j
    const Grid &gr = wf.grid();
    const auto &nucleus = wf.nucleus();
    const int Z = nucleus.z();
    const std::vector<double> rho_nuc = nuclear_density(gr, nucleus);

    std::vector<double> Z_r(gr.num_points(), 0.0); // enclosed charge
    const std::vector<double> fac_Z(gr.num_points(), 4.0 * M_PI * Z);
    std::vector<double> r_squared(gr.num_points(), 0.0);

    for (std::size_t i = 0; i < gr.num_points(); ++i) {
      r_squared[i] = gr.r(i) * gr.r(i);
    }

    NumCalc::additivePIntegral<NumCalc::zero_to_r>(Z_r, fac_Z, r_squared,
                                                   rho_nuc, gr);

    std::vector<double> Ncore_r(gr.num_points(), 0.0);
    for (const auto &Fc : wf.core()) {

      const std::vector<double> fac_Ncore(gr.num_points(), Fc.twojp1());

      NumCalc::additivePIntegral<NumCalc::zero_to_r>(Ncore_r, fac_Ncore, Fc.f(),
                                                     Fc.f(), gr);
      NumCalc::additivePIntegral<NumCalc::zero_to_r>(Ncore_r, fac_Ncore, Fc.g(),
                                                     Fc.g(), gr);
    } // end for loop over core states

    //double a_e = PhysConst::alpha / (2.0 * M_PI);
    //double fac = a_e * PhysConst::alpha / 2.0;
    double global_fac = PhysConst::alpha2 / (4.0 * M_PI);

    std::vector<double> H_AMM(gr.num_points(), 0.0);

    for (std::size_t i = 0; i < gr.num_points(); ++i) {
      H_AMM[i] = global_fac * (Z_r[i] - Ncore_r[i]) / (gr.r(i) * gr.r(i));
    }

    return H_AMM;
  }

  static std::vector<double> form_DeltaH(const Wavefunction &wf) {

    const Grid &gr = wf.grid();
    std::vector<double> Ncore_r(gr.num_points(), 0.0);
    for (const auto &Fc : wf.core()) {

      const std::vector<double> fac_Ncore(gr.num_points(), Fc.twojp1());

      NumCalc::additivePIntegral<NumCalc::zero_to_r>(Ncore_r, fac_Ncore, Fc.f(),
                                                     Fc.f(), gr);
      NumCalc::additivePIntegral<NumCalc::zero_to_r>(Ncore_r, fac_Ncore, Fc.g(),
                                                     Fc.g(), gr);
    } // end for loop over core states

    //double a_e = PhysConst::alpha / (2.0 * M_PI);
    //double fac = a_e * PhysConst::alpha / 2.0;
    double global_fac = PhysConst::alpha2 / (4.0 * M_PI);
    std::vector<double> DeltaH(gr.num_points(), 0.0);
    for (std::size_t i = 0; i < gr.num_points(); ++i) {
      DeltaH[i] = -global_fac * Ncore_r[i] / (gr.r(i) * gr.r(i));
    }
    return DeltaH;
  }

  std::unique_ptr<TensorOperator> clone() const override final {
    return std::make_unique<h2>(*this);
  }

  static std::unique_ptr<TensorOperator> generate(const IO::InputBlock &input,
                                                  const Wavefunction &wf) {
    input.check({{"", "No input options"}});
    if (input.has_option("help"))
      return nullptr;
    return std::make_unique<h2>(wf);
  }

private:
  const std::vector<double> m_r; // store radial vector (copy)
};

} // namespace DiracOperator