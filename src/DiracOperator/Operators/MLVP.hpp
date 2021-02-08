#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <iostream>

// struct to pass the parameters to the GSL function
struct MLVP_params {
    // simple struct that only stores current point on the radial grid and the nuclear radius
    double r, rN;
};

// added magnetic loop function outside the name space
// Note that this is ball nuclear magnetisation MLVP function which multiplies 1/r**2
// hence has to be executed with "pointlike" F(r)
inline double gslfunc_ML_added(double t, void *p) {

    // obtain the radial grid parameters and convert to relativistic units
    auto [r, r_N] = *(static_cast<MLVP_params *>(p));
    r = r/PhysConst::alpha;
    r_N = r_N/PhysConst::alpha;

    // form the integrand (all calculations are preformed in relativistic units)
    double prefactor = 2./3.*PhysConst::alpha/M_PI*3./(16.*r_N*r_N*r_N);
    double t_part = std::sqrt(t*t-1.)/(t*t*t)*(1.+1./(2.*t*t));
    double plus_r = (4.*r*r_N+2./t*(r_N+r)+1./(t*t))*std::exp(-2.*t*(r+r_N)); 
    double minus_r = (4.*r*r_N-2./t*std::fabs(r_N-r)-1./(t*t))*std::exp(-2.*t*std::fabs(r_N-r));  
 
    return prefactor*t_part*(plus_r+minus_r);

}

namespace DiracOperator {


class MLVP final : public TensorOperator {

public: // constructor
  MLVP(const TensorOperator *const h0, const Grid &rgrid, double rN)
      : TensorOperator(
            h0->rank(), h0->parity() == 1 ? Parity::even : Parity::odd,
            h0->getc(), MLVP_func(rgrid, rN, h0->getv()), h0->get_d_order(),
            h0->imaginaryQ() ? Realness::imaginary : Realness::real,
            h0->freqDependantQ),
        m_h0(h0)
  // , m_a(a), m_b(b)
  {}

  std::string name() const override final {
    return m_h0->name() + "_MLVP";
  }
  std::string units() const override final { return m_h0->units(); }

  double angularF(const int ka, const int kb) const override final {
    return m_h0->angularF(ka, kb);
  }

  double angularCff(int ka, int kb) const override final {
    return m_h0->angularCff(ka, kb);
  }
  double angularCgg(int ka, int kb) const override final {
    return m_h0->angularCgg(ka, kb);
  }
  double angularCfg(int ka, int kb) const override final {
    return m_h0->angularCfg(ka, kb);
  }
  double angularCgf(int ka, int kb) const override final {
    return m_h0->angularCgf(ka, kb);
  }

  // Have m_h0 pointer, so delete copy/asign constructors
  MLVP(const DiracOperator::MLVP &) = delete;
  MLVP &operator=(const DiracOperator::MLVP &) = delete;

private:
  const TensorOperator *const m_h0;

public:
  static std::vector<double> MLVP_func(const Grid &rgrid, double rN, std::vector<double> v = {}) {

    const double a0 = PhysConst::alpha;
    if (v.empty()) {
      // If v is empty, means it should be {1,1,1,1,...}
      v.resize(rgrid.num_points, 1.0);
    }
    

    // variables needed for gsl
    static constexpr double abs_err_lim = 0.0;
    static constexpr double rel_err_lim = 1.0e-8;
    static constexpr unsigned long max_num_subintvls = 750; //?
    gsl_set_error_handler_off(); //?

    // allocate the memory for required variables
    double int_result = 0;
    double abs_err = 0;
    
    // compute the integral at each radial grid point
    for (auto i = 0ul; i < rgrid.num_points; ++i) {
      
      
      // seems to diverge if not set to zero? 
      int_result = 0.0;
      abs_err = 0.0;


      MLVP_params params = {rgrid.r[i], rN}; // get the point on the radial grid
      
      // intialise gsl
      gsl_function f_gsl;
      f_gsl.function = &gslfunc_ML_added; // trying to declate gsl_function
      f_gsl.params = &params;
      
      // integrate using gsl
      gsl_integration_workspace *gsl_int_wrk =
          gsl_integration_workspace_alloc(max_num_subintvls + 1);

      gsl_integration_qagiu(&f_gsl, 1.0, abs_err_lim, rel_err_lim,
                            max_num_subintvls, gsl_int_wrk, &int_result, &abs_err);
      gsl_integration_workspace_free(gsl_int_wrk);

      // multiply the operator
      v[i] *= int_result;
    
    }
    return v;
  }
};


} // namespace DiracOperator
