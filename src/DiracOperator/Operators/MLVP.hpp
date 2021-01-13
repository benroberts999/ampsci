#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <iostream>

// added magnetic loop function outside the name space
inline double gslfunc_ML_added(double t, void *p) {
    const auto r =  *(static_cast<double *>(p)); // get the value of r from outside
    
    // form the integrant
    double factor = 2.0/3.0*PhysConst::alpha/M_PI;
    double a = std::sqrt(1.0-0.5/(t*t));
    double b = (1.0+0.5/(t*t))*1.0/t;
    double c = std::exp(-2.0*r*t/PhysConst::alpha)*(2.0*r*t/PhysConst::alpha+1.0);
    
    return factor*a*b*c;

}

namespace DiracOperator {


class MLVP final : public TensorOperator {

public: // constructor
  MLVP(const TensorOperator *const h0, const Grid &rgrid)
      : TensorOperator(
            h0->rank(), h0->parity() == 1 ? Parity::even : Parity::odd,
            h0->getc(), MLVP_func(rgrid, h0->getv()), h0->get_d_order(),
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
  static std::vector<double> MLVP_func(const Grid &rgrid, std::vector<double> v = {}) {

    const double a0 = PhysConst::alpha;
    if (v.empty()) {
      // If v is empty, means it should be {1,1,1,1,...}
      v.resize(rgrid.num_points, 1.0);
    }
    

    // variables needed for gsl
    static constexpr double abs_err_lim = 0.0;
    static constexpr double rel_err_lim = 1.0e-10;
    static constexpr unsigned long max_num_subintvls = 1750; //?
    gsl_set_error_handler_off(); //?

    // compute the integral in parallel
    #pragma omp parallel for
    for (auto i = 0ul; i < rgrid.num_points; ++i) {
      
      // allocate the memory for required variables
      double *r_val = (double*) malloc(sizeof(double));
      double *int_result = (double*) malloc(sizeof(double));
      double *abs_err = (double*) malloc(sizeof(double));
      
      // seems to diverge if not set to zero? 
      *int_result = 0.0;
      *abs_err = 0.0;


      *r_val = rgrid.r[i]; // get the point on the radial grid
      
      // intialise gsl
      gsl_function f_gsl;
      f_gsl.function = &gslfunc_ML_added; // trying to declate gsl_function
      f_gsl.params = r_val;
      
      // integrate using gsl
      gsl_integration_workspace *gsl_int_wrk =
          gsl_integration_workspace_alloc(max_num_subintvls + 1);

      gsl_integration_qagiu(&f_gsl, 1.0, abs_err_lim, rel_err_lim,
                            max_num_subintvls, gsl_int_wrk, int_result, abs_err);
      gsl_integration_workspace_free(gsl_int_wrk);

      // multiply the operator
      v[i] *= *int_result;
    
      // free the memory
      free(int_result);
      free(abs_err);
      free(r_val);  
  }
    return v;
  }
};


} // namespace DiracOperator
