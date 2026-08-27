#pragma once
#include "Angular/Wigner369j.hpp"
#include "DiracOperator/Operators/Vee.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace Vee_pots {
DiracSpinor V_Fv(const std::vector<DiracSpinor> &core, const bool eN,
                 const DiracSpinor &Fv, const std::string type,
                 const int kappa_n, const double y, const bool contact,
                 const double mu);
DiracSpinor calc_Fdir(const std::string type, const DiracSpinor &Fa,
                      const DiracSpinor &Fv, const bool contact,
                      const double mu);
DiracSpinor calc_Fexch(const std::string type, const int kappa_n,
                       const DiracSpinor &Fa, const DiracSpinor &Fv,
                       const bool contact, const double mu);

DiracSpinor Bk_ab_v(const int k, const bool contact, const double mu,
                    const bool betaalpha, const std::string type,
                    const DiracSpinor &Fa, const DiracSpinor &Fc,
                    const DiracSpinor &Fv);

double u_anav_contact(const Wavefunction &wf, const DiracSpinor &Fn,
                      const DiracSpinor &Fa, const DiracSpinor &Fv);
double u_anva_contact(const Wavefunction &wf, const DiracSpinor &Fn,
                      const DiracSpinor &Fa, const DiracSpinor &Fv);
double R_AV_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
                 const DiracSpinor &Fc, const DiracSpinor &Fd);
DiracSpinor g0(const DiracSpinor &Fa);
DiracSpinor g5(const DiracSpinor &Fa);
DiracSpinor old_ig5(const DiracSpinor &Fa);
DiracSpinor old_i_g0_g5(const DiracSpinor &Fa);
double mod_sph_bessel_i(double n, double x);
double mod_sph_bessel_k(double n, double x);

} // namespace Vee