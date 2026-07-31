#pragma once
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace Module {
void BSM_Vee(const IO::InputBlock &input, const Wavefunction &wf);
void test_CI(const IO::InputBlock &input, const Wavefunction &wf);
void find_fierz(const IO::InputBlock &input, const Wavefunction &wf);
double CI_edm(const Wavefunction &wf, const std::string basis_string,
              const double mu, const bool contact, const int J_V,
              const int pi_V, const int i_V, const bool verbose = true);
double calc_ME(const std::string op, const std::string type, const bool contact,
               const double mu, const bool tdhf, const double omega,
               const Wavefunction &wf, const DiracSpinor &Fw,
               const DiracSpinor &Fv);
double calc_Vwv(const std::string type, const bool contact, const double mu,
                const bool tdhf, const double omega, const Wavefunction &wf,
                const DiracSpinor &Fw, const DiracSpinor &Fv);
double calc_Dwv(const std::string type, const bool contact, const double mu,
                const bool tdhf, const double omega, const Wavefunction &wf,
                const DiracSpinor &Fw, const DiracSpinor &Fv);
// double calc_Dwv(const std::string type, const bool contact, const double mu,
//                 const bool tdhf, const Wavefunction &wf, const DiracSpinor &Fv);
// double calc_Dwv(const std::string type, const bool contact, const double mu,
//                 const bool tdhf, const Wavefunction &wf);
void ee_isotope_shift(const std::string int_type, const IO::InputBlock &input,
                      const Wavefunction &wf);
double dE(const double mu, const std::string int_type,
          const std::vector<DiracSpinor> &core, const DiracSpinor Fv);

void matrix_elements(const Wavefunction &wf, const DiracSpinor &Fv,
                     const DiracSpinor &Fw, const std::string op,
                     const std::string type, const bool contact,
                     const bool tdhf = true, const double omega = 0.0,
                     const double min_mu = 1e-04, const double max_mu = 1e4,
                     const int N_mu = 100);
// ############## OLD FUNCTIONS #################

double d_ab(const Grid &gr, const DiracSpinor &Fa, const DiracSpinor &Fb);
double V_nv_direct(const bool contact, const std::vector<DiracSpinor> core,
                   const DiracSpinor &Fv, const DiracSpinor &Fn, const double y,
                   const double mu, const bool g0_both = true);
double R_abcd_contact(const double mu, const DiracSpinor &Fa,
                      const DiracSpinor &Fb, const DiracSpinor &Fc,
                      const DiracSpinor &Fd);
double Rk_abcd(const double k, const double mu, const DiracSpinor &Fa,
               const DiracSpinor &Fb, const DiracSpinor &Fc,
               const DiracSpinor &Fd, const std::string int_type = "sp",
               const bool g0_both = true);
double Rk_abcd_massless(const double k, const DiracSpinor &Fa,
                        const DiracSpinor &Fb, const DiracSpinor &Fc,
                        const DiracSpinor &Fd);
std::vector<double> Bk_ab(const double k, const double mu,
                          const DiracSpinor &Fa, const DiracSpinor &Fb);

void sps_testing(const Wavefunction &wf, const bool contact);
} // namespace Module
