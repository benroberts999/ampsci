#pragma once
#include "Angular/Angular.hpp"
#include <vector>
class DiracSpinor;
class DiracOperator;
class Wavefunction;

// Allow it to do RPA daigram-style as well!
//(needs a basis)

enum class dPsiType { X, Y };

class ExternalField {
public:
  ExternalField(const DiracOperator *const h,
                const std::vector<DiracSpinor> &core,
                const std::vector<double> &vl, const double alpha);

private:
  // dPhi = X exp(-iwt) + Y exp(+iwt)
  // (H - e - w)X = -(h + dV - de)Phi
  // (H - e + w)Y = -(h* + dV* - de)Phi
  // X_c = sum_x X_x,
  // j(x)=j(c)-k,...,j(c)+k.  And: pi(x) = pi(c)*pi(h)
  std::vector<std::vector<DiracSpinor>> m_X; // X[core_state][kappa_x]
  std::vector<std::vector<DiracSpinor>> m_Y;

  const DiracOperator *const m_h; //??
  const std::vector<DiracSpinor> *const p_core;
  const std::vector<double> m_vl;
  const double m_alpha;
  // const double m_omega;
  // const bool static_fieldQ;
  const int m_rank;
  const int m_pi;
  const bool m_imag;

  Angular::SixJ m_6j; // = Angular::SixJ(m_rank, 7); // XXX temp!

public:
  std::vector<DiracSpinor> &get_dPsis(const DiracSpinor &phic, dPsiType XorY);
  const DiracSpinor &get_dPsi_x(const DiracSpinor &phic, dPsiType XorY,
                                const int kappa_x);

  void solve_TDHFcore(const double omega, int max_its = 100);
  void solve_TDHFcore_matrix(const Wavefunction &wf, const double omega,
                             const int max_its = 30);

  // does it matter if a or b is in the core?
  double dV_ab(const DiracSpinor &phia, const DiracSpinor &phib,
               bool conj = false);
  // double dV_ab_Y(const DiracSpinor &phia, const DiracSpinor &phib);
  DiracSpinor dV_ab_rhs(const DiracSpinor &phia, const DiracSpinor &phib,
                        bool conj = false);
  // DiracSpinor dV_ab_Y_rhs(const DiracSpinor &phi_alpha,
  //                         const DiracSpinor &phi_a);

  double dX_nm_bbe_rhs(const DiracSpinor &phi_n, const DiracSpinor &phi_m,
                       const DiracSpinor &phi_b, const DiracSpinor &X_beta);
  double dY_nm_bbe_rhs(const DiracSpinor &phi_n, const DiracSpinor &phi_m,
                       const DiracSpinor &phi_b, const DiracSpinor &Y_beta);

  void print() const;

private:
  std::size_t core_index(const DiracSpinor &phic);
};
