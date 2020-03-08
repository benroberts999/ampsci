#pragma once
#include "Angular/Angular.hpp"
#include <vector>
class DiracSpinor;
class DiracOperator;
class Wavefunction;

enum class dPsiType { X, Y }; // add conj?

class ExternalField {
public:
  ExternalField(const DiracOperator *const h,
                const std::vector<DiracSpinor> &core,
                const std::vector<double> &vl, const double alpha);

  ExternalField &operator=(const ExternalField &) = delete;
  ExternalField(const ExternalField &) = default;
  ~ExternalField() = default;

private:
  // dPhi = X exp(-iwt) + Y exp(+iwt)
  // (H - e - w)X = -(h + dV - de)Phi
  // (H - e + w)Y = -(h* + dV* - de)Phi
  // X_c = sum_x X_x,
  // j(x)=j(c)-k,...,j(c)+k.  And: pi(x) = pi(c)*pi(h)
  std::vector<std::vector<DiracSpinor>> m_X = {};
  std::vector<std::vector<DiracSpinor>> m_Y = {};
  // can just write these to disk! Read them in, continue as per normal

  const DiracOperator *const m_h; //??
  const std::vector<DiracSpinor> *const p_core;
  const std::vector<double> m_vl; // Add H_mag ?
  const double m_alpha;
  const int m_rank;
  const int m_pi;
  const bool m_imag;

  // Angular::SixJ m_6j; // used?

public:
  void solve_TDHFcore(const double omega, int max_its = 100,
                      const bool print = true);
  void solve_TDHFcore_matrix(const Wavefunction &wf, const double omega,
                             const int max_its = 25);
  void reZero();

  // does it matter if a or b is in the core?
  double dV_ab(const DiracSpinor &Fa, const DiracSpinor &Fb, bool conj,
               const DiracSpinor *const Fexcl = nullptr) const;
  double dV_ab(const DiracSpinor &Fa, const DiracSpinor &Fb) const;

  DiracSpinor dV_ab_rhs(const DiracSpinor &Fa, const DiracSpinor &Fb,
                        bool conj = false,
                        const DiracSpinor *const Fexcl = nullptr) const;

  // make const? Private?
  const std::vector<DiracSpinor> &get_dPsis(const DiracSpinor &Fc,
                                            dPsiType XorY) const;
  const DiracSpinor &get_dPsi_x(const DiracSpinor &Fc, dPsiType XorY,
                                const int kappa_x) const;

  // For "full" matrix version: not ready...
  double dX_nm_bbe_rhs(const DiracSpinor &Fn, const DiracSpinor &Fm,
                       const DiracSpinor &Fb, const DiracSpinor &X_beta) const;
  double dY_nm_bbe_rhs(const DiracSpinor &Fn, const DiracSpinor &Fm,
                       const DiracSpinor &Fb, const DiracSpinor &Y_beta) const;

  void print() const;

private:
  std::size_t core_index(const DiracSpinor &Fc) const;
};
