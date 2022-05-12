#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include <utility>
#include <vector>

namespace HF {

namespace hidden {
struct Breit_Bk_ba; // forward decl
}

//==============================================================================
//! Breit (Hartree-Fock Breit) interaction potential
class Breit {
public:
  //! Contains ptr to core: careful if updating core (e.g., in HF)
  Breit(const std::vector<DiracSpinor> &in_core, double in_scale = 1.0)
      : p_core(&in_core), m_scale(in_scale) {}
  // nb: During HF, must update orbitals each time

  //! () operator: returns VbrFa(Fa)
  DiracSpinor operator()(const DiracSpinor &Fa) const { return VbrFa(Fa); }
  //! Calculates V_br*Fa = \sum_b\sum_k B^k_ba F_b [Breit part of HF-Breit pot.]
  DiracSpinor VbrFa(const DiracSpinor &Fa) const;

  //! dV_b*Fa, dV_b is the exchange RPA correction arising due to Fb -> Fb + dFb
  //! @details
  //! K is multipolarity of RPA operator, Fb is core state, with Xbeta and Ybeta
  //! perturbations. "reduced rhs"
  DiracSpinor dVbrX_Fa(int kappa, int K, const DiracSpinor &Fa,
                       const DiracSpinor &Fb, const DiracSpinor &Xbeta,
                       const DiracSpinor &Ybeta) const;
  //! Direct RPA correction arising due to Fb -> Fb + dFb
  DiracSpinor dVbrD_Fa(int kappa, int K, const DiracSpinor &Fa,
                       const DiracSpinor &Fb, const DiracSpinor &Xbeta,
                       const DiracSpinor &Ybeta) const;

private:
  const std::vector<DiracSpinor> *const p_core;
  const double m_scale;

  // Calculates \sum_k B^k_ba F_b (single core contr. to V_brFa)
  void BkbaFb(DiracSpinor *BFb, const DiracSpinor &Fa,
              const DiracSpinor &Fb) const;

  void MOPk_ij_Fc(DiracSpinor *BFc, const double Cang,
                  const hidden::Breit_Bk_ba &Bkab, int k, int ki, int kj,
                  const DiracSpinor &Fc) const;
  void Nk_ij_Fc(DiracSpinor *BFc, const double Cang,
                const hidden::Breit_Bk_ba &Bkij, int k, int ki, int kj,
                const DiracSpinor &Fc) const;

  // Internal angular coefs
  double eta(int k, int ka, int kb) const;
  std::pair<double, double> Mk(int k) const;
  double Nkba(int k, int kb, int ka) const;
  std::pair<double, double> Ok(int k) const;
  double Pk(int k) const;

public:
  Breit &operator=(const Breit &) = delete;
  Breit(const Breit &) = delete;
  ~Breit() = default;
};

//==============================================================================
namespace hidden {

struct Breit_Bk_ba {
  // Class to hold the Breit-Coulomb integrals
public:
  Breit_Bk_ba(const DiracSpinor &Fb, const DiracSpinor &Fa);

  const std::size_t max_k;
  std::vector<std::vector<double>> bk_0{};
  std::vector<std::vector<double>> bk_inf{};
  std::vector<std::vector<double>> gk_0{};
  std::vector<std::vector<double>> gk_inf{};
};
} // namespace hidden

} // namespace HF
