#pragma once
#include "CorrelationPotential.hpp"
#include "Maths/Grid.hpp"
#include <memory>
#include <string>
#include <vector>
// class Grid;
namespace HF {
class HartreeFock;
}

namespace MBPT {

enum class States { core, excited, both };

enum class GrMethod { Green, basis };

inline std::string_view ParseEnum(GrMethod method) {
  switch (method) {
  case GrMethod::Green:
    return "Green";
  case GrMethod::basis:
    return "basis";
  }
  return "unkown";
}

enum class ExchangeMethod { Goldstone, w1w2, w1, none };

inline std::string_view ParseEnum(ExchangeMethod method) {
  switch (method) {
  case ExchangeMethod::Goldstone:
    return "Goldstone";
  case ExchangeMethod::w1w2:
    return "w1w2";
  case ExchangeMethod::w1:
    return "w1";
  case ExchangeMethod::none:
    return "none";
  }
  return "unkown";
}

//==============================================================================

class FeynmanSigma final : public CorrelationPotential {
public:
  FeynmanSigma(const HF::HartreeFock *const in_hf,
               const std::vector<DiracSpinor> &basis, const Sigma_params &sigp,
               const rgrid_params &subgridp,
               // const std::vector<DiracSpinor> &valence,
               const std::string &atom, bool verbose = true,
               ExchangeMethod m_ex_method = ExchangeMethod::Goldstone);

  FeynmanSigma &operator=(const FeynmanSigma &) = delete;
  FeynmanSigma(const FeynmanSigma &) = delete;
  ~FeynmanSigma() = default;

  void formSigma(int kappa, double en, int n = 0) override final;

protected:
public:
  //!@brief
  // Calculates (radial) Hartree-Fock Greens function G_kappa(er + i*ei).
  // states = States::{core, excited, both}; method = GrMethod::{Green, basis}
  [[nodiscard]] ComplexGMatrix Green(int kappa, ComplexDouble en,
                                     States states = States::both,
                                     GrMethod method = GrMethod::Green) const;

  // Make this private?
  //! Takes Gr = G(e_r) and e_i, returns G(e_r + i e_i). nb: Must be FULL green
  [[nodiscard]] ComplexGMatrix GreenAtComplex(const ComplexGMatrix &Gr,
                                              double e_imag) const;

  //! Calculates radial polarisation operator
  [[nodiscard]] ComplexGMatrix Polarisation_k(int k, ComplexDouble omega,
                                              GrMethod method) const;

  //! Returns (reference to) q^k (radial) matrix. Note: includes dri*drj
  [[nodiscard]] const ComplexGMatrix &get_qk(int k) const;

  //! Returns (ref to) radial exchange matrix Vx_kappa. Nb: includes dri*drj
  [[nodiscard]] const GMatrix &get_Vx_kappa(int kappa) const;

  [[nodiscard]] const ComplexGMatrix &get_dri() const { return *m_dri; }
  [[nodiscard]] const ComplexGMatrix &get_drj() const { return *m_drj; }

  //! Calculates direct Sigma using Feynman method
  [[nodiscard]] GMatrix FeynmanDirect(int kv, double env, int k = -1) const;

  [[nodiscard]] GMatrix FeynmanEx_w1w2(int kv, double en) const;
  //! Calculates exchange Sigma using Feynman method [w_1 version]
  [[nodiscard]] GMatrix FeynmanEx_1(int kv, double env) const;

  //! Contructs G_a Green-like fn for single state: returns f*|ket><bra|
  [[nodiscard]] ComplexGMatrix G_single(const DiracSpinor &ket,
                                        const DiracSpinor &bra,
                                        const ComplexDouble f) const;

private:
  // Calculates initial data needed for Feynman (|a><a|, w grids etc)
  void prep_Feynman();
  // Calculates + stores the (radial) q^k matrix, for each k (includes dri*drj)
  void form_Q_dr();
  // Calculates and stores radial exchange matrix Vx for each kappa
  void form_Vx();
  // Forms Vx for given kappa, exhange operator matrix (includes dri,drj)
  [[nodiscard]] GMatrix calculate_Vx_kappa(int kappa) const;
  GMatrix calculate_Vhp(const DiracSpinor &Fa) const;

  // Calculates and stores radial projection operators for core state |a><a|
  void form_Pa_core();
  // Sets up imaginary frequency (omega) grids for integrations
  void setup_omega_grid();

  // Calculate "Core" Green fn by direct summation over only core states
  [[nodiscard]] ComplexGMatrix Green_core(int kappa, ComplexDouble en) const;
  // Calculates "excited" Greens function, by: G_ex = G - G_core
  [[nodiscard]] ComplexGMatrix Green_ex(int kappa, ComplexDouble en,
                                        GrMethod method = GrMethod::Green,
                                        const DiracSpinor *Fc_hp = nullptr,
                                        int k_hp = 0) const;

  // Force Gk to be orthogonal to the core states
  void makeGOrthogCore(ComplexGMatrix *Gk, int kappa) const;

  // // Calculates Hartree-Fock Green function (including exchange), for real en
  // [[nodiscard]] ComplexGMatrix Green_hf_real(int kappa, double en) const;
  // Calculates HF Green function (including exchange), for Complex en
  [[nodiscard]] ComplexGMatrix Green_hf(int kappa, ComplexDouble en,
                                        const DiracSpinor *Fc_hp = nullptr,
                                        int k_hp = 0) const;

  // Calculate HF Greens function (complex en), using basis expansion
  [[nodiscard]] ComplexGMatrix Green_hf_basis(int kappa, ComplexDouble en,
                                              bool ex_only = false) const;
  // Contructs G0(w) (no exchange) Greens fn; x0,xI homog. solns reg at 0,infty
  [[nodiscard]] GMatrix MakeGreensG0(const DiracSpinor &x0,
                                     const DiracSpinor &xI,
                                     const double w) const;

  // X = [1 + i pi*q]^{-1}
  ComplexGMatrix X_PiQ(const ComplexGMatrix &pik,
                       const ComplexGMatrix &qk) const;

  std::vector<std::vector<ComplexGMatrix>> make_pi_wk(int max_k,
                                                      GrMethod pol_method,
                                                      double omre,
                                                      const Grid &wgrid) const;

  // ComplexGMatrix form_QPQ_wk(const ComplexGMatrix &PiQ) const;
  std::vector<std::vector<ComplexGMatrix>> form_QPQ_wk(int max_k,
                                                       GrMethod pol_method,
                                                       double en_re,
                                                       const Grid &wgrid) const;
  std::vector<std::vector<ComplexGMatrix>>
  form_Greens_kapw(int max_kappa_index, GrMethod method, double omre,
                   const Grid &wgrid) const;

  // exchange: sum_kl gA*qk*ql*(c1 * pa*gxBm + c2 * gxBp*pa)
  [[nodiscard]] GMatrix
  sumkl_GQPGQ(const ComplexGMatrix &gA, const ComplexGMatrix &gxBm,
              const ComplexGMatrix &gxBp, const ComplexGMatrix &pa, int kv,
              int kA, int kB, int ka,
              const std::vector<ComplexGMatrix> *const = nullptr) const;

  [[nodiscard]] GMatrix sumkl_gqgqg(const ComplexGMatrix &gA,
                                    const ComplexGMatrix &gB,
                                    const ComplexGMatrix &gG, int kv, int kA,
                                    int kB, int kG, int kmax) const;

  double Lkl_abcd(int k, int l, int ka, int kb, int kc, int kd) const;

  // result += Real{ sum_ij [ factor * a1j * bij * cj2 * (d_1i * e_i2) ] }
  void tensor_5_product(GMatrix *result, const ComplexDouble &factor,
                        const ComplexGMatrix &a, const ComplexGMatrix &b,
                        const ComplexGMatrix &c, const ComplexGMatrix &d,
                        const ComplexGMatrix &e) const;

  // Better solution than this!
  GMatrix Exchange_Goldstone(const int kappa, const double en,
                             const std::vector<double> &v_fk) const;

private:
  const bool m_screen_Coulomb;
  const bool m_holeParticle;

  const double m_omre;
  const double m_w0;
  const double m_w_ratio;

  const GrMethod m_Green_method;
  const GrMethod m_Pol_method;

  const HF::HartreeFock *const p_hf;
  const int m_min_core_n;
  int m_max_kappaindex_core;
  int m_max_kappaindex;

  std::vector<ComplexGMatrix> m_qhat{};
  std::vector<ComplexGMatrix> m_Pa{}; // |a><a| for each core state
  std::vector<GMatrix> m_Vxk{};       // one each kappa in core

  std::unique_ptr<ComplexGMatrix> m_dri = nullptr;
  std::unique_ptr<ComplexGMatrix> m_drj = nullptr;
  std::unique_ptr<Grid> m_wgridD = nullptr;
  // only use every nth point on Im(w) grid for exchange
  std::size_t m_wX_stride{1}; // XXX input?

  std::vector<std::vector<ComplexGMatrix>> m_qpq_wk{};

  int m_k_cut = 10; // XXX Make input?

  const bool m_print_each_k = false;

  ExchangeMethod m_ex_method;
  // ExchangeMethod m_ex_method = ExchangeMethod::Goldstone;
  // ExchangeMethod m_ex_method = ExchangeMethod::w1;
  // ExchangeMethod m_ex_method = ExchangeMethod::none;
  // ExchangeMethod m_ex_method = ExchangeMethod::w1w2;

  // for calculating fk (if not explicitely given)
  std::unique_ptr<const FeynmanSigma> m_S0{nullptr};
  std::unique_ptr<const FeynmanSigma> m_Sx{nullptr};

  std::string m_ladder_file;
  bool m_verbose;
};

} // namespace MBPT
