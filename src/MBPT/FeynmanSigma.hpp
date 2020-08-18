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

//! Many-body perturbation theory
namespace MBPT {

enum class States { core, excited, both };
enum class GrMethod { Green, basis };

//******************************************************************************
/*!
@brief
@details
*/
class FeynmanSigma final : public CorrelationPotential {
public:
  FeynmanSigma(const HF::HartreeFock *const in_hf,
               const std::vector<DiracSpinor> &basis, const Sigma_params &sigp,
               const rgrid_params &subgridp, const std::vector<double> &en_list,
               const std::string &atom);

  FeynmanSigma &operator=(const FeynmanSigma &) = delete;
  FeynmanSigma(const FeynmanSigma &) = delete;
  ~FeynmanSigma() = default;

protected:
  // Fills the correlation potential for given kappa
  void fill_Sigma_k(GMatrix *Sigma, const int kappa, const double en) final;

public:
  //!@brief
  // Calculates (radial) Hartree-Fock Greens function G_kappa(er + i*ei).
  // states = States::{core, excited, both}; method = GrMethod::{Green, basis}
  ComplexGMatrix Green(int kappa, double en_re, double en_im = 0.0,
                       States states = States::both,
                       GrMethod method = GrMethod::Green) const;

  //! Takes Gr = G(e_r) and e_i, returns G(e_r + i e_i). nb: Must be FULL green
  ComplexGMatrix GreenAtComplex(const ComplexGMatrix &Gr, double e_imag) const;

  //! Calculates polarisation operator, all core states
  ComplexGMatrix Polarisation(int ka, int kA, double om_re, double om_im,
                              GrMethod method = GrMethod::basis) const;

  //! Returns (reference to) q^k (radial) matrix. Note: includes dri*drj
  const ComplexGMatrix &get_qk(int k) const;

  //! Screens Coulomb: q_scr = q * [1-pi*q]^-1
  ComplexGMatrix screenedCoulomb(const ComplexGMatrix &q,
                                 const ComplexGMatrix &pi) const;

  //! Returns (ref to) radial exchange matrix Vx_kappa. Nb: includes dri*drj
  const GMatrix &get_Vx_kappa(int kappa) const;

  const ComplexGMatrix &get_dri() const { return *m_dri; }
  const ComplexGMatrix &get_drj() const { return *m_drj; }

  //! Calculates direct Sigma using Feynman method
  GMatrix FeynmanDirect(int kv, double env);
  //! Calculates exchange Sigma using Feynman method [w_1 version]
  GMatrix FeynmanEx_1(int kv, double env);

  // Contructs G_a Green-like fn for single state: returns f*|ket><bra|
  ComplexGMatrix G_single(const DiracSpinor &ket, const DiracSpinor &bra,
                          const ComplexDouble f) const;

private:
  // Calculates initial data needed for Feynman (|a><a|, w grids etc)
  void prep_Feynman();
  // Calculates + stores the (radial) q^k matrix, for each k (includes dri*drj)
  void form_Q_dr();
  // Calculates and stores radial exchange matrix Vx for each kappa
  void form_Vx();
  // Forms Vx for given kappa, exhange operator matrix (includes dri,drj)
  GMatrix calculate_Vx_kappa(int kappa) const;
  // Calculates and stores radial projection operators for core states |a><a|
  void form_Pa_core();
  // Sets up imaginary frequency (omega) grids for integrations
  void setup_omega_grid();

  // Calculate "Core" Green fn by direct summation over only core states
  ComplexGMatrix Green_core(int kappa, double en_re, double en_im) const;
  // Calculates "excited" Greens function, by: G_ex = G - G_core
  ComplexGMatrix Green_ex(int kappa, double en_re, double en_im,
                          GrMethod method = GrMethod::Green) const;
  // Calculates Hartree-Fock Green function (including exchange), for real en
  ComplexGMatrix Green_hf(int kappa, double en) const;
  // Calculates HF Green function (including exchange), for Complex en
  ComplexGMatrix Green_hf(int kappa, double en_re, double en_imag) const;
  // Calculate HF Greens function (complex en), using basis expansion
  ComplexGMatrix Green_hf_basis(int kappa, double en_re, double en_im,
                                bool ex_only = false) const;
  // Contructs G0(w) (no exchange) Greens fn; x0,xI homog. solns reg at 0,infty
  GMatrix MakeGreensG0(const DiracSpinor &x0, const DiracSpinor &xI,
                       const double w) const;

  // Takes a G(e) and de_r, de_i, returns G(e + de_r + i*de_i) - needs testing
  ComplexGMatrix GreenAtComplexShift(const ComplexGMatrix &Gr, double de_re,
                                     double om_imag) const;

  // Calculates polarisation operator, only single core state (pa = |a><a|)
  ComplexGMatrix Polarisation_a(const ComplexGMatrix &pa, double ena,
                                int k_alpha, double om_re, double om_im,
                                GrMethod method = GrMethod::basis) const;

  // direct: sum_k [ck qk * pi(w) * qk], ck angular factor
  GMatrix sumk_cGQPQ(int kv, int ka, int kalpha, int kbeta,
                     const ComplexGMatrix &g_beta,
                     const ComplexGMatrix &pi_aalpha) const;
  // exchange: sum_kl gA*qk*ql*(c1 * pa*gxBm + c2 * gxBp*pa)
  GMatrix sumkl_GQPGQ(const ComplexGMatrix &gA, const ComplexGMatrix &gxBm,
                      const ComplexGMatrix &gxBp, const ComplexGMatrix &pa,
                      int kv, int kA, int kB, int ka) const;

  // result += Real{ sum_ij [ factor * a1j * bij * cj2 * (d_1i * e_i2) ] }
  void tensor_5_product(GMatrix *result, const ComplexDouble &factor,
                        const ComplexGMatrix &a, const ComplexGMatrix &b,
                        const ComplexGMatrix &c, const ComplexGMatrix &d,
                        const ComplexGMatrix &e) const;

private:
  const bool screen_Coulomb;
  const double m_omre;
  const bool basis_for_Green;
  const bool basis_for_Pol;
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
  std::unique_ptr<Grid> m_wgridX = nullptr;

  const bool exclude_exchange = true; // for testing!
};

} // namespace MBPT
