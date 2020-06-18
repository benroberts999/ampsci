#pragma once
#include "CorrelationPotential.hpp"
#include <memory>
#include <string>
#include <vector>
class Grid;
namespace HF {
class HartreeFock;
}

//! Many-body perturbation theory
namespace MBPT {

//******************************************************************************
/*!
@brief
@details
*/
class FeynmanSigma : public CorrelationPotential {
public:
  FeynmanSigma(const HF::HartreeFock *const in_hf,
               const std::vector<DiracSpinor> &basis, const Sigma_params &sigp,
               const rgrid_params &subgridp, const std::vector<double> &en_list,
               const std::string &atom);

  FeynmanSigma &operator=(const FeynmanSigma &) = delete;
  FeynmanSigma(const FeynmanSigma &) = delete;
  ~FeynmanSigma() = default;

protected:
  void fill_Sigma_k(GMatrix *Sigma, const int kappa, const double en) override;

private:
  // Calculate "Core" Green fn by direct summation over only core states
  ComplexGMatrix Green_core(int kappa, double en_re, double en_im) const;
  // Calculates "excited" Greens function, by: G_ex = G - G_core
  ComplexGMatrix Green_ex(int kappa, double en_re, double en_im) const;
  // Calculates Hartree-Fock Green function (including exchange), for real en
  ComplexGMatrix Green_hf(int kappa, double en) const;
  // Calculates HF Green function (including exchange), for Complex en
  ComplexGMatrix Green_hf(int kappa, double en_re, double en_imag) const;
  // Calculate HF Greens function (complex en), using basis expansion
  ComplexGMatrix Green_hf_basis(int kappa, double en_re, double en_im,
                                bool ex_only) const;
  // Takes a G(e_r) and e_i, returns G(e_r + i*e_i)
  ComplexGMatrix ComplexG(const ComplexGMatrix &Gr, double om_imag) const;
  // Takes a G(e) and de_r, de_i, returns G(e + de_r + i*de_i) - needs testing
  ComplexGMatrix ComplexG(const ComplexGMatrix &Gr, double de_re,
                          double om_imag) const;
  // Calculates polarisation operator, all core states
  ComplexGMatrix Polarisation(int kappa_a, int kappa_alpha, double om_re,
                              double om_im) const;
  // Calculates polarisation operator, only single core state (pa = |a><a|)
  ComplexGMatrix Polarisation_a(const ComplexGMatrix &pa, double ena,
                                int k_alpha, double om_re, double om_im) const;
  // Screens Coulomb: q_scr = q * [1-pi*q]^-1
  ComplexGMatrix screenedCoulomb(const ComplexGMatrix &q,
                                 const ComplexGMatrix &pi) const {
    // not checked!
    // XXX Are the dr_i, dr_j parts correct??
    return q * ((-1.0 * pi * q).plusIdent(1.0).invert());
  }

  // Calculates initial data needed for Feynman (|a><a|, w grids etc)
  void prep_Feynman();

  // Contructs G0(w) (no exchange) Greens fn; x0,xI homog. solns reg at 0,infty
  GMatrix MakeGreensG(const DiracSpinor &x0, const DiracSpinor &xI,
                      const double w) const;

  // Contructs G_a Green-like fn for single state: returns f*|ket><bra|
  ComplexGMatrix G_single(const DiracSpinor &ket, const DiracSpinor &bra,
                          const ComplexDouble f) const;
  // Forms Vx, exhange operator matrix (includes dri,drj)
  GMatrix Make_Vx(int kappa) const;

  // Calculates direct Sigma using Feynman method
  GMatrix FeynmanDirect(int kv, double env);
  // Calculates exchange Sigma using Feynman method [w_1 version]
  GMatrix FeynmanEx_1(int kv, double env);

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
};

} // namespace MBPT
