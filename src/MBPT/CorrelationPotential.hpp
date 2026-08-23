#pragma once
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "MBPT/Feynman.hpp"
#include "MBPT/Goldstone.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <optional>
#include <vector>

namespace MBPT {

struct SigmaData {
  int kappa;
  double en;
  SpinorMatrix<double> Sigma;
  int n{0};
  double lambda{1.0};
  //! fk screening factors used for this Sigma (empty if none)
  std::vector<double> fk{};
};

enum class SigmaMethod { Goldstone, Feynman };

struct rgrid_params {
  double r0{1.0e-4};
  double rmax{30.0};
  std::size_t stride{4};
};

//==============================================================================

class CorrelationPotential {
  const HF::HartreeFock *m_HF;
  std::vector<DiracSpinor> m_basis; // so we can delay goldstone construction
  std::vector<SigmaData> m_Sigmas{};
  double m_r0, m_rmax;
  std::size_t m_stride;
  std::size_t m_i0, m_size; // need?

  SigmaMethod m_method;
  int m_n_min_core;
  // int m_n_min_core_F;
  bool m_includeG;
  bool m_includeBreit_b2;
  int m_n_max_breit;

  std::optional<Goldstone> m_Gold{};

  FeynmanOptions m_Foptions;
  bool m_calculate_fk; // if not, need fk and etak
  std::vector<double> m_fk;
  std::vector<double> m_etak;
  // Apply fk to both Coulomb lines of the exchange diagrams (else outer only)
  bool m_fk_both_lines;
  // Exchange diagrams by Feynman method (else Goldstone with fk)
  bool m_feynman_exchange;

  std::optional<Feynman> m_Fy{};

  // These are only for calculating fk and eta
  std::optional<Feynman> m_Fy0{};
  std::optional<Feynman> m_FyX{};
  std::optional<Feynman> m_FyH{};

  std::string m_fname{};

  // Ladder correction, Sigma_L: read from file (produced by the Ladder{}
  // block/driver - see MBPT::ladder). Stored separately from the base Sigma
  // (may have its own sub-grid and include_G); added in SigmaFv, and scaled
  // by the same lambda as the base Sigma.
  std::vector<SigmaData> m_Sigma_L{};
  std::string m_ladder_file{};

  // Energy derivative, dSigma/dE (forward difference; see m_delta_en).
  // Formed alongside Sigma when m_form_derivative is set; appended to the
  // sigma file (older files simply have none)
  bool m_form_derivative{false};
  std::vector<SigmaData> m_dSigma{};
  static constexpr double m_delta_en = 0.01;

public:
  CorrelationPotential(
    const std::string &fname, const HF::HartreeFock *vHF,
    const std::vector<DiracSpinor> &basis, double r0, double rmax,
    std::size_t stride, int n_min_core, SigmaMethod method,
    bool include_g = false, bool include_Breit_b2 = false, int n_max_breit = 0,
    const FeynmanOptions &Foptions = {}, bool calculate_fk = true,
    const std::vector<double> &fk = {}, const std::vector<double> &etak = {},
    const std::string &ladder_file = "", bool form_derivative = false,
    bool fk_both_lines = false, bool feynman_exchange = false);

  // // not thread safe!
  // void formSigma(int kappa, double en, int n = 0) {}
  // not thread safe!
  void formSigma(int kappa, double ev, int n, const DiracSpinor *Fv = nullptr);

  bool empty() const { return m_Sigmas.empty(); }

  const GMatrix *getSigma(int kappa, int n = 0) const;

  double getLambda(int kappa, int n = 0) const;

  void clear() { m_Sigmas.clear(); }

  //! returns Spinor: Sigma|Fv>
  //! @details If Sigma for kappa_v doesn't exist, returns |0>.
  DiracSpinor SigmaFv(const DiracSpinor &Fv) const;
  DiracSpinor operator()(const DiracSpinor &Fv) const { return SigmaFv(Fv); }

  //! For each valence state, prints MBPT(2) energy correction (direct,
  //! exchange, total), and <v|Sigma|v> using the stored Sigma matrix.
  void print_de(const std::vector<DiracSpinor> &valence);

  //! True if any dSigma/dE matrices are present (see form_derivative)
  bool has_derivative() const { return !m_dSigma.empty(); }

  //! Pointer to dSigma/dE data for given kappa (and n); nullptr if not present
  const SigmaData *get_derivative(int kappa, int n = 0) const;

  /*!
    @brief Returns lambda * dSigma/dE |Fv>; returns |0> if no derivative
    exists for this kappa.
    @details
    The energy derivative of Sigma, formed by finite difference alongside
    Sigma itself (option form_derivative), so it corresponds to the actual
    method used (Goldstone/Feynman, screening, etc.). Scaled by the same
    lambda as the base Sigma.

    @note The ladder correction (Sigma_L) has no energy derivative: it is
          included in SigmaFv() but not here.
  */
  DiracSpinor dSigmaFv(const DiracSpinor &Fv) const;

  //! Stores scaling factors, lambda, for each kappa (Sigma -> lamda*Sigma)
  void scale_Sigma(const std::vector<double> &lambdas);

  // if n=0, scales _all_
  void scale_Sigma(double lambda, int kappa, int n = 0);

  //! Prints the scaling factors to screen
  void print_scaling() const;

  //! Prints the scaling factors to screen
  void print_info() const;

  //! Short string identifying the method (main differences only), e.g.,
  //! "Feynman, all-order" or "Goldstone". Intended for file-cache keys.
  std::string method_string() const;

  //! Fitting (scaling) factors as a string, "kappa=value," rounded to 4 dp;
  //! empty if none scaled. For file-cache keys (e.g., the ci-file hash).
  std::string lambda_string() const;

  /*!
    @brief Average of the stored fk screening factors over the lowest Sigma
    of each l up to @p l_max; empty if none stored.
    @details
    For re-use of the Sigma_1 screening in Sigma_2 (CI), where there is one
    effective fk per Coulomb line but no unique state: the low-l factors are
    the relevant ones for the CI valence space (higher-l factors can differ
    significantly, e.g., f/g at k = 0, but such states rarely contribute).
    Weighted by l, not kappa: the fine-structure pair of each l is averaged
    first, then the mean is taken over the ls (so s counts the same as p,
    d, ...). Truncated to the shortest stored list.

    @param l_max  Maximum l included in the average (CI uses 2: s, p, d).
  */
  std::vector<double> average_fk(int l_max = 2) const;

  //! Prints the sub-grid parameters to screen
  void print_subGrid() const;

  void write(const std::string &fname) { read_write(fname, IO::FRW::write); }

  //! Pointer to the stored Sigma data (matrix, energy formed at, lambda) for
  //! given kappa (and n); nullptr if not present
  const SigmaData *get(int kappa, int n = 0) const;

private:
  bool read_write(const std::string &fname, IO::FRW::RoW rw);
  void setup_Feynman();
  std::vector<double> calculate_fk(double ev, const DiracSpinor &v) const;
  std::vector<double> calculate_etak(double ev, const DiracSpinor &v) const;
  const SigmaData *get_ladder(int kappa, int n = 0) const;

  // given_fk: screening factors to use (from state_fk); if nullptr,
  // calculated internally (if applicable). print = false: silent (e.g., the
  // extra evaluation for the derivative - only the used Sigma is reported)
  GMatrix formSigma_F(int kappa, double ev, const DiracSpinor *Fv = nullptr,
                      const std::vector<double> *given_fk = nullptr,
                      bool print = true);
  GMatrix formSigma_G(int kappa, double ev, const DiracSpinor *Fv = nullptr,
                      bool print = true);

  // Calculates the fk screening factors for this state, once (with print;
  // stores first set in m_fk). nullopt if not applicable (not
  // Feynman-screening, fk given manually, or no Fv)
  std::optional<std::vector<double>> state_fk(double ev, const DiracSpinor *Fv);

  // Forms dSigma/dE for given kappa by forward finite difference (step
  // m_delta_en): one extra Sigma evaluation at ev + delta, re-using the
  // base Sigma matrix and the same fk; stores in m_dSigma
  void form_derivative(int kappa, double ev, int n, const DiracSpinor *Fv,
                       const GMatrix &Sigma0,
                       const std::vector<double> *given_fk = nullptr);

public:
  CorrelationPotential &operator=(const CorrelationPotential &) = default;
  CorrelationPotential(const CorrelationPotential &) = default;
  ~CorrelationPotential() = default;

  //
};

} // namespace MBPT
