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
};

enum class SigmaMethod { Goldstone, Feynman };

struct rgrid_params {
  double r0{1.0e-4};
  double rmax{30.0};
  std::size_t stride{4};
};

// method:
/*

  - Goldstone/Feynman
    - Screening
    - hp
  - Include g?

  - Effective screening: input, or calculate?
    - If given, use, else: calculate


*/

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

  std::optional<Feynman> m_Fy{};

  // These are only for calculating fk and eta
  std::optional<Feynman> m_Fy0{};
  std::optional<Feynman> m_FyX{};
  std::optional<Feynman> m_FyH{};

public:
  CorrelationPotential(const std::string &fname, const HF::HartreeFock *vHF,
                       const std::vector<DiracSpinor> &basis, double r0,
                       double rmax, std::size_t stride, int n_min_core,
                       SigmaMethod method, bool include_g = false,
                       bool include_Breit_b2 = false, int n_max_breit = 0,
                       const FeynmanOptions &Foptions = {},
                       bool calculate_fk = true,
                       const std::vector<double> &fk = {},
                       const std::vector<double> &etak = {});

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

  //! Stores scaling factors, lambda, for each kappa (Sigma -> lamda*Sigma)
  void scale_Sigma(const std::vector<double> &lambdas);

  // if n=0, scales _all_
  void scale_Sigma(double lambda, int kappa, int n = 0);

  //! Prints the scaling factors to screen
  void print_scaling() const;

  //! Prints the scaling factors to screen
  void print_info() const;

  //! Prints the sub-grid parameters to screen
  void print_subGrid() const;

  void write(const std::string &fname) { read_write(fname, IO::FRW::write); }

private:
  bool read_write(const std::string &fname, IO::FRW::RoW rw);
  void setup_Feynman();
  std::vector<double> calculate_fk(double ev, const DiracSpinor &v) const;
  std::vector<double> calculate_etak(double ev, const DiracSpinor &v) const;
  const SigmaData *get(int kappa, int n = 0) const;

  GMatrix formSigma_F(int kappa, double ev, const DiracSpinor *Fv = nullptr);
  GMatrix formSigma_G(int kappa, double ev, const DiracSpinor *Fv = nullptr);

public:
  CorrelationPotential &operator=(const CorrelationPotential &) = default;
  CorrelationPotential(const CorrelationPotential &) = default;
  ~CorrelationPotential() = default;

  //
};

} // namespace MBPT
