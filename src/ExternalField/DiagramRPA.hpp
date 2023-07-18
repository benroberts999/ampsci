#pragma once
#include "CorePolarisation.hpp"
#include "Coulomb/QkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>
class Wavefunction;
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}
namespace HF {
class HartreeFock;
}

namespace ExternalField {

//! RPA correction to matrix elements, using Diagram technique
class DiagramRPA : public CorePolarisation {

private:
  const HF::HartreeFock *p_hf;
  std::vector<DiracSpinor> holes{};
  std::vector<DiracSpinor> excited{};
  double eps_targ = 1.0e-10;

  // t0's never change
  // NO! They change if omega is updated (frequency dependent operator!)
  std::vector<std::vector<double>> t0am{};
  std::vector<std::vector<double>> t0ma{};
  // t's updated each solve_core itteration
  std::vector<std::vector<double>> tam{};
  std::vector<std::vector<double>> tma{};

  // Note: W's depend on rank (also parity)! Can re-use!?
  // These are probably an excellent candidate for unordered_map?
  std::vector<std::vector<std::vector<std::vector<double>>>> Wanmb{};
  std::vector<std::vector<std::vector<std::vector<double>>>> Wabmn{};
  std::vector<std::vector<std::vector<std::vector<double>>>> Wmnab{};
  std::vector<std::vector<std::vector<std::vector<double>>>> Wmban{};

  // // nb: much slower to use Qk table
  // static constexpr bool m_USE_QK = false;
  // Coulomb::QkTable m_qk{};

public:
  //! Normal constructor: needs core to split basis: only uses basis.
  DiagramRPA(const DiracOperator::TensorOperator *const h,
             const std::vector<DiracSpinor> &basis,
             const HF::HartreeFock *in_hf, const std::string &atom = "Atom");

  //! Second constructor: copies over W matrices (depend only on k/pi)
  DiagramRPA(const DiracOperator::TensorOperator *const h,
             const DiagramRPA *const drpa);

public:
  //! Itterates the RPA equations for core electrons
  virtual void solve_core(const double omega, int max_its = 200,
                          const bool print = true) override final;

  //! Calculates RPA correction to matrix element: <A||dV||B>
  virtual double dV(const DiracSpinor &Fa,
                    const DiracSpinor &Fb) const override final;

  double dV_diagram(const DiracSpinor &Fa, const DiracSpinor &Fb) const;

  //! @brief Clears the t_am and t_ma RPA ME's [RPA ME's for hole-excited]
  //! @Details If a previous run failed, can clear t_am's + re-try
  virtual void clear() override final;

  //! Updates lowest-order t_am matrix elements and resets RPA (+updates operator)
  void update_t0s(const DiracOperator::TensorOperator *const h = nullptr);

  //! Copies the tam (and tma) values across from different RPAD. If two
  //! operators are similar, this can save time on the itterations.
  void grab_tam(const DiagramRPA *const drpa) {
    tam = drpa->tam;
    tma = drpa->tma;
  }

private:
  // Note: only writes W (depends on k/pi, and basis). Do not write t's, since
  // they depend on operator. This makes it very fast when making small changes
  // to operator (don't need to re-calc W)
  // Note: doesn't depend on grid!
  bool read_write(const std::string &fname, IO::FRW::RoW rw);

  void fill_W_matrix(const DiracOperator::TensorOperator *const h);
  void setup_ts(const DiracOperator::TensorOperator *const h);

public:
  DiagramRPA &operator=(const DiagramRPA &) = delete;
  DiagramRPA(const DiagramRPA &) = default;
  ~DiagramRPA() = default;
};

} // namespace ExternalField
