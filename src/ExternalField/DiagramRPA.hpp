#pragma once
#include "CorePolarisation.hpp"
#include "Coulomb/QkTable.hpp"
#include "HF/Breit.hpp"
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
  std::optional<HF::Breit> m_Br{};
  std::vector<DiracSpinor> m_holes{};
  std::vector<DiracSpinor> m_excited{};

  // t0 (no RPA). These stay constant for frequency-independent operators
  // But change for f-dependent. Set during solve_core
  std::vector<std::vector<double>> m_t0am{};
  std::vector<std::vector<double>> m_t0ma{};
  // t's updated each solve_core itteration
  std::vector<std::vector<double>> m_tam{};
  std::vector<std::vector<double>> m_tma{};

  // Note: W's depend on rank (also parity)! Can re-use!?
  // These are probably an excellent candidate for unordered_map?
  std::vector<std::vector<std::vector<std::vector<double>>>> m_Wanmb{};
  std::vector<std::vector<std::vector<std::vector<double>>>> m_Wabmn{};
  std::vector<std::vector<std::vector<std::vector<double>>>> m_Wmnab{};
  std::vector<std::vector<std::vector<std::vector<double>>>> m_Wmban{};

public:
  //! Normal constructor: needs core to split basis: only uses basis.
  DiagramRPA(const DiracOperator::TensorOperator *const h,
             const std::vector<DiracSpinor> &basis,
             const HF::HartreeFock *in_hf, const std::string &atom = "",
             bool print = true);

  //! Second constructor: copies over W matrices (depend only on k/pi)
  [[deprecated]]
  DiagramRPA(const DiracOperator::TensorOperator *const h,
             const DiagramRPA *const drpa);

public:
  //! Itterates the RPA equations for core electrons
  void solve_core(const double omega, int max_its = 200,
                  const bool print = true) override final;

  //! Returns RPA method
  Method method() const override final { return Method::diagram; }

  //! Calculates RPA correction to matrix element: <A||dV||B>
  double dV(const DiracSpinor &Fa, const DiracSpinor &Fb) const override final;

  //! @brief Clears the t_am and t_ma RPA ME's [RPA ME's for hole-excited]
  //! @Details If a previous run failed, can clear t_am's + re-try
  void clear() override final;

  //! Updates lowest-order t_am matrix elements and resets RPA (+updates operator)
  [[deprecated]]
  void update_t0s(const DiracOperator::TensorOperator *const h = nullptr);

  //! Copies the tam (and tma) values across from different RPAD. If two
  //! operators are similar, this can save time on the itterations.
  [[deprecated]]
  void grab_tam(const DiagramRPA *const drpa) {
    m_tam = drpa->m_tam;
    m_tma = drpa->m_tma;
  }

private:
  // Note: only writes W (depends on k/pi, and basis). Do not write t's, since
  // they depend on operator. This makes it very fast when making small changes
  // to operator (don't need to re-calc W)
  // Note: doesn't depend on grid!
  bool read_write(const std::string &fname, IO::FRW::RoW rw);

  // Calculates all required W^k integrals (uses rank, parity of h)
  void fill_W_matrix(const DiracOperator::TensorOperator *const h, bool print);
  // Calculates all t_am matrix elements (without RPA)
  void setup_ts(const DiracOperator::TensorOperator *const h);

public:
  DiagramRPA &operator=(const DiagramRPA &) = delete;
  DiagramRPA(const DiagramRPA &) = default;
  ~DiagramRPA() = default;
};

} // namespace ExternalField
