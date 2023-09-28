#pragma once
#include "CorePolarisation.hpp"
#include "Coulomb/YkTable.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>

namespace ExternalField {

//! RPA correction to matrix elements, using Diagram technique
class DiagramRPA0_jL : public CorePolarisation {

private:
  const HF::HartreeFock *p_hf;
  std::vector<DiracSpinor> holes{};
  std::vector<DiracSpinor> excited{};
  Coulomb::YkTable m_yk{};

public:
  DiagramRPA0_jL(const DiracOperator::TensorOperator *const h,
                 const std::vector<DiracSpinor> &basis,
                 const HF::HartreeFock *in_hf, int max_l)
      : CorePolarisation(h), p_hf(in_hf) {

    if (p_hf == nullptr || h == nullptr) {
      std::cout << "\nFAIL:25 in DiagramRPA0_jL - hf cannot be null\n"
                << std::flush;
    }

    // Set up basis:
    const auto &core = p_hf->core();
    for (const auto &Fi : basis) {
      const bool inCore = std::find(cbegin(core), cend(core), Fi) != cend(core);
      if (inCore) {
        holes.push_back(Fi);
      } else {
        excited.push_back(Fi);
      }
    }

    m_yk.calculate(basis);
    // need Ck up to continuum states; depends on max_l
    const auto max_tj = DiracSpinor::max_tj(basis) + 2 * max_l;
    m_yk.extend_angular(max_tj); //??
  }

public:
  //! Itterates the RPA equations for core electrons
  virtual void solve_core(const double, int, const bool) override final {
    std::cout << "ERROR: solve_core() not implemented\n";
  }

  //! Returns RPA method
  virtual Method method() const override final { return Method::Error; }

  //! Calculates RPA correction to matrix element: <A||dV||B>
  virtual double dV(const DiracSpinor &,
                    const DiracSpinor &) const override final {
    std::cout << "\nWARNING: not implemented. Use dV_diagram_jL()\n"
              << std::flush;
    std::abort();
    return 0.0;
  }

  //----------------------------------------------------------------------------
  // nb: Fv MUST appear within basis!
  double dV_diagram_jL(const DiracSpinor &Fw, const DiracSpinor &Fv,
                       const DiracOperator::jL *jl, std::size_t L,
                       double q) const {

    if (holes.empty() || excited.empty())
      return 0.0;

    if (Fv.en() > Fw.en()) {
      const auto sj = ((Fv.twoj() - Fw.twoj()) % 4 == 0) ? 1 : -1;
      const auto si = m_imag ? -1 : 1;
      return (sj * si) * dV_diagram_jL(Fv, Fw, jl, L, q);
    }

    const auto ww = 0.0;
    const auto iL = int(L);

    const auto f = (1.0 / (2 * iL + 1));

    std::vector<double> sum_a(holes.size());
#pragma omp parallel for
    for (std::size_t ia = 0; ia < holes.size(); ia++) {
      const auto &Fa = holes[ia];
      const auto s1 = ((Fa.twoj() - Fw.twoj() + 2 * iL) % 4 == 0) ? 1 : -1;
      for (std::size_t im = 0; im < excited.size(); im++) {
        const auto &Fm = excited[im];
        // if (t0am[ia][im] == 0.0)
        if (jl->is_zero(Fa, Fm, L))
          continue;
        const auto s2 = ((Fa.twoj() - Fm.twoj()) % 4 == 0) ? 1 : -1;
        // const auto Wwmva = Coulomb::Wk_abcd(Fw, Fm, Fv, Fa, L);
        // const auto Wwavm = Coulomb::Wk_abcd(Fw, Fa, Fv, Fm, L);
        const auto Wwmva = m_yk.W(iL, Fw, Fm, Fv, Fa);
        const auto Wwavm = m_yk.W(iL, Fw, Fa, Fv, Fm);
        auto ttam = jl->rme(Fa, Fm, L, q);

        const auto ttma = jl->symm_sign(Fa, Fm) * ttam;
        const auto A = ttam * Wwmva / (Fa.en() - Fm.en() - ww);
        const auto B = Wwavm * ttma / (Fa.en() - Fm.en() + ww);
        sum_a[ia] += s1 * (A + s2 * B);
      }
    }
    return f * std::accumulate(begin(sum_a), end(sum_a), 0.0);
  }

  void clear() override final { return; }

public:
  DiagramRPA0_jL &operator=(const DiagramRPA0_jL &) = delete;
  DiagramRPA0_jL(const DiagramRPA0_jL &) = default;
  ~DiagramRPA0_jL() = default;
};

} // namespace ExternalField
