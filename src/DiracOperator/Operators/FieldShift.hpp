#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace DiracOperator {

//! @brief Field shift operator, F_fs, defined dV/<r^2>
/*! @details 
*/
class F_fs final : public ScalarOperator {

private:
  double m_dr2;

public:
  F_fs(const Grid &rgrid, const Nuclear::Nucleus &nuc1,
       const Nuclear::Nucleus &nuc2)
      : ScalarOperator(
            Parity::even,
            PhysConst::Hartree_GHz /
                (nuc2.r_rms() * nuc2.r_rms() - nuc1.r_rms() * nuc1.r_rms()),
            qip::overloads::operator-(Nuclear::formPotential(nuc2, rgrid.r()),
                                      Nuclear::formPotential(nuc1, rgrid.r()))),
        m_dr2(nuc2.r_rms() * nuc2.r_rms() - nuc1.r_rms() * nuc1.r_rms()) {}

  std::string name() const override final { return "dV_FS"; }
  std::string units() const override final { return "GHz/fm^2"; }
  double dr2() const { return m_dr2; }
};

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Ffs(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"", "Field shift operator"},
       {"drrms", {"Delta in r_rms (in fm) - must not be zero [0.01]"}},
       {"", "The following are normally not set:"},
       {"dt", {"Delta in skin-thickness t (in fm) [0.0]"}},
       {"dbeta", {"Delta in quadrupole deformation parameter, beta [0.0]"}},
       {"print", "Print details? [true]"}});
  if (input.has_option("help")) {
    return nullptr;
  }

  const auto nuc1 = wf.nucleus();

  const auto drrms = input.get("drrms", 0.01);
  const auto dt = input.get("dt", 0.0);
  const auto dbeta = input.get("dbeta", 0.0);

  auto nuc2 = nuc1;
  nuc2.set_rrms(nuc1.r_rms() + drrms);
  nuc2.t() = nuc1.t() + dt;
  nuc2.beta() = nuc1.beta() + dbeta;

  if (input.get("print", true)) {
    std::cout << "Field shift operator:\n";
    std::cout << "Nuc1: " << nuc1 << "\n";
    std::cout << "Nuc2: " << nuc2 << "\n";
    const auto dr2 = nuc2.r_rms() * nuc2.r_rms() - nuc1.r_rms() * nuc1.r_rms();
    std::cout << "d<r^2> = " << dr2 << " fm^2";
  }

  return std::make_unique<F_fs>(wf.grid(), nuc1, nuc2);
}

} // namespace DiracOperator
