#include "MBPT/Feynman.hpp"
#include "DiracODE/DiracODE.hpp"
#include "MBPT/FeynmanSigma.hpp"
#include "MBPT/Sigma2.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include <algorithm>
#include <cassert>
#include <complex>
#include <numeric>

//==============================================================================
TEST_CASE("MBPT: Feynman", "[MBPT][Feynman]") {

  Wavefunction wf({1000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear"},
                  {"Na", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Ne]");
  wf.formBasis({"35spdfg", 40, 7, 1.0e-3, 1.0e-4, 40.0, false});
  wf.solve_valence("4sp3d");

  const std::size_t stride = 4;
  const auto i0 = wf.grid().getIndex(1.0e-4);
  const auto size = (wf.grid().getIndex(30.0) - i0) / stride + 1;

  MBPT::Feynman Fy(wf.vHF(), i0, stride, size,
                   {MBPT::Screening::exclude, MBPT::HoleParticle::exclude, 4,
                    -0.33 * wf.energy_gap(), 0.05, 2.0},
                   1);

  std::cout << "Vex:\n";
  for (auto &v : wf.valence()) {

    const auto vx = Fy.get_Vx_kappa(v.kappa());

    const auto dv = vx * v;
    const auto dex = qip::inner_product(v.f(), dv.f()) / double(Fy.stride());
    const auto dex0 = v * (wf.vHF()->vexFa(v));

    std::cout << v << " " << v.en() << " " << dex0 << " " << dex << "\n";
  }

  std::cout << "\n";

  const auto [holes, excited] =
      DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel());
  Coulomb::YkTable yk(holes, excited);

  std::cout << "\nSigma(2):\n";
  for (auto &v : wf.valence()) {
    const auto Sd = Fy.Sigma_direct(v.kappa(), v.en());

    const auto de1 = v * (Sd * v);

    const auto de0 = MBPT::Sigma_vw(v, v, yk, holes, excited);

    std::cout << v << " " << v.en() << " " << de0 << " " << de1 << " "
              << "\n";
  }

  //----------------------------------------------------------------------------
  MBPT::Feynman Fy_sc(wf.vHF(), i0, stride, size,
                      {MBPT::Screening::include, MBPT::HoleParticle::exclude, 4,
                       -0.33 * wf.energy_gap(), 0.05, 2.0},
                      1);
  // MBPT::FeynmanSigma Fy2_sc(
  //     wf.vHF(), wf.basis(),
  //     MBPT::Sigma_params{MBPT::Method::Feynman, 1, false, 4, false, false,
  //                        -0.33 * wf.energy_gap(), 0.05, 2.0, true, false},
  //     {1.0e-4, 30.0, 4}, "");

  MBPT::Feynman Fy_scX(wf.vHF(), i0, stride, size,
                       {MBPT::Screening::only, MBPT::HoleParticle::exclude, 4,
                        -0.33 * wf.energy_gap(), 0.05, 2.0},
                       1);

  std::cout << "\nsc:\n";
  for (auto &v : wf.valence()) {
    const auto Sd = Fy.Sigma_direct(v.kappa(), v.en());
    const auto Sd_sc = Fy_sc.Sigma_direct(v.kappa(), v.en());

    const auto Sd_scX = Fy_scX.Sigma_direct(v.kappa(), v.en());

    const auto sc1 = v * (Sd_sc * v) - v * (Sd * v);

    const auto sc3 = v * (Sd_scX * v);

    std::cout << v << " " << sc1 << " "
              << " " << sc3 << "\n";
  }

  //----------------------------------------------------------------------------

  MBPT::Feynman Fy_hp(wf.vHF(), i0, stride, size,
                      {MBPT::Screening::exclude, MBPT::HoleParticle::include, 4,
                       -0.33 * wf.energy_gap(), 0.05, 2.0},
                      1);
  // MBPT::FeynmanSigma Fy2_hp(
  //     wf.vHF(), wf.basis(),
  //     MBPT::Sigma_params{MBPT::Method::Feynman, 1, false, 4, false, false,
  //                        -0.33 * wf.energy_gap(), 0.05, 2.0, screen, true},
  //     {1.0e-4, 30.0, 4}, "");

  std::cout << "\nhp:\n";
  for (auto &v : wf.valence()) {
    const auto Sd = Fy.Sigma_direct(v.kappa(), v.en());
    const auto Sd_hp = Fy_hp.Sigma_direct(v.kappa(), v.en());

    const auto hp1 = v * (Sd * v) - v * (Sd_hp * v);

    // const auto hp2 = v * (Fy2.FeynmanDirect(v.kappa(), v.en()).drj() * v) -
    //                  v * (Fy2_hp.FeynmanDirect(v.kappa(), v.en()).drj() * v);

    std::cout << v << " " << hp1 << " "
              << "\n";
  }

  //----------------------------------------------------------------------------

  MBPT::Feynman Fy_ao(wf.vHF(), i0, stride, size,
                      {MBPT::Screening::include, MBPT::HoleParticle::include, 4,
                       -0.33 * wf.energy_gap(), 0.05, 2.0},
                      1);
  // MBPT::FeynmanSigma Fy2_sc(
  //     wf.vHF(), wf.basis(),
  //     MBPT::Sigma_params{MBPT::Method::Feynman, 1, false, 4, false, false,
  //                        -0.33 * wf.energy_gap(), 0.05, 2.0, true, false},
  //     {1.0e-4, 30.0, 4}, "");

  MBPT::Feynman Fy_higher(wf.vHF(), i0, stride, size,
                          {MBPT::Screening::only, MBPT::HoleParticle::include,
                           4, -0.33 * wf.energy_gap(), 0.05, 2.0},
                          1);

  std::cout << "\nsc only:\n";
  for (auto &v : wf.valence()) {
    const auto Sd_2o = Fy.Sigma_direct(v.kappa(), v.en());
    const auto Sd_ao = Fy_ao.Sigma_direct(v.kappa(), v.en());

    const auto Sd_higher = Fy_higher.Sigma_direct(v.kappa(), v.en());

    const auto sc1 = v * (Sd_ao * v) - v * (Sd_2o * v);

    const auto sc3 = v * (Sd_higher * v);

    std::cout << v << " " << sc1 << " "
              << " " << sc3 << "\n";
  }
}
