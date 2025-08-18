#include "Modules/exampleModule.hpp"
#include "DiracOperator/DiracOperator.hpp" //For E1 operator
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/InputBlock.hpp"
#include "LinAlg/LinAlg.hpp"
#include "MBPT/Feynman.hpp"
#include "MBPT/RDMatrix.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <vector>

namespace Module {

// want this function to construct the VBr coordinate operator
std::vector<MBPT::GMatrix> form_VBr(const IO::InputBlock &input,
                                    const Wavefunction &wf) {
  // returns pointer to HF class that will have core and valence orbitals for the atom I give to the input file
  auto vHF = wf.vHF();

  const auto r0 = input.get({"Grid"}, "r0", 1.0e-6);
  const auto rmax = input.get({"Grid"}, "rmax", 120.0);

  const auto default_stride = [&]() {
    // By default, choose stride such that there is 150 points over [1e-4,30]
    const auto dstride =
        int(wf.grid().getIndex(30.0) - wf.grid().getIndex(1.0e-4)) / 150;
    return (dstride <= 2) ? 2 : dstride;
  }();
  const auto stride = input.get({"Correlations"}, "stride", default_stride);
  const auto i0 = vHF->grid().getIndex(r0);
  const auto size = (vHF->grid().getIndex(rmax) - i0) / stride + 1;

  auto max_ki = 2 * DiracSpinor::max_l(wf.core());

  // ---- CLASS ACTUALLY STARTS DOING STUFF HERE ---- //
  auto Vx_kappa = std::vector<MBPT::GMatrix>(
      std::size_t(max_ki + 1), {i0, stride, size, false, vHF->grid_sptr()});

  // initialises the hydrogen wave function object
  Wavefunction wfH(vHF->grid_sptr(), {"1", 0, "Ball"});

  wfH.set_HF();
  wfH.solve_core(false);
  wfH.formBasis(SplineBasis::Parameters(
      "90spdfghi", 90, 9, 1e-4, 0.0, 90.0,
      true)); // forms the hydrogen wavefunction basis; NOTE: the last boolean option is for including the negative energy states in the basis

  // constructs the Breit matrix as VBr(r1,r2)=\sum_i[VBrF_i](r1)F_i^†(r2)
  for (const auto &Fn : wfH.basis()) {
    auto VxFn = vHF->VBr(Fn);
    if (Fn.k_index() < max_ki + 1) {
      Vx_kappa[Fn.k_index()].add(
          VxFn, Fn,
          1.0); // multiplies VBrFa(r1) by Fa^†(r2) from the right
    }
  }

  // includes integration measures
  for (auto &Vx : Vx_kappa) {
    Vx.dri_in_place();
    Vx.drj_in_place();
  }
  return Vx_kappa;
}

/*
// want this function to calculate matrix elements of the VBr operator with respect to core and valence HF states for the atom i give into the input file
void exampleModule(const IO::InputBlock &input, const Wavefunction &wf) {

  // returns pointer to HF class that will have core and valence orbitals for the atom I give to the input file
  auto vHF = wf.vHF();

  const auto r0 = input.get({"Grid"}, "r0", 1.0e-6);
  const auto rmax = input.get({"Grid"}, "rmax", 120.0);

  const auto default_stride = [&]() {
    // By default, choose stride such that there is 150 points over [1e-4,30]
    const auto dstride =
        int(wf.grid().getIndex(30.0) - wf.grid().getIndex(1.0e-4)) / 150;
    return (dstride <= 2) ? 2 : dstride;
  }();
  const auto stride = input.get({"Correlations"}, "stride", default_stride);
  const auto i0 = vHF->grid().getIndex(r0);
  const auto size = (vHF->grid().getIndex(rmax) - i0) / stride + 1;

  auto max_ki = 2 * DiracSpinor::max_l(wf.core());

  // construct coordinate VBr matrix
  std::vector<MBPT::GMatrix> VBr = form_VBr(input, wf);

  std::cout << "State         Regular      Coordinate matrix" << std::endl;

  // loop over all core orbitals of the HF core for the atom I give to the input file, not the hydrogen core
  for (const auto &Fi : wf.valence()) {
    auto expval_gold = Fi * vHF->VBr(Fi);
    auto expval_feyn =
        Fi * (VBr.at(std::size_t(Angular::indexFromKappa(Fi.kappa()))) * Fi);

    std::cout << Fi.symbol() << "      " << expval_gold << "       "
              << expval_feyn << std::endl;
  }
}
*/

// // this is my attempt at writing the code for forming the BCHF basis using Derevianko's method of constructing the B-spline basis
// void exampleModule(const IO::InputBlock &input, const Wavefunction &wf) {

//   // extracting only the s orbital basis states
//   auto s_basis =
//       qip::select_if(wf.basis(), [](const auto &f) { return f.kappa() == -1; });

//   int N = s_basis.size(); // s-orbital basis size

//   // allows us to calculate [VBrFa](r) without BCHF basis from input file
//   HF::Breit Br(1.0);

//   std::cout << "Basis size is " << N << std::endl;
//   auto vHF = wf.vHF();  // pointer to HF
//   auto basis = s_basis; // vector of spline basis states

//   // extracts energy of highest core electron, Fermi energy
//   //auto Fermi_en = vHF->core()[vHF->core().size() - 1].en();

//   //std::cout << "Fermi energy is " << Fermi_en << std::endl;

//   // initialise matrix to store (dU)ik + ei where dU = U_BCHF - U_CHF = VBr
//   LinAlg::Matrix<double> dUik_plus_ei{N, N};

// #pragma omp parallel for collapse(2)
//   for (int i = 0; i < N; i++) {
//     for (int k = 0; k < N; k++) {

//       // extracting the ith and kth basis state
//       const auto &Fi = basis[i];
//       const auto &Fk = basis[k];

//       auto dU_ik = Fi * Br.VbrFa(Fk, wf.core()); // (VBr)_ik

//       dUik_plus_ei[i][k] = dU_ik;

//       if (i == k) {
//         dUik_plus_ei[i][k] += Fi.en();
//       }
//     }
//   }
//   // solves eigenvalue equation for BCHF eigenvalue varepsilonj
//   auto [e_values, e_vectors] = LinAlg::symmhEigensystem(dUik_plus_ei);

//   std::cout << std::endl;

//   for (int j = 0; j < s_basis.size(); j++) {
//     // only print energies of states with energies above Fermi level?
//     //if (e_values.size() < Fermi_en) {
//     std::cout << e_values[j] * PhysConst::Hartree_invcm << " "
//               << s_basis.at(j).en() * PhysConst::Hartree_invcm << std::endl;
//     //}
//   }
// }

// comparing my frequency-dependent Breit radial integrals against the full ``analytical" integrals without using the single for loop tricks
void exampleModule(const IO::InputBlock &input, const Wavefunction &wf) {

  const auto HFcore = wf.core();

  const auto Fa = HFcore[16];
  const auto Fb = HFcore[3];
  int k = 6;

  const auto num_points = Fa.grid().num_points();

  const auto maxi = std::min({Fa.max_pt(), Fb.max_pt(), num_points}); // ok?
  const auto irmax = (maxi == 0 || maxi > num_points) ? num_points : maxi;
  // nb bmax may be num_points
  const auto bmax = irmax;

  // Quadrature integration weights:
  const auto weights = [=](std::size_t i) {
    if (i < NumCalc::Nquad)
      return NumCalc::dq_inv * NumCalc::cq[i];
    if (i < num_points - NumCalc::Nquad)
      return 1.0;
    return NumCalc::dq_inv * NumCalc::cq[num_points - i - 1];
  };

  auto X = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) + Fa.g(i) * Fb.f(i));
  };

  auto Y = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) - Fa.g(i) * Fb.f(i));
  };

  // fills vector with spherical Bessel functions of first kind j_k(\alpha*\omega*r)
  const auto jkwr = SphericalBessel::fillBesselVec_kr(k, 1.0, Fa.grid().r());

  // fills vector with spherical Bessel functions of second kind y_k(\alpha*\omega*r)
  const auto ykwr =
      SphericalBessel::fillSecondBesselVec_kr(k, 1.0, Fa.grid().r());

  double w = pow(PhysConst::alpha2, 6);

  //----------- TESTING g INTEGRALS

  // calculates g0 and ginf using the single for loop trick
  std::vector<double> g0;
  std::vector<double> ginf;
  Coulomb::gk_ab_freqw(k, Fa, Fb, g0, ginf, maxi, 1.0);

  // initialises vector that will store the g^{k(ac)}_{bd}(r,w) integrals calculated using the full integral definition
  std::vector<double> gexact_0(num_points, 0.0);
  std::vector<double> gexact_inf(num_points, 0.0);

  // does the two separate integrals from 0 to r and then from r to infty in the naive double for loop method
  for (int i = 1; i < irmax; i++) {
    for (int j = 0; j <= i - 1; j++) {
      gexact_0[i] += ykwr[i] * jkwr[j] * X(j) * weights(j) * Fa.grid().drdu(j);
    }
    for (int j = i; j < irmax; j++) {
      gexact_inf[i] +=
          jkwr[i] * ykwr[j] * X(j) * weights(j) * Fa.grid().drdu(j);
    }
    gexact_0[i] = -w * (2 * k + 1.0) * gexact_0[i] * Fa.grid().du();
    gexact_inf[i] = -w * (2 * k + 1.0) * gexact_inf[i] * Fa.grid().du();
  }

  // std::cout << "Comparing the g integrals:" << std::endl;
  // for (int i = 0; i < (int)num_points / 4; i++) {
  //   std::cout << Fa.grid()(4 * i) << "   "
  //             << (g0[4 * i] - gexact_0[4 * i]) / gexact_0[4 * i]
  //             << "            "
  //             << (ginf[4 * i] - gexact_inf[4 * i]) / gexact_inf[4 * i]
  //             << std::endl;
  // }

  //-------TESTING h INTEGRALS

  // calculates g0 and ginf using the single for loop trick
  std::vector<double> h0;
  std::vector<double> hinf;
  Coulomb::hk_ab_freqw(k, Fa, Fb, h0, hinf, maxi, 1.0);

  // initialises vector that will store the g^{k(ac)}_{bd}(r,w) integrals calculated using the full integral definition
  std::vector<double> hexact_0(num_points, 0.0);
  std::vector<double> hexact_inf(num_points, 0.0);

  // does the two separate integrals from 0 to r and then from r to infty in the naive double for loop method
  for (int i = 1; i < irmax; i++) {
    for (int j = 0; j <= i - 1; j++) {
      hexact_0[i] += ykwr[i] * jkwr[j] * Y(j) * weights(j) * Fa.grid().drdu(j);
    }
    for (int j = i; j < irmax; j++) {
      hexact_inf[i] +=
          jkwr[i] * ykwr[j] * Y(j) * weights(j) * Fa.grid().drdu(j);
    }
    hexact_0[i] = -w * (2 * k + 1.0) * hexact_0[i] * Fa.grid().du();
    hexact_inf[i] = -w * (2 * k + 1.0) * hexact_inf[i] * Fa.grid().du();
  }

  // std::cout << "Comparing the h integrals:" << std::endl;
  // for (int i = 0; i < (int)num_points / 4; i++) {
  //   std::cout << Fa.grid()(4 * i) << "   "
  //             << (h0[4 * i] - hexact_0[4 * i]) / hexact_0[4 * i]
  //             << "            "
  //             << (hinf[4 * i] - hexact_inf[4 * i]) / hexact_inf[4 * i]
  //             << std::endl;
  // }

  //------- TESTING THE v^{k(ac)}_{bd} radial integrals -------//

  auto Xij = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) + Fa.g(i) * Fb.f(i));
  };

  auto Yij = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) - Fa.g(i) * Fb.f(i));
  };

  std::vector<double> Pkbd(num_points, 0.0);
  std::vector<double> Qkbd(num_points, 0.0);

  for (int i = 0; i < num_points; i++) {
    Pkbd[i] = ((Fa.kappa() - Fb.kappa()) / k) * Xij(i) - Yij(i);
    Qkbd[i] = ((Fa.kappa() - Fb.kappa()) / (k + 1.0)) * Xij(i) + Yij(i);
  }

  std::vector<double> v1;
  std::vector<double> v2;
  std::vector<double> v3;
  std::vector<double> v4;

  w = PhysConst::alpha2;

  Coulomb::vk_ab_freqw(k, Fa, Fb, Fa.grid(), v1, v2, v3, v4, maxi, w);

  std::vector<double> v1_exact(num_points, 0.0);
  std::vector<double> v2_exact(num_points, 0.0);
  std::vector<double> v3_exact(num_points, 0.0);
  std::vector<double> v4_exact(num_points, 0.0);

  const auto &r = Fa.grid().r();
  // fills vector with spherical Bessel functions of first kind j_k(\alpha*\omega*r)
  const auto jkpluswr = SphericalBessel::fillBesselVec_kr(k + 1, w, r);
  const auto jkminuswr = SphericalBessel::fillBesselVec_kr(k - 1, w, r);

  // fills vector with spherical Bessel functions of second kind y_k(\alpha*\omega*r)
  const auto ykpluswr = SphericalBessel::fillSecondBesselVec_kr(k + 1, w, r);
  const auto ykminuswr = SphericalBessel::fillSecondBesselVec_kr(k - 1, w, r);

  double A1 = 0.0, A2 = 0.0, A3 = 0.0;

  for (int i = 1; i < irmax; i++) {
    for (int j = 0; j <= i - 1; j++) {
      // v1_exact[i] += -2.0 * w * ykpluswr[i] * jkminuswr[j] * Pkbd[j] *
      //                    weights(j) * Fa.grid().drdu(j) * Fa.grid().du() -
      //                2.0 * ((2 * k + 1.0) / (w * w)) *
      //                    (pow(r[j], k - 1) / pow(r[i], k + 2)) * Pkbd[j] *
      //                    weights(j) * Fa.grid().drdu(j) * Fa.grid().du();

      //=================== for testing purposes just split up into first and second terms

      v1_exact[i] += -2.0 * w * ykpluswr[i] * jkminuswr[j] * Pkbd[j] *
                     weights(j) * Fa.grid().drdu(j) * Fa.grid().du();
      v1_exact[i] += -2.0 * ((2 * k + 1.0) / (w * w)) *
                     (pow(r[j], k - 1) / pow(r[i], k + 2)) * Pkbd[j] *
                     weights(j) * Fa.grid().drdu(j) * Fa.grid().du();

      v2_exact[i] +=
          ykminuswr[i] * jkpluswr[j] * Qkbd[j] * weights(j) * Fa.grid().drdu(j);
    }
    // v1_exact[i] = 2.0 * w * v1_exact[i] * Fa.grid().du();
    // v1_exact[i] = ((2.0 * (2 * k + 1.0)) / w * w) * v1_exact[i] *
    //               Fa.grid().du() * (1.0 / pow(r[i], k + 1));
    //v1_exact[i] = -2.0 * v1_exact[i] * Fa.grid().du();
    v2_exact[i] = -2.0 * w * v2_exact[i] * Fa.grid().du();
  }

  // std::cout << "Comparing the v integrals that go up to r1:" << std::endl;
  // for (int i = 0; i < (int)irmax / 4; i++) {
  //   std::cout << Fa.grid()(4 * i) << "   "
  //             << (v1_exact[4 * i] - v1[4 * i]) / v1_exact[4 * i] << "          "
  //             << (v2_exact[4 * i] - v2[4 * i]) / v2_exact[4 * i] << std::endl;
  // }

  double B1 = 0.0, B2 = 0.0, B3 = 0.0;

  for (int i = 1; i < irmax; i++) {
    for (int j = i; j <= irmax; j++) {
      // first term
      v3_exact[i] += -2.0 * w * ykpluswr[j] * jkminuswr[i] * Qkbd[j] *
                     weights(j) * Fa.grid().drdu(j) * Fa.grid().du();

      //second term
      v3_exact[i] += -2.0 * ((2 * k + 1.0) / (w * w)) *
                     (pow(r[i], k - 1) / pow(r[j], k + 2)) * Qkbd[j] *
                     weights(j) * Fa.grid().drdu(j) * Fa.grid().du();

      v4_exact[i] +=
          ykminuswr[j] * jkpluswr[i] * Pkbd[j] * weights(j) * Fa.grid().drdu(j);
    }
    // v1_exact[i] = 2.0 * w * v1_exact[i] * Fa.grid().du();
    // v1_exact[i] = ((2.0 * (2 * k + 1.0)) / w * w) * v1_exact[i] *
    //               Fa.grid().du() * (1.0 / pow(r[i], k + 1));
    //v1_exact[i] = -2.0 * v1_exact[i] * Fa.grid().du();
    v4_exact[i] = -2.0 * w * v4_exact[i] * Fa.grid().du();
  }

  // std::cout << "Comparing the v integrals that go up to infinity:" << std::endl;
  // for (int i = 0; i < (int)irmax / 4; i++) {
  //   std::cout << Fa.grid()(4 * i) << "   "
  //             << (v3_exact[4 * i] - v3[4 * i]) / v3_exact[4 * i] << "          "
  //             << (v4_exact[4 * i] - v4[4 * i]) / v4_exact[4 * i] << std::endl;
  // }

  // std::vector<double> energies;
  // for (const auto &Fb : HFcore) {
  //   energies.push_back(Fb.en());
  //   std::cout << Fb.en() << std::endl;
  // }

  //---------------- TESTING ALTERNATE DEFINITION OF BESSEL FUNCTIONS
  // There is a problem with the code in that if k=0, we will have to evaluate Bessel functions of order -1, which doesn't work with the default GSL functions for the Bessel functions

  // testing whether my two functions of the first and second kind Bessel functions give the same results, which they do

  std::vector<double> jkwr_test = SphericalBessel::fillBesselVec_kr(2, 2.0, r);
  std::vector<double> jkwr_test_alt =
      SphericalBessel::fillBesselVec_kr_alt(2, 2.0, r);

  std::vector<double> ykwr_test =
      SphericalBessel::fillSecondBesselVec_kr(2, 2.0, r);
  std::vector<double> ykwr_test_alt =
      SphericalBessel::fillSecondBesselVec_kr_alt(2, 2.0, r);

  // difference between the second kind Bessel functions seem large for small values but the spherical Bessel functions of second kind are huge at small r so I think this is fine
  // for (int i = 0; i < num_points / 4; i++) {
  //   std::cout << jkwr_test[i] - jkwr_test_alt[i] << "            "
  //             << SphericalBessel::exactGSL_YL(2, 2.0 * r[i]) << std::endl;
  // }

  // now test that my function can calculate spherical Bessel functions of negative order
  // for (int i = 0; i < (int)num_points / 4; i++) {
  //   std::cout << r[4 * i] << "          "
  //             << SphericalBessel::exactGSL_JL_alt(-2, 1000.0 * w * r[4 * i])
  //             << "          "
  //             << SphericalBessel::exactGSL_YL_alt(-2, 1000.0 * w * r[4 * i])
  //             << std::endl;
  // }

  //--------- TEST TO SEE IF THE BESSEL FUNCTION THINGS REDUCE TO THE RIGHT THINGS IN THE APPROPRIATE LIMIT

  // will store the parts of the frequency-independent integral I want to compare against
  std::vector<double> g1;
  std::vector<double> g2;
  std::vector<double> b1;
  std::vector<double> b2;
  std::vector<double> vec;
  w = PhysConst::alpha2;

  Coulomb::gk_ab(k - 1, Fa, Fb, g1, vec, maxi);
  Coulomb::gk_ab(k + 1, Fa, Fb, g2, vec, maxi);
  Coulomb::bk_ab(k - 1, Fa, Fb, b1, vec, maxi);
  Coulomb::bk_ab(k + 1, Fa, Fb, b2, vec, maxi);

  std::vector<double> v0_freq_ind(num_points, 0.0);
  for (int i = 0; i < num_points; i++) {
    v0_freq_ind[i] =
        ((Fa.kappa() - Fb.kappa()) / k) * (g1[i] - g2[i]) - b1[i] + b2[i];
  }

  Coulomb::vk_ab_freqw(k, Fa, Fb, Fa.grid(), v1, v2, v3, v4, maxi, w);

  // std::cout << "Comparing the frequency-ind. and frequency-dep. v "
  //              "integrals that go up to r1 in limit of w->0:"
  //           << std::endl;
  // for (int i = 0; i < (int)irmax / 4; i++) {
  //   std::cout << Fa.grid()(4 * i) << "   " << v1[4 * i] << "   "
  //             << v1_exact[4 * i] << "          " << v0_freq_ind[4 * i]
  //             << std::endl;
  // }

  //=================== TESTING IF LIMIT EVEN WORKS FOR A SINGLE POINT
  // testing v integrands at single point

  const auto r1 = 48.45643;
  const auto r2 = 14.897;
  k = 6;
  w = pow(PhysConst::alpha2, 1);

  // std::cout << "Freq-dep. equation in w->0 limit: "
  //           << -2.0 * (w * SphericalBessel::exactGSL_JL_alt(k - 1, w * r2) *
  //                          SphericalBessel::exactGSL_YL_alt(k + 1, w * r1) +
  //                      ((double(2 * k + 1)) / (w * w)) *
  //                          (pow(r2, k - 1) / pow(r1, k + 2)))
  //           << std::endl;
  // std::cout << "Freq-ind. equation: "
  //           << (pow(r2, k - 1) / pow(r1, k)) -
  //                  ((pow(r2, k + 1) / pow(r1, k + 2)))
  //           << std::endl;

  // testing the r^k/r^k+1 -> -w(2k + 1)jk(wr)yk(wr) equivalence in low-w limit
  // std::cout << "Freq-dep. equation in w->0 limit: "
  //           << -w * (2 * 6 + 1.0) *
  //                  SphericalBessel::exactGSL_JL_alt(6, w * 14.897) *
  //                  SphericalBessel::exactGSL_YL_alt(6, w * 48.45643)
  //           << std::endl;
  // std::cout << "Freq-ind. equation: " << (pow(14.897, 6) / pow(48.45643, 6 + 1))
  //           << std::endl;

  // testing term that should go to zero
  // std::cout
  //     << w * SphericalBessel::exactGSL_JL_alt(6 + 1, 0.000001 * w * 1000.0) *
  //            SphericalBessel::exactGSL_YL_alt(6 - 1, 0.000001 * w * 0.000045643)
  //     << std::endl;

  //======================== COMPARING MATRIX ELEMENTS IN w->0 LIMIT

  DiracSpinor Fi = HFcore[12];
  DiracSpinor Fj = HFcore[3];
  DiracSpinor Fk = HFcore[3];
  DiracSpinor Fl = HFcore[12];
  k = 3;
  const auto br = HF::Breit(1.0);

  std::cout << "Frequency-independent Breit integral: "
            << br.Bk_abcd(k, Fi, Fj, Fk, Fl) << std::endl;
  std::cout << "Frequency-dependent Breit integral: "
            << br.Bk_abcd_eac_freqw(k, Fi, Fj, Fk, Fl) << std::endl;
}

// =================================================================================================================== //
/*
void exampleModule(const IO::InputBlock &input, const Wavefunction &wf) {
  // This is an example module, designed to help you write new modules

  // In this example, we will solve a new wavefunction, assuming a different
  // nuclear charge distribution, and see the difference in the energies and E1
  // matrix elements this produces.

  // Read in some optional input options: A, rms, and type
  // "checkBlock" is optional, but helps prevent errors on input options:
  input.check({{"A", "Atomic mass number [0]"},
               {"rrms", "root-mean-square nuclear radii [-1]"},
               {"type", "Fermi/ball/pointlike [Fermi]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  // If no A is specified, use 0 (i.e., poinlike)
  auto A = input.get("A", 0);
  // if rrms<0, the default value (for the given A) will be looked up
  auto rrms = input.get("rrms", -1.0);
  // Will assume a Fermi nucleus, unless A or r =0
  auto nuc_type = input.get<std::string>("type", "Fermi");
  // Set nuc. type explicitely to 'pointlike' if A=0, or r_rms = 0.0
  if (A == 0 || rrms == 0.0)
    nuc_type = "pointlike";
  if (nuc_type == "pointlike") {
    A = 0;
    rrms = 0.0;
  }

  // Use the same Grid and alpha, but different nuclear parameters (except Z)
  Wavefunction wfpt(wf.grid().params(),
                    {wf.Znuc(), A, nuc_type, rrms, Nuclear::default_t},
                    wf.alpha() / PhysConst::alpha);

  std::cout << "\n";
  IO::print_line();
  std::cout << "Calculating finite-nuclear size corrections for \n"
            << wf.atom() << ", " << wf.nucleus() << "\nvs.\n"
            << wfpt.atom() << ", " << wfpt.nucleus() << "\n\n";

  // Note: Assume only Hartree-Fock approximation, no Breit, QED, Sigma etc.

  // Solve Hartree-Fock core for new wavefuinction:
  wfpt.solve_core("HartreeFock", 0.0, wf.coreConfiguration());
  // print the new core energies:
  wfpt.printCore();

  // Solve Hartree-Fock valence
  wfpt.solve_valence(DiracSpinor::state_config(wf.valence()));
  wfpt.printValence();

  // Calculate the energy shifts in atomic units, and GHz
  std::cout << "\nFinite nuclear charge energy shifts:\n";
  std::cout << "  state            au           GHz\n";
  for (auto i = 0ul; i < wf.valence().size(); ++i) {
    const auto del_e = wf.valence()[i].en() - wfpt.valence()[i].en();
    printf("%7s  %12.5e  %12.5e\n", wf.valence()[i].symbol().c_str(), del_e,
           del_e * PhysConst::Hartree_GHz);
  }

  // Now, we calculate E1 matrix elements.
  std::cout << "\nFinite nuclear charge shifts to E1 reduced MEs (no RPA):\n";

  std::cout << "                      A=" << wf.Anuc()
            << "         A=" << wfpt.Anuc() << "         Shift\n";

  // 1) Create E1 operator:
  const auto hE1 = DiracOperator::E1(wf.grid());

  // 2) Loop through each pair of valence states, calc E1 matrix elements:
  for (auto a = 0ul; a < wf.valence().size(); ++a) {
    const auto &Fa = wf.valence()[a];
    const auto &F0a = wfpt.valence()[a]; // pointlike orbital
    for (auto b = 0ul; b < a; ++b) {
      const auto &Fb = wf.valence()[b];
      const auto &F0b = wfpt.valence()[b]; // pointlike orbital

      // Skip the MEs which are zero due to selection rules:
      if (hE1.isZero(Fa.kappa(), Fb.kappa()))
        continue;

      const auto e1 = hE1.reducedME(Fa, Fb);
      const auto e10 = hE1.reducedME(F0a, F0b);
      const auto delta = e1 - e10;

      std::cout << Fa << " - " << Fb << ": ";
      printf("%12.5e  %12.5e  %12.5e\n", e1, e10, delta);
    }
  }
}
*/

} // namespace Module
