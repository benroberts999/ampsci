#include "Modules/DCPdiagrams.hpp"
#include "Angular/CkTable.hpp"
#include "Angular/Wigner369j.hpp"
#include "Angular/include.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/meTable.hpp"
#include "DiracOperator/include.hpp" //For E1 operator
#include "ExternalField/TDHF.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"
#include <optional>
#include <string>
#include <vector>

//Module: Run through runModules.cpp

//constexpr int jindex_kappa(int ka) { return (ka > 0) ? ka - 1 : -ka - 1; }
//! returns 2j given kappa
//constexpr int twoj_k(int ka) { return (ka > 0) ? 2 * ka - 1 : -2 * ka - 1; }
namespace Module {

//Functions computing each of 12 pairs of diagrams contributing to the induced amplitude at second order in perturbation theory

//Inputs:
// Rank of irreducible tensor formed from coupling of two irreducible tensor operators: k_1,k_2 -> K
// Two operators for external fields with ranks k_1 and k_2 respectively
// Core and excited basis states for atomic structure calculations

// Q^k table?-- later on first just call functions

//valence states v and w corresponding to the relevant transition

// Including core polarisation/RPA (Modified Hartree Fock potential dV) as an optional parameter

// Compute all matrix elements between core and excited states. First element of vector is table for operator 1 (time dependent)

// Output: table of matrix elements for given pair of operators (optional input is TDHF object for inclusion of RPA)
std::vector<Coulomb::meTable<double>>
compute_me(const DiracOperator::TensorOperator *const h1,
           const DiracOperator::TensorOperator *const h2,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited, const double omega,
           const ExternalField::CorePolarisation *dV1,
           const ExternalField::CorePolarisation *dV2) {
  //create tables
  Coulomb::meTable<double> matrix_elements1;
  Coulomb::meTable<double> matrix_elements2;

  // initialise vector
  std::vector<Coulomb::meTable<double>> ME_tables;

  //in general, for the evaluation of DCP diagrams, core-core, core-valence and valence-valence matrix elements are required

  //core core

  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    for (auto ib = 0ul; ib < core.size(); ++ib) {
      const auto &b = core[ib];
      //compute matrix elements with RPA/core polarisation corrections
      const auto d1_ab = h1->reducedME(a, b) + (dV1 ? dV1->dV(a, b) : 0.0);
      const auto d2_ab = h2->reducedME(a, b) + (dV2 ? dV2->dV(a, b) : 0.0);

      matrix_elements1.add(a, b, d1_ab);
      matrix_elements1.add(b, a, d1_ab * h1->symm_sign(a, b));

      matrix_elements2.add(a, b, d2_ab);
      matrix_elements2.add(b, a, d2_ab * h2->symm_sign(a, b));
    }
  }

  //valence-valence
  for (auto im = 0ul; im < excited.size(); ++im) {
    const auto &m = excited[im];
    for (auto in = 0ul; in < excited.size(); ++in) {
      const auto &n = excited[in];
      //compute matrix elements with RPA/core polarisation corrections
      const auto d1_mn = h1->reducedME(m, n) + (dV1 ? dV1->dV(m, n) : 0.0);
      const auto d2_mn = h2->reducedME(m, n) + (dV2 ? dV2->dV(m, n) : 0.0);

      matrix_elements1.add(m, n, d1_mn);
      matrix_elements1.add(n, m, d1_mn * h1->symm_sign(m, n));

      matrix_elements2.add(m, n, d2_mn);
      matrix_elements2.add(n, m, d2_mn * h2->symm_sign(m, n));
    }
  }

  //core-excited
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    for (auto im = 0ul; im < excited.size(); ++im) {
      const auto &m = excited[im];
      //compute matrix elements with RPA/core polarisation corrections
      const auto d1_am = h1->reducedME(a, m) + (dV1 ? dV1->dV(a, m) : 0.0);
      const auto d2_am = h2->reducedME(a, m) + (dV2 ? dV2->dV(a, m) : 0.0);

      matrix_elements1.add(a, m, d1_am);
      matrix_elements1.add(m, a, d1_am * h1->symm_sign(a, m));

      matrix_elements2.add(a, m, d2_am);
      matrix_elements2.add(m, a, d2_am * h2->symm_sign(a, m));
    }
  }

  ME_tables.push_back(matrix_elements1);
  ME_tables.push_back(matrix_elements2);

  return ME_tables;
}

// Input

// First pair of diagrams
double A1(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega) {
  //initialise A1
  double A1 = 0.0;

  //get tables
  Coulomb::meTable<double> met_1 = ME_tables.at(0);
  Coulomb::meTable<double> met_2 = ME_tables.at(1);

  //sum over core states, then excited states
  for (auto ia = 0ul; ia < core.size(); ++ia) {

    for (auto im = 0ul; im < excited.size(); ++im) {

      for (auto in = 0ul; in < excited.size(); ++in) {
        // std::cout<<__LINE__<<" "<<ia<<"\n";
        //get states
        const auto &a = core[ia];
        const auto &m = excited[im];
        const auto &n = excited[in];
        //compute energy denominator
        double e_denom =
            1 / (a.en() - m.en() + omega) * 1 / (a.en() - n.en() - omega);
        //matrix elements (Can include RPA)

        const auto dv_ma = met_1.getv(m, a);
        const auto dv_an = met_2.getv(a, n);

        //phase
        int kw = w.kappa();
        int kn = n.kappa();
        auto tjw = Angular::twoj_k(kw);
        auto tjn = Angular::twoj_k(kn);
        auto phase = Angular::neg1pow_2(tjw + tjn);

        //6j symbol
        int km = m.kappa();
        auto tjm = Angular::twoj_k(km);
        int ka = a.kappa();
        auto tja = Angular::twoj_k(ka);
        double sixj = Angular::sixj_2(2 * K, 2 * k1, 2 * k2, tja, tjn, tjm);

        //W^k
        double Wk = Coulomb::Wk_abcd(K, w, n, v, m);

        //contribution to sum
        auto cont = e_denom * dv_ma * dv_an * sixj * Wk * phase;
        A1 += cont;
      }
    }
  }
  A1 = (1 / (sqrt(2 * K + 1))) * pow(-1, k1 + k2) * A1;
  return A1;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Second pair of diagrams
double A2(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega) {
  //initialise A1
  double A2 = 0.0;

  //get tables
  Coulomb::meTable<double> met_1 = ME_tables.at(0);
  Coulomb::meTable<double> met_2 = ME_tables.at(1);

  //sum over core states, then excited states
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    //get states
    const auto &a = core[ia];
    for (auto ib = 0ul; ib < core.size(); ++ib) {
      const auto &b = core[ib];
      for (auto im = 0ul; im < excited.size(); ++im) {
        // std::cout<<__LINE__<<" "<<ia<<"\n";
        //get states
        const auto &m = excited[im];
        //compute energy denominator
        double e_denom =
            1 / (a.en() - m.en() + omega) * 1 / (b.en() - m.en() - omega);
        //matrix elements (Can include RPA)

        const auto dv_ma = met_1.getv(m, a);
        const auto dv_bm = met_2.getv(b, m);

        //phase
        int kw = w.kappa();
        int ka = a.kappa();
        auto tjw = Angular::twoj_k(kw);
        auto tja = Angular::twoj_k(ka);
        auto phase = Angular::neg1pow_2(tjw + tja);

        //6j symbol
        int km = m.kappa();
        auto tjm = Angular::twoj_k(km);
        int kb = b.kappa();
        auto tjb = Angular::twoj_k(kb);
        double sixj = Angular::sixj_2(2 * K, 2 * k1, 2 * k2, tjm, tjb, tja);

        //W^k
        double Wk = Coulomb::Wk_abcd(K, w, a, v, b);

        //contribution to sum
        auto cont = e_denom * dv_ma * dv_bm * sixj * Wk * phase;
        A2 += cont;
      }
    }
  }
  A2 = -1 * (1 / (sqrt(2 * K + 1))) * pow(-1, K) * A2;
  return A2;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Third pair of diagrams
double A3(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega) {
  //initialise A3
  double A3 = 0.0;

  //get tables
  Coulomb::meTable<double> met_1 = ME_tables.at(0);
  Coulomb::meTable<double> met_2 = ME_tables.at(1);

  //sum over core states, then excited states
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    //get states
    const auto &a = core[ia];
    for (auto ib = 0ul; ib < core.size(); ++ib) {
      const auto &b = core[ib];
      for (auto im = 0ul; im < excited.size(); ++im) {
        // std::cout<<__LINE__<<" "<<ia<<"\n";
        //get states
        const auto &m = excited[im];
        //compute energy denominator
        double e_denom = 1 / (a.en() - m.en() + omega) * 1 / (b.en() - m.en());
        //matrix elements (Can include RPA)

        const auto dv_ma = met_1.getv(m, a);
        const auto dv_ab = met_2.getv(a, b);

        //phase
        int kw = w.kappa();
        int kb = b.kappa();
        auto tjw = Angular::twoj_k(kw);
        auto tjb = Angular::twoj_k(kb);
        auto phase = Angular::neg1pow_2(tjw + tjb);

        //6j symbol
        int km = m.kappa();
        auto tjm = Angular::twoj_k(km);
        int ka = a.kappa();
        auto tja = Angular::twoj_k(ka);
        double sixj = Angular::sixj_2(2 * K, 2 * k1, 2 * k2, tja, tjb, tjm);

        //W^k
        double Wk = Coulomb::Wk_abcd(K, w, b, v, m);

        //contribution to sum
        auto cont = e_denom * dv_ma * dv_ab * sixj * Wk * phase;
        A3 += cont;
      }
    }
  }
  A3 = -1 * (1 / (sqrt(2 * K + 1))) * pow(-1, k1 + k2) * A3;
  return A3;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Fourth pair of diagrams
double A4(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega) {
  //initialise A4
  double A4 = 0.0;

  //get tables
  Coulomb::meTable<double> met_1 = ME_tables.at(0);
  Coulomb::meTable<double> met_2 = ME_tables.at(1);

  //sum over core states, then excited states
  for (auto im = 0ul; im < excited.size(); ++im) {
    //get states
    const auto &m = excited[im];
    for (auto in = 0ul; in < excited.size(); ++in) {
      const auto &n = excited[in];
      for (auto ia = 0ul; ia < core.size(); ++ia) {
        // std::cout<<__LINE__<<" "<<ia<<"\n";
        //get states
        const auto &a = core[ia];
        //compute energy denominator
        double e_denom = 1 / (a.en() - m.en() + omega) * 1 / (a.en() - n.en());
        //matrix elements (Can include RPA)

        const auto dv_ma = met_1.getv(m, a);
        const auto dv_nm = met_2.getv(n, m);

        //phase
        int kw = w.kappa();
        int ka = a.kappa();
        auto tjw = Angular::twoj_k(kw);
        auto tja = Angular::twoj_k(ka);
        auto phase = Angular::neg1pow_2(tjw + tja);

        //6j symbol
        int km = m.kappa();
        auto tjm = Angular::twoj_k(km);
        int kn = n.kappa();
        auto tjn = Angular::twoj_k(kn);
        double sixj = Angular::sixj_2(2 * K, 2 * k1, 2 * k2, tjm, tjn, tja);

        //W^k
        double Wk = Coulomb::Wk_abcd(K, w, a, v, n);

        //contribution to sum
        auto cont = e_denom * dv_ma * dv_nm * sixj * Wk * phase;
        A4 += cont;
      }
    }
  }
  A4 = (1 / (sqrt(2 * K + 1))) * pow(-1, K) * A4;
  return A4;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Fifth pair of diagrams
double A5(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega) {
  //initialise A5
  double A5 = 0.0;

  //get tables
  Coulomb::meTable<double> met_1 = ME_tables.at(0);
  Coulomb::meTable<double> met_2 = ME_tables.at(1);

  //sum over core states, then excited states
  for (auto im = 0ul; im < excited.size(); ++im) {
    //get states
    const auto &m = excited[im];
    for (auto ib = 0ul; ib < core.size(); ++ib) {
      const auto &b = core[ib];
      for (auto ia = 0ul; ia < core.size(); ++ia) {
        // std::cout<<__LINE__<<" "<<ia<<"\n";
        //get states
        const auto &a = core[ia];
        //compute energy denominator
        double e_denom = 1 / (a.en() - m.en() - omega) * 1 / (b.en() - m.en());
        //matrix elements (Can include RPA)

        const auto dv_ab = met_1.getv(a, b);
        const auto dv_bm = met_2.getv(b, m);

        //phase
        int kw = w.kappa();
        int km = m.kappa();

        auto tjw = Angular::twoj_k(kw);
        auto tjm = Angular::twoj_k(km);
        auto phase = Angular::neg1pow_2(tjw + tjm);

        //6j symbol
        int ka = a.kappa();
        auto tja = Angular::twoj_k(ka);
        int kb = b.kappa();
        auto tjb = Angular::twoj_k(kb);
        double sixj = Angular::sixj_2(2 * K, 2 * k1, 2 * k2, tjb, tjm, tja);

        //W^k
        double Wk = Coulomb::Wk_abcd(K, w, m, v, a);

        //contribution to sum
        auto cont = e_denom * dv_ab * dv_bm * sixj * Wk * phase;
        A5 += cont;
      }
    }
  }
  A5 = -1 * (1 / (sqrt(2 * K + 1))) * pow(-1, k1 + k2) * A5;
  return A5;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Sixth pair of diagrams
double A6(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega) {
  //initialise A6
  double A6 = 0.0;

  //get tables
  Coulomb::meTable<double> met_1 = ME_tables.at(0);
  Coulomb::meTable<double> met_2 = ME_tables.at(1);

  //sum over core states, then excited states
  for (auto im = 0ul; im < excited.size(); ++im) {
    //get states
    const auto &m = excited[im];
    for (auto in = 0ul; in < excited.size(); ++in) {
      const auto &n = excited[in];
      for (auto ia = 0ul; ia < core.size(); ++ia) {
        // std::cout<<__LINE__<<" "<<ia<<"\n";
        //get states
        const auto &a = core[ia];
        //compute energy denominator
        double e_denom = 1 / (a.en() - m.en() - omega) * 1 / (a.en() - n.en());
        //matrix elements (Can include RPA)

        const auto dv_nm = met_1.getv(n, m);
        const auto dv_an = met_2.getv(a, n);

        //phase
        int kw = w.kappa();
        int km = m.kappa();
        auto tjw = Angular::twoj_k(kw);
        auto tjm = Angular::twoj_k(km);
        auto phase = Angular::neg1pow_2(tjw + tjm);

        //6j symbol
        int ka = a.kappa();
        auto tja = Angular::twoj_k(ka);
        int kn = n.kappa();
        auto tjn = Angular::twoj_k(kn);
        double sixj = Angular::sixj_2(2 * K, 2 * k1, 2 * k2, tjn, tja, tjm);

        //W^k
        double Wk = Coulomb::Wk_abcd(K, w, m, v, a);

        //contribution to sum
        auto cont = e_denom * dv_nm * dv_an * sixj * Wk * phase;
        A6 += cont;
      }
    }
  }
  A6 = (1 / (sqrt(2 * K + 1))) * pow(-1, K) * A6;
  return A6;
}

// seventh pair of diagrams
double A7(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega) {
  //initialise A7
  double A7 = 0.0;

  //get tables
  Coulomb::meTable<double> met_1 = ME_tables.at(0);
  Coulomb::meTable<double> met_2 = ME_tables.at(1);

  //sum over core states, then excited states
  for (auto ia = 0ul; ia < core.size(); ++ia) {

    for (auto im = 0ul; im < excited.size(); ++im) {

      for (auto in = 0ul; in < excited.size(); ++in) {
        // std::cout<<__LINE__<<" "<<ia<<"\n";
        //get states
        const auto &a = core[ia];
        const auto &m = excited[im];
        const auto &n = excited[in];
        //compute energy denominator
        double e_denom = 1 / (a.en() - m.en()) * 1 / (a.en() - n.en() - omega);
        //matrix elements (Can include RPA)

        const auto dv_ma = met_2.getv(m, a);
        const auto dv_an = met_1.getv(a, n);

        //phase
        int kw = w.kappa();
        int kn = n.kappa();
        auto tjw = Angular::twoj_k(kw);
        auto tjn = Angular::twoj_k(kn);
        auto phase = Angular::neg1pow_2(tjw + tjn);

        //6j symbol
        int km = m.kappa();
        auto tjm = Angular::twoj_k(km);
        int ka = a.kappa();
        auto tja = Angular::twoj_k(ka);
        double sixj = Angular::sixj_2(2 * K, 2 * k1, 2 * k2, tja, tjm, tjn);

        //W^k
        double Wk = Coulomb::Wk_abcd(K, w, n, v, m);

        //contribution to sum
        auto cont = e_denom * dv_ma * dv_an * sixj * Wk * phase;
        A7 += cont;
      }
    }
  }
  A7 = (1 / (sqrt(2 * K + 1))) * pow(-1, K) * A7;
  return A7;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Eighth pair of diagrams
double A8(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega) {
  //initialise A1
  double A8 = 0.0;

  //get tables
  Coulomb::meTable<double> met_1 = ME_tables.at(0);
  Coulomb::meTable<double> met_2 = ME_tables.at(1);

  //sum over core states, then excited states
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    //get states
    const auto &a = core[ia];
    for (auto ib = 0ul; ib < core.size(); ++ib) {
      const auto &b = core[ib];
      for (auto im = 0ul; im < excited.size(); ++im) {
        // std::cout<<__LINE__<<" "<<ia<<"\n";
        //get states
        const auto &m = excited[im];
        //compute energy denominator
        double e_denom = 1 / (a.en() - m.en()) * 1 / (b.en() - m.en() - omega);
        //matrix elements (Can include RPA)

        const auto dv_ma = met_2.getv(m, a);
        const auto dv_bm = met_1.getv(b, m);

        //phase
        int kw = w.kappa();
        int ka = a.kappa();
        auto tjw = Angular::twoj_k(kw);
        auto tja = Angular::twoj_k(ka);
        auto phase = Angular::neg1pow_2(tjw + tja);

        //6j symbol
        int km = m.kappa();
        auto tjm = Angular::twoj_k(km);
        int kb = b.kappa();
        auto tjb = Angular::twoj_k(kb);
        double sixj = Angular::sixj_2(2 * K, 2 * k1, 2 * k2, tjm, tja, tjb);

        //W^k
        double Wk = Coulomb::Wk_abcd(K, w, a, v, b);

        //contribution to sum
        auto cont = e_denom * dv_ma * dv_bm * sixj * Wk * phase;
        A8 += cont;
      }
    }
  }
  A8 = -1 * (1 / (sqrt(2 * K + 1))) * pow(-1, k1 + k2) * A8;
  return A8;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Ninth pair of diagrams
double A9(const int K, const int k1, const int k2,
          std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
          const DiracSpinor &v, double omega) {
  //initialise A3
  double A9 = 0.0;

  //get tables
  Coulomb::meTable<double> met_1 = ME_tables.at(0);
  Coulomb::meTable<double> met_2 = ME_tables.at(1);

  //sum over core states, then excited states
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    //get states
    const auto &a = core[ia];
    for (auto ib = 0ul; ib < core.size(); ++ib) {
      const auto &b = core[ib];
      for (auto im = 0ul; im < excited.size(); ++im) {
        // std::cout<<__LINE__<<" "<<ia<<"\n";
        //get states
        const auto &m = excited[im];
        //compute energy denominator
        double e_denom = 1 / (a.en() - m.en()) * 1 / (b.en() - m.en() + omega);
        //matrix elements (Can include RPA)

        const auto dv_ma = met_2.getv(m, a);
        const auto dv_ab = met_1.getv(a, b);

        //phase
        int kw = w.kappa();
        int kb = b.kappa();
        auto tjw = Angular::twoj_k(kw);
        auto tjb = Angular::twoj_k(kb);
        auto phase = Angular::neg1pow_2(tjw + tjb);

        //6j symbol
        int km = m.kappa();
        auto tjm = Angular::twoj_k(km);
        int ka = a.kappa();
        auto tja = Angular::twoj_k(ka);
        double sixj = Angular::sixj_2(2 * K, 2 * k1, 2 * k2, tja, tjm, tjb);

        //W^k
        double Wk = Coulomb::Wk_abcd(K, w, b, v, m);

        //contribution to sum
        auto cont = e_denom * dv_ma * dv_ab * sixj * Wk * phase;
        A9 += cont;
      }
    }
  }
  A9 = -1 * (1 / (sqrt(2 * K + 1))) * pow(-1, K) * A9;
  return A9;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Tenth pair of diagrams
double A10(const int K, const int k1, const int k2,
           std::vector<Coulomb::meTable<double>> ME_tables,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
           const DiracSpinor &v, double omega) {
  //initialise A10
  double A10 = 0.0;

  //get tables
  Coulomb::meTable<double> met_1 = ME_tables.at(0);
  Coulomb::meTable<double> met_2 = ME_tables.at(1);

  //sum over core states, then excited states
  for (auto im = 0ul; im < excited.size(); ++im) {
    //get states
    const auto &m = excited[im];
    for (auto in = 0ul; in < excited.size(); ++in) {
      const auto &n = excited[in];
      for (auto ia = 0ul; ia < core.size(); ++ia) {
        // std::cout<<__LINE__<<" "<<ia<<"\n";
        //get states
        const auto &a = core[ia];
        //compute energy denominator
        double e_denom = 1 / (a.en() - m.en()) * 1 / (a.en() - n.en() + omega);
        //matrix elements (Can include RPA)

        const auto dv_ma = met_2.getv(m, a);
        const auto dv_nm = met_1.getv(n, m);

        //phase
        int kw = w.kappa();
        int ka = a.kappa();
        auto tjw = Angular::twoj_k(kw);
        auto tja = Angular::twoj_k(ka);
        auto phase = Angular::neg1pow_2(tjw + tja);
        //6j symbol
        int km = m.kappa();
        auto tjm = Angular::twoj_k(km);
        int kn = n.kappa();
        auto tjn = Angular::twoj_k(kn);
        double sixj = Angular::sixj_2(2 * K, 2 * k1, 2 * k2, tjm, tja, tjn);

        //W^k
        double Wk = Coulomb::Wk_abcd(K, w, a, v, n);

        //contribution to sum
        auto cont = e_denom * dv_ma * dv_nm * sixj * Wk * phase;
        A10 += cont;
      }
    }
  }
  A10 = (1 / (sqrt(2 * K + 1))) * pow(-1, k1 + k2) * A10;
  return A10;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// eleventh pair of diagrams
double A11(const int K, const int k1, const int k2,
           std::vector<Coulomb::meTable<double>> ME_tables,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
           const DiracSpinor &v, double omega) {
  //initialise A11
  double A11 = 0.0;

  //get tables
  Coulomb::meTable<double> met_1 = ME_tables.at(0);
  Coulomb::meTable<double> met_2 = ME_tables.at(1);

  //sum over core states, then excited states
  for (auto im = 0ul; im < excited.size(); ++im) {
    //get states
    const auto &m = excited[im];
    for (auto ib = 0ul; ib < core.size(); ++ib) {
      const auto &b = core[ib];
      for (auto ia = 0ul; ia < core.size(); ++ia) {
        // std::cout<<__LINE__<<" "<<ia<<"\n";
        //get states
        const auto &a = core[ia];
        //compute energy denominator
        double e_denom =
            1 / (a.en() - m.en() - omega) * 1 / (b.en() - m.en() - omega);
        //matrix elements (Can include RPA)

        const auto dv_ab = met_2.getv(a, b);
        const auto dv_bm = met_1.getv(b, m);

        //phase
        int kw = w.kappa();
        int km = m.kappa();
        auto tjw = Angular::twoj_k(kw);
        auto tjm = Angular::twoj_k(km);
        auto phase = Angular::neg1pow_2(tjw + tjm);

        //6j symbol
        int ka = a.kappa();
        auto tja = Angular::twoj_k(ka);
        int kb = b.kappa();
        auto tjb = Angular::twoj_k(kb);
        double sixj = Angular::sixj_2(2 * K, 2 * k1, 2 * k2, tjb, tja, tjm);

        //W^k
        double Wk = Coulomb::Wk_abcd(K, w, m, v, a);

        //contribution to sum
        auto cont = e_denom * dv_ab * dv_bm * sixj * Wk * phase;
        A11 += cont;
      }
    }
  }
  A11 = -1 * (1 / (sqrt(2 * K + 1))) * pow(-1, K) * A11;
  return A11;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// twelth pair of diagrams
double A12(const int K, const int k1, const int k2,
           std::vector<Coulomb::meTable<double>> ME_tables,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited, const DiracSpinor &w,
           const DiracSpinor &v, double omega) {
  //initialise A6
  double A12 = 0.0;

  //get tables
  Coulomb::meTable<double> met_1 = ME_tables.at(0);
  Coulomb::meTable<double> met_2 = ME_tables.at(1);

  //sum over core states, then excited states
  for (auto im = 0ul; im < excited.size(); ++im) {
    //get states
    const auto &m = excited[im];
    for (auto in = 0ul; in < excited.size(); ++in) {
      const auto &n = excited[in];
      for (auto ia = 0ul; ia < core.size(); ++ia) {
        // std::cout<<__LINE__<<" "<<ia<<"\n";
        //get states
        const auto &a = core[ia];
        //compute energy denominator
        double e_denom =
            1.0 / (a.en() - m.en() - omega) * 1.0 / (a.en() - n.en() - omega);
        //matrix elements (Can include RPA)

        const auto dv_nm = met_2.getv(n, m);
        const auto dv_an = met_1.getv(a, n);

        //phase
        int kw = w.kappa();
        int km = m.kappa();
        auto tjw = Angular::twoj_k(kw);
        auto tjm = Angular::twoj_k(km);
        auto phase = Angular::neg1pow_2(tjw + tjm);

        //6j symbol
        int ka = a.kappa();
        auto tja = Angular::twoj_k(ka);
        int kn = n.kappa();
        auto tjn = Angular::twoj_k(kn);
        double sixj = Angular::sixj_2(2 * K, 2 * k1, 2 * k2, tjn, tjm, tja);

        //W^k
        double Wk = Coulomb::Wk_abcd(K, w, m, v, a);

        //contribution to sum
        auto cont = e_denom * dv_nm * dv_an * sixj * Wk * phase;
        A12 += cont;
      }
    }
  }
  A12 = (1 / (sqrt(2 * K + 1))) * pow(-1, k1 + k2) * A12;
  return A12;
}

// No DCP for polarisability only
// <w||h1||n><n||h2||v>/dE + swap
double AK_pol(const int K, const int k1, const int k2,
              std::vector<Coulomb::meTable<double>> ME_tables,
              const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
              const DiracSpinor &v, double omega) {

  //get tables
  Coulomb::meTable<double> met_1 = ME_tables.at(0);
  Coulomb::meTable<double> met_2 = ME_tables.at(1);

  //sum excited states
  //initialise AK excited
  double AK_e = 0.0;
  for (auto im = 0ul; im < spectrum.size(); ++im) {
    //get states
    const auto &m = spectrum[im];

    // std::cout<<__LINE__<<" "<<ia<<"\n";
    //get states
    //compute energy denominator term
    double e_denom2 = 1 / (w.en() - m.en());
    double e_denom1 = pow(-1, K) / (v.en() - m.en());

    //matrix elements (Can include RPA)

    //Matrix elements
    const auto d1_wm = met_1.getv(w, m);
    const auto d1_mv = met_2.getv(m, v);

    const auto d2_wm = met_2.getv(w, m);
    const auto d2_mv = met_1.getv(m, v);

    //phase
    int kw = w.kappa();
    auto tjw = Angular::twoj_k(kw);
    auto phase = Angular::neg1pow_2(tjw);

    //6j symbol
    int kv = v.kappa();
    auto tjv = Angular::twoj_k(kv);
    int km = m.kappa();
    auto tjm = Angular::twoj_k(km);
    double sixj = Angular::sixj_2(tjm, tjw, 2 * k2, 2 * K, 2 * k1, tjv);

    //contribution to sum
    auto cont = e_denom1 * d1_wm * d1_mv * sixj * phase +
                e_denom2 * d2_wm * d2_mv * sixj * phase;
    AK_e += cont;
  }

  double AK = sqrt(2 * K + 1) * (AK_e);
  return AK;
}

//----------------------------------------------------------------------------------------------------------------------------------------------
// Corrections to relevant physical quantities depend on some additional scaling factor
void Polarisability(double diagrams, std::string static_operator, int K,
                    double omega) {
  //for pnc E1 amplitude
  if (static_operator == "pnc") {
    double threej1 = Angular::threej_2(0.0, 2.0, 2.0, 0.0, 0.0, 0.0);
    double threej2 = Angular::threej_2(1, 2, 1, 1, 0, -1);

    double PNC_amplitude_correction = -1.0 * threej1 * threej2 * diagrams;
    //double PNC_amplitude_correction = diagrams;
    std::cout << "The DCP correction to the PNC amplitude is "
              << PNC_amplitude_correction << " iea_B(-Q_W/N) 10^{-11}.\n";
  }
  //For Stark induced amplitude: transition so omega non zero
  if (static_operator == "E1") {
    //for scalaer case
    if (K == 0) {
      if (omega == 0.0) {
        double alpha_correction = (1.0 / sqrt(6)) * diagrams;
        std::cout << "The scalar static dipole "
                     "polarisability is "
                  << alpha_correction << " a_B^3.\n";

      } else {
        double alpha_correction = -(1.0 / sqrt(6.0)) * diagrams;
        std::cout << "The  scalar transition polarisability is "
                  << alpha_correction << " a_B^3.\n";
      }
    } else if (K == 1) {
      double beta_correction = (1.0 / 2.0) * diagrams;
      std::cout << "The vector transition polarisability is " << beta_correction
                << " a_B^3.\n";
    }
  }
};

//===========================================================================================================================================================

//Calculates contirbution from all 24 DCP diagrams at third order in the residual Coulomb interaction and sums to obtain
// total amplitude

double delta_A_KQ(const int K, const int k1, const int k2,
                  std::vector<Coulomb::meTable<double>> Table,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const DiracSpinor &Fw, const DiracSpinor &Fv, double omega) {

  double A1_eval = A1(K, k1, k2, Table, core, excited, Fw, Fv, omega);

  double A2_eval = A2(K, k1, k2, Table, core, excited, Fw, Fv, omega);

  double A3_eval = A3(K, k1, k2, Table, core, excited, Fw, Fv, omega);

  double A4_eval = A4(K, k1, k2, Table, core, excited, Fw, Fv, omega);

  double A5_eval = A5(K, k1, k2, Table, core, excited, Fw, Fv, omega);

  double A6_eval = A6(K, k1, k2, Table, core, excited, Fw, Fv, omega);

  double A7_eval = A7(K, k1, k2, Table, core, excited, Fw, Fv, omega);

  double A8_eval = A8(K, k1, k2, Table, core, excited, Fw, Fv, omega);

  double A9_eval = A9(K, k1, k2, Table, core, excited, Fw, Fv, omega);

  double A10_eval = A10(K, k1, k2, Table, core, excited, Fw, Fv, omega);

  double A11_eval = A11(K, k1, k2, Table, core, excited, Fw, Fv, omega);

  double A12_eval = A12(K, k1, k2, Table, core, excited, Fw, Fv, omega);

  double sum = A1_eval + A2_eval + A3_eval + A4_eval + A5_eval + A6_eval +
               A7_eval + A8_eval + A9_eval + A10_eval + A11_eval + A12_eval;
  return sum;
}

//----------------------------------------------------------------------------------------------------------------------------------------------
// Corrections to relevant physical quantities depend on some additional scaling factor
void DCP_correction(double diagrams, std::string static_operator, int K,
                    double omega) {
  //for pnc E1 amplitude
  if (static_operator == "pnc") {
    double threej1 = Angular::threej_2(0.0, 2.0, 2.0, 0.0, 0.0, 0.0);
    double threej2 = Angular::threej_2(1, 2, 1, 1, 0, -1);

    double PNC_amplitude_correction = -1 * threej1 * threej2 * diagrams;
    //double PNC_amplitude_correction = diagrams;
    std::cout << "The DCP correction to the PNC amplitude is "
              << PNC_amplitude_correction << " iea_B(-Q_W/N) 10^{-11}.\n";
  }
  //For Stark induced amplitude
  if (static_operator == "E1") {
    if (K == 0) {
      if (omega == 0.0) {
        double threej3 = Angular::threej_2(2.0, 2.0, 0.0, 0.0, 0.0, 0.0);
        double threej4 = Angular::threej_2(1, 0, 1, 1, 0, -1);
        double alpha_correction = (1 / sqrt(6)) * diagrams;
        std::cout << "The DCP correction to the scalar static dipole "
                     "polarisability is "
                  << alpha_correction << " a_B^3.\n"
                  << threej3 * threej4;

      } else {
        double alpha_correction = -(1.0 / sqrt(6.0)) * diagrams;
        std::cout
            << "The DCP correction to the scalar transition polarisability is "
            << alpha_correction << " a_B^3.\n";
      }
    } else if (K == 1) {
      double beta_correction = (1.0 / 2.0) * diagrams;
      std::cout
          << "The DCP correction to the vector transition polarisability is "
          << beta_correction << " a_B^3.\n";
    }
  }
};

void DCPdiagrams(const IO::InputBlock &input, const Wavefunction &wf) {
  // This is a module to compute the diagrams contributing to the reduced amplitude due to the double core polarisation
  // in the presence of two external fields.

  //Takes as input:
  // Input block and wavefunctions

  //Computes DCP correction to alpha and beta using Diagram method

  input.check(
      {{"transition", "List. states (e.g., 6s,6s) []"},
       {"operator_1", "e.g., E1, hfs (see ampsci -o for available operators)"},
       {"operator_2", "e.g., E1, hfs (see ampsci -o for available operators)"},
       {"K", "Irreducible componenet K=0,1,2"},
       {"rpa", "Include RPA? [true]"},
       {"omega", "frequency (for single w) [default: transition freq.]"},
       {"Qk_file",
        "true/false/filename - SR: filename for QkTable file. If blank will "
        "not use QkTable; if exists, will read it in; if doesn't exist, will "
        "create it and write to disk. If 'true' will use default filename. "
        "Save time (10x) at cost of memory. Note: Using QkTable "
        "implies splines used for diagram legs"}});
  // If we are just requesting 'help', don't run module:
  double test = Angular::threej_2(0.0, 2.0, 2.0, 0.0, 0.0, 0.0);
  std::cout << test << "\n\n\n\n";

  double test2 = Angular::threej_2(1.0, 2.0, 1.0, 1.0, 0.0, -1.0);
  std::cout << test2 << "\n\n\n\n";

  if (input.has_option("help")) {
    return;
  }

  const auto K = input.get("K", 0);

  //get states
  const auto states = input.get("transition", std::vector<std::string>{});
  if (states.size() != 2) {
    std::cout << "Error 491 in transitionPolarisability(): transition option "
                 "must have exactly two states comma-separated\n";
    return;
  }
  const auto pv = wf.getState(states.at(0));
  const auto pw = wf.getState(states.at(1));
  if (!pv) {
    std::cout << "Error: Couldn't find state: " << states.at(0) << "?\n";
    return;
  }
  if (!pw) {
    std::cout << "Error: Couldn't find state: " << states.at(1) << "?\n";
    return;
  }
  const auto &Fv = *pv;
  const auto &Fw = *pw;

  const auto omega_default = std::abs(Fv.en() - Fw.en());
  std::cout << omega_default << "\n\n\n";

  // To compute static dipole corrections we have omega=0;
  const auto omega = input.get("omega", omega_default);

  //separate core and excited orbitals using split by energy (fermi)
  auto basis = wf.basis();
  auto spectrum = wf.spectrum();

  double fermi = wf.FermiLevel();

  std::pair<std::vector<DiracSpinor>, std::vector<DiracSpinor>> split =
      DiracSpinor::split_by_energy(basis, fermi, 1);

  std::vector<DiracSpinor> core = split.first;
  std::vector<DiracSpinor> excited = split.second;

  //get ranks of operators
  const auto oper1 = input.get<std::string>("operator_1", "");

  // Get optional 'options' for operator
  auto h1_options = IO::InputBlock(oper1, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h1_options = *tmp_opt;
  }
  const auto h1 = DiracOperator::generate(oper1, h1_options, wf);

  const auto oper2 = input.get<std::string>("operator_2", "");
  // Get optional 'options' for operator
  auto h2_options = IO::InputBlock(oper2, {});
  const auto tmp_opt2 = input.getBlock("options");
  if (tmp_opt2) {
    h2_options = *tmp_opt2;
  }
  const auto h2 = DiracOperator::generate(oper2, h2_options, wf);

  const auto k1 = h1->rank();
  const auto k2 = h2->rank();

  int kv = Fv.kappa();
  //auto jv = 0.5 * Angular::twoj_k(kv);

  // If including RPA, solve TDHF equations for polarisation of core electrons
  const auto rpaQ = input.get("rpa", true);
  // Construct TDHF operator
  std::cout << "test.... " << pow(-1, 1) << "\n";

  if (rpaQ) {

    std::cout
        << "Including RPA/core polarisation correction to matrix elements: ";

    auto dV1 = ExternalField::TDHF(h1.get(), wf.vHF());
    auto dV2 = ExternalField::TDHF(h2.get(), wf.vHF());
    std::cout << "Solving TDHF equations of core using operator 1.... ";
    dV1.solve_core(omega, 99);
    std::cout << "Complete. \n";
    std::cout << "Solving TDHF equations of core using operator 2.... ";
    dV2.solve_core(0.0, 99);
    std::cout << "Complete. \n";

    // Compute and store matrix elements with RPA
    std::cout << "Computing and storing all matrix elements.... ";
    std::vector<Coulomb::meTable<double>> Table =
        compute_me(h1.get(), h2.get(), core, excited, omega, &dV1, &dV2);
    std::cout << "Complete. \n";

    // Compute Diagram DCP corrections

    double DCP = delta_A_KQ(K, k1, k2, Table, core, excited, Fw, Fv, omega);
    DCP_correction(DCP, oper2, K, omega);

    //Compute quantity without dcp correction to check
    std::cout << "\n\n\n Calculating without DCP correction... \n\n";
    double Ak_pol = AK_pol(K, k1, k2, Table, spectrum, Fw, Fv, omega);
    Polarisability(Ak_pol, oper2, K, omega);
  } else {
    // Compute and store matrix elements without RPA
    std::cout << "Computing and storing all matrix elements without RPA.... ";
    std::vector<Coulomb::meTable<double>> Table =
        compute_me(h1.get(), h2.get(), spectrum, wf.valence(), omega);
    std::cout << "Complete. \n";

    // Compute Diagram DCP corrections

    double DCP = delta_A_KQ(K, k1, k2, Table, core, excited, Fw, Fv, omega);
    DCP_correction(DCP, oper2, K, omega);

    //Compute quantity without dcp correction to check
    std::cout << "\n\n\n Calculating without DCP correction... \n\n";
    double Ak_pol = AK_pol(K, k1, k2, Table, spectrum, Fw, Fv, omega);
    Polarisability(Ak_pol, oper2, K, omega);
  }
}

// sum contribution from all diagrams

// Compute contribution to scalar/vector transition polarisability through appropriate scaling factor
} // namespace Module
// namespace Module
// namespace Module
