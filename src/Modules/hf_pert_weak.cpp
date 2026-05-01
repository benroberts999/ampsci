//Module for computing the third order hyperfine-weak induced electric dipole tranisition matrix element required for extraction
// of the nuclear anapole moment from atomic parity violation experiments
#include "Modules/hf_pert_weak.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/meTable.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "DiracOperator/include.hpp" //For operators
#include "DiracOperator/include.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/matrixElements.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/NuclearData.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Physics/PhysConst_constants.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <optional>
#include <string>
#include <vector>

namespace Module {

// First we compute all of the required reduced matrix elements needed
// This includes matrix elements corresponding to the dipole, hyperfine in NSI weak interaction operators
// Matrix elements filling the tables can be computed with or wirhout core polarisation RPA corrections

// Brueckner correlations can be included to second/all orders with fitting throught the spectrum and valence states

// Inputs to compute matrix elements

//operators
// TDHF objects
//spectrum
//valence transition states

std::vector<Coulomb::meTable<double>>
compute_me_3f(const DiracOperator::TensorOperator *hpnc,
              const DiracOperator::TensorOperator *he1,
              const DiracOperator::TensorOperator *const hfs,
              const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
              const DiracSpinor &v,
              const ExternalField::CorePolarisation *dV_pnc,
              const ExternalField::CorePolarisation *dV_e1,
              const ExternalField::DiagramRPA *dV_hf) {

  //Need matrix elements between valence states in transition and intermediate states from the spectrum

  //create tables for matrix elements
  Coulomb::meTable<double> pnc_me; //PNC
  Coulomb::meTable<double> e1_me;  //E1
  Coulomb::meTable<double> hf_me;  // Magnetic dipole hyperfine

  // initialise vector of tables for output
  std::vector<Coulomb::meTable<double>> ME_tables;

  //fill tables: need matrix elements

  //summing over all states in the spectrum

  //compute all matrix elements -> including RPA if dV object is passed based on input options
  for (const auto &i : spectrum) {

    //first we consider matrix elements involving the initial or final states
    //compute me
    //pnc matrix elements--------------------------------------------------------------------
    const auto pnc_wi =
        hpnc->reducedME(w, i) + (dV_pnc ? dV_pnc->dV(w, i) : 0.0);
    const auto pnc_iv =
        hpnc->reducedME(i, v) + (dV_pnc ? dV_pnc->dV(i, v) : 0.0);
    const auto w_wf = w.shortSymbol();
    const auto i_wf = i.shortSymbol();
    const auto v_wf = v.shortSymbol();

    //output for test
    /*
    if (i.n() < 16) {
      std::cout << "Reduced NSI PNC matrix element  " << " <<" << w_wf
                << "||h_pnc||" << i_wf << ">> = " << pnc_wi
                << "(-Q_w/N)x10^-11 \n.";
    }
    const auto pnc_iv =
        hpnc->reducedME(i, v) + (dV_pnc ? dV_pnc->dV(i, v) : 0.0);
    if (i.n() < 16) {
      std::cout << "Reduced NSI PNC matrix element  " << " <<" << i_wf
                << "||h_pnc||" << v_wf << ">> = " << pnc_iv
                << "(-Q_w/N)x10^-11\n.";
    }
                */
    //add to table
    pnc_me.add(w, i, pnc_wi);
    pnc_me.add(i, v, pnc_iv);

    //e1 matrix elements--------------------------------------------------------------------
    const auto e1_wi = he1->reducedME(w, i) + (dV_e1 ? dV_e1->dV(w, i) : 0.0);
    const auto e1_iv = he1->reducedME(i, v) + (dV_e1 ? dV_e1->dV(i, v) : 0.0);
    //add to table
    e1_me.add(w, i, e1_wi);
    e1_me.add(i, v, e1_iv);

    //hfs matrix elements--------------------------------------------------------------------
    const auto hfs_wi = hfs->reducedME(w, i) + (dV_hf ? dV_hf->dV(w, i) : 0.0);
    const auto hfs_iv = hfs->reducedME(i, v) + (dV_hf ? dV_hf->dV(i, v) : 0.0);
    //add to table
    hf_me.add(w, i, hfs_wi);
    hf_me.add(i, v, hfs_iv);

    //output for test
    /*
    if (i.n() < 16) {
      std::cout << "Reduced hyperfine matrix element  " << " <<" << w_wf
                << "||h_hfs||" << i_wf << ">> = " << hfs_wi << "MHz \n.";
    }
    if (i.n() < 16) {
      std::cout << "Reduced hyperfine matrix element " << " <<" << i_wf
                << "||h_hfs||" << v_wf << ">> = " << hfs_iv << "MHz \n.";
    }
*/
    //now include additional loop over spectrum
    for (const auto &j : spectrum) {
      const auto pnc_ij =
          hpnc->reducedME(i, j) + (dV_pnc ? dV_pnc->dV(i, j) : 0.0);
      const auto e1_ij = he1->reducedME(i, j) + (dV_e1 ? dV_e1->dV(i, j) : 0.0);
      const auto hfs_ij =
          hfs->reducedME(i, j) + (dV_hf ? dV_hf->dV(i, j) : 0.0);

      //add to table
      //hpnc->rme3js(i.twoj(), i.twoj();
      pnc_me.add(i, j, pnc_ij);
      e1_me.add(i, j, e1_ij);
      hf_me.add(i, j, hfs_ij);
    }
  }
  //finally create vector of matrix element tables for single element output
  ME_tables.push_back(pnc_me);
  ME_tables.push_back(e1_me);
  ME_tables.push_back(hf_me);
  return ME_tables;
}

// To compute the reduced matrix element we split the total quantity into two summations

// <w||h_[w+hf]|v>=h_1+h_2
//Inputs
// All matrix elements for each operator
// initial and final wavefunction
//spectrum
//nuclear spin
//Total initial and final angular momenta (F=I+J)
double h1(std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
          const DiracSpinor &v, int I2, int Fv2, int Fw2, int two_k) {
  double h1 = 0.0;
  //get angular momenta of valence electronic states from wavefunctions (we actually extract 2j) from kappa
  int kw = w.kappa();
  int kv = v.kappa();
  auto tjw = Angular::twoj_k(kw);
  auto tjv = Angular::twoj_k(kv);

  //get energies for valence states
  double e_v = v.en();
  double e_w = w.en();

  //get tables
  Coulomb::meTable<double> pnc_me = ME_tables.at(0);
  Coulomb::meTable<double> e1_me = ME_tables.at(1);
  Coulomb::meTable<double> hf_me = ME_tables.at(2);

  for (const auto &j : spectrum) {
    // make sure j is not equal to v
    if (j != v) {

      //first compute angular factor
      auto tjj = Angular::twoj_k(j.kappa());

      auto phase = Angular::neg1pow_2(tjv - tjw + 2);
      double sixj1 = Angular::sixj_2(Fw2, Fv2, 2, tjj, tjw, I2);
      double sixj2 = Angular::sixj_2(I2, I2, two_k, tjj, tjv, Fv2);
      double angular = phase * sixj1 * sixj2;

      //initialise sums over i for given j
      double sum_1 = 0.0;
      double sum_2 = 0.0;
      double sum_3 = 0.0;

      //get energy for state j
      double e_j = j.en();

      //second loop
      for (const auto &i : spectrum) {
        if (i != w) {
          //get energy for state i
          double e_i = i.en();

          //first term
          double e_denom1 = (e_j - e_v) * (e_i - e_w);
          double t1 = pnc_me.getv(w, i) * e1_me.getv(i, j) * hf_me.getv(j, v) /
                      e_denom1;

          //second term
          double e_denom2 = (e_j - e_v) * (e_i - e_v);
          double t2 = pnc_me.getv(i, j) * e1_me.getv(w, i) * hf_me.getv(j, v) /
                      e_denom2;

          //third term
          double e_denom3 = (e_j - e_v) * (e_i - e_w);
          double t3 = pnc_me.getv(i, v) * e1_me.getv(w, j) * hf_me.getv(j, i) /
                      e_denom3;

          sum_1 = sum_1 + t1;
          sum_2 = sum_2 + t2;
          sum_3 = sum_3 + t3;
        }
      }
      //fourth term (no summation over i)
      double e_denom4 = (e_j - e_v) * (e_j - e_v);
      double t4 =
          pnc_me.getv(j, v) * e1_me.getv(w, j) * hf_me.getv(v, v) / e_denom4;
      h1 = h1 + angular * (sum_1 + sum_2 + sum_3 - t4);
    }
  }
  return h1;
}

double h2(std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
          const DiracSpinor &v, int I2, int Fv2, int Fw2, int two_k) {
  double h2 = 0.0;
  //get angular momenta of valence electronic states from wavefunctions (we actually extract 2j) from kappa
  int kw = w.kappa();
  int kv = v.kappa();
  auto tjw = Angular::twoj_k(kw);
  auto tjv = Angular::twoj_k(kv);

  //get energies for valence states
  double e_v = v.en();
  double e_w = w.en();

  //get tables
  Coulomb::meTable<double> pnc_me = ME_tables.at(0);
  Coulomb::meTable<double> e1_me = ME_tables.at(1);
  Coulomb::meTable<double> hf_me = ME_tables.at(2);

  for (const auto &j : spectrum) {
    // make sure j is not equal to v
    if (j != w) {

      //first compute angular factor
      auto tjj = Angular::twoj_k(j.kappa());

      auto phase = Angular::neg1pow_2(Fv2 - Fw2 + 2);
      double sixj1 = Angular::sixj_2(Fw2, Fv2, 2, tjv, tjj, I2);
      double sixj2 = Angular::sixj_2(I2, I2, two_k, tjj, tjw, Fw2);
      double angular = phase * sixj1 * sixj2;

      //initialise sums over i for given j
      double sum_1 = 0.0;
      double sum_2 = 0.0;
      double sum_3 = 0.0;

      //get energy for state j
      double e_j = j.en();

      //second loop
      for (const auto &i : spectrum) {
        if (i != v) {
          //get energy for state i
          double e_i = i.en();

          //first term
          double e_denom1 = (e_i - e_v) * (e_j - e_w);
          double t1 = pnc_me.getv(i, v) * e1_me.getv(j, i) * hf_me.getv(w, j) /
                      e_denom1;

          //second term
          double e_denom2 = (e_j - e_w) * (e_i - e_w);
          double t2 = pnc_me.getv(w, i) * e1_me.getv(j, v) * hf_me.getv(i, j) /
                      e_denom2;

          //third term
          double e_denom3 = (e_j - e_w) * (e_i - e_w);
          double t3 = pnc_me.getv(j, i) * e1_me.getv(i, v) * hf_me.getv(w, j) /
                      e_denom3;

          sum_1 = sum_1 + t1;
          sum_2 = sum_2 + t2;
          sum_3 = sum_3 + t3;
        }
      }
      //fourth term (no summation over i)
      double e_denom4 = (e_j - e_w) * (e_j - e_w);
      double t4 =
          pnc_me.getv(w, j) * e1_me.getv(j, v) * hf_me.getv(w, w) / e_denom4;
      h2 = h2 + angular * (sum_1 + sum_2 + sum_3 - t4);
    }
  }
  return h2;
}

void hf_pert_weak(const IO::InputBlock &input, const Wavefunction &wf) {
  // This function takes input and computes the weak perturbed hyperfine transition amplitude

  //Takes as input:
  // Input block and wavefunctions

  //check the input options
  input.check({{"transition", "List. states (e.g., 6s,6s) []"},
               {"rpa", "Include RPA? [true]"},
               {"two_I", "two times the nuclear spin (integer)"},
               {"two_Fw", "two times total angular momentum of final state w"},
               {"two_Fv", "two times total angular momentum of final state v"},
               {"hfs_options{}", "Options for HFS operator (see -o hfs)"},
               {"two_k", "two times hyperfine multipolarity"}});
  // If we are just requesting 'help', don't run module:

  if (input.has_option("help")) {
    return;
  }

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

  //   const auto omega = Fw.en() - Fv.en();

  double tjw = Fw.twoj();
  std::cout << sqrt(tjw + 1) << "\n\n\n\n";
  //get spectrum
  // We should use _spectrum_ for the sos - but if it is empty, just use basis
  auto spectrum = wf.spectrum().empty() ? wf.basis() : wf.spectrum();

  //generate operators

  // input: nuc parameters for rho for pnc operator:
  const auto c_dflt =
      Nuclear::c_hdr_formula_rrms_t(Nuclear::find_rrms(wf.Znuc(), wf.Anuc()));
  const auto t = input.get("t", Nuclear::default_t);
  const auto c = input.get("c", c_dflt);

  //parameters for hyperfine

  // Generate operators:
  const auto N_nuc = wf.Anuc() - wf.Znuc();
  //NSI weak interaction
  DiracOperator::PNCnsi hpnc(c, t, wf.grid(), -N_nuc, "i(Qw/-N)*e-11");

  //E1
  DiracOperator::E1 he1(wf.grid());

  // Generate hyperfine operator with default options for now (magnetic dipole)
  const auto hfs_options = input.getBlock("hfs_options");

  // Generate hyperfine operator (including BW effect from input)
  const auto hfs = DiracOperator::generate(
      "hfs", hfs_options ? *hfs_options : IO::InputBlock{}, wf);

  // If including RPA, solve TDHF equations for polarisation of core electrons

  const auto rpaQ = input.get("rpa", true);
  if (rpaQ) {
    //for E1 and pnc we use TDHF method for rpa
    auto dVE1 = ExternalField::TDHF(&he1, wf.vHF());
    auto dVpnc = ExternalField::TDHF(&hpnc, wf.vHF());
    const auto omega_dflt = Fw.en() - Fv.en();
    const auto omega = input.get("omega", omega_dflt);
    auto E1_it = input.get("E1_rpa_it", 99);
    auto pnc_it = input.get("pnc_rpa_it", 99);
    dVE1.solve_core(omega, E1_it);
    dVpnc.solve_core(0.0, pnc_it);
    const auto two_k = input.get("two_k", 0);

    /*
    if(two_k==4){
      auto dVhf = ExternalField::TDHF(hfs.get(), wf.vHF());
      auto hfs_it = input.get("hfs_rpa_it", 99);
      dVpnc.solve_core(0.0, hfs_it);

    }
      */

      std::cout << "\nIncluding RPA (diagram method) - must have basis\n";
    std::unique_ptr<ExternalField::DiagramRPA> rpa_hf;
    rpa_hf = std::make_unique<ExternalField::DiagramRPA>(
        hfs.get(), wf.basis(), wf.vHF(), wf.identity());
    rpa_hf->solve_core(0.0, 100, true);

    // Compute and store matrix elements with RPA
    std::cout << "Computing and storing all matrix elements.... ";
    std::vector<Coulomb::meTable<double>> Table = compute_me_3f(
        &hpnc, &he1, hfs.get(), spectrum, Fw, Fv, &dVpnc, &dVE1, rpa_hf.get());
    std::cout << "Complete. \n";

    //use matrix elements to compute the amplitude
    int I2 = input.get("two_I", 0);
    const auto Fv2 = input.get("two_Fv", 0);
    const auto Fw2 = input.get("two_Fw", 0);
    const auto two_k = input.get("two_k", 0);

    double h_1 = h1(Table, spectrum, Fw, Fv, I2, Fv2, Fw2, two_k);

    double h_2 = h2(Table, spectrum, Fw, Fv, I2, Fv2, Fw2, two_k);

    double I = 0.5 * I2;

    double h_sum =
        sqrt(I * (I + 1) * (I2 + 1) * (Fv2 + 1) * (Fw2 + 1)) * (h_1 + h_2);

    //need additional factor to account for different definition of the reduced matrix element (phase+3j symbol)
    // double tjw = Fw.twoj();
    double h_final = h_sum;
    std::cout << "value is " << h_final << "  MHz .with RPA\n\n";
    std::cout << "In atomic unuts value is " << h_final / PhysConst::Hartree_MHz
              << "  au .";

  } else {
    // Compute and store matrix elements without RPA
    std::cout << "Computing and storing all matrix elements without RPA.... ";
    std::vector<Coulomb::meTable<double>> Table =
        compute_me_3f(&hpnc, &he1, hfs.get(), spectrum, Fw, Fv);
    std::cout << "Complete. \n";

    //use matrix elements to compute the amplitude
    int I2 = input.get("two_I", 0);
    const auto Fv2 = input.get("two_Fv", 0);
    const auto Fw2 = input.get("two_Fw", 0);
    const auto two_k = input.get("two_k", 0);

    double h_1 = h1(Table, spectrum, Fw, Fv, I2, Fv2, Fw2, two_k);

    double h_2 = h2(Table, spectrum, Fw, Fv, I2, Fv2, Fw2, two_k);

    double I = 0.5 * I2;

    double h_sum =
        sqrt(I * (I + 1) * (I2 + 1) * (Fv2 + 1) * (Fw2 + 1)) * (h_1 + h_2);

    //need additional factor to account for different definition of the reduced matrix element (phase+3j symbol)
    // double tjw = Fw.twoj();
    double h_final = h_sum;
    std::cout << "value is " << h_final << "  MHz .without RPA\n\n";
    std::cout << "In atomic unuts value is " << h_final / PhysConst::Hartree_MHz
              << "  au .";
  }

  std::cout << "\n\n\n test succesful!... ";
}

// namespace Module

// namespace Module
} // namespace Module