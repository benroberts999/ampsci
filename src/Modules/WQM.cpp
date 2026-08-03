//Module for computing pnc amplitude induced by WQM
#include "Modules/Modules.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/meTable.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "DiracOperator/include.hpp" //For operators
#include "DiracOperator/include.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
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
#include <iostream>
#include <algorithm>

namespace Module {

// First we compute all of the required reduced matrix elements needed
// This includes matrix elements corresponding to the dipole and sd weak interaction matrix elements
// Matrix elements filling the tables can be computed with or wirhout core polarisation RPA corrections

// Brueckner correlations can be included to second/all orders with fitting throught the spectrum and valence states

// Inputs to compute matrix elements

//operators
// TDHF objects
//spectrum
//valence transition states


std::vector<Coulomb::meTable<double>>
compute_me(const DiracOperator::TensorOperator *hpnc,
              const DiracOperator::TensorOperator *he1,
              const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
              const DiracSpinor &v,
              const ExternalField::CorePolarisation *dV_pnc=nullptr,
              const ExternalField::CorePolarisation *dV_e1=nullptr)
 {

  //Need matrix elements between valence states and intermediate states

  //create tables for matrix elements
  Coulomb::meTable<double> pnc_me; //PNC
  Coulomb::meTable<double> e1_me;  //E1

  // initialise vector of tables for output
  std::vector<Coulomb::meTable<double>> ME_tables;

  //fill tables: need matrix elements

  //summing over all states in the spectrum

  



  //compute all matrix elements -> including RPA if dV object is passed based on input options
  for (const auto &n : spectrum) {
  


  //specify pnc matrix element definition -- 
    //pnc matrix elements--------------------------------------------------------------------
    const auto pnc_wn =
        hpnc-> reducedME(w, n) + (dV_pnc ? dV_pnc->dV(w, n) : 0.0);
    const auto pnc_nv =
        hpnc-> reducedME(n, v) + (dV_pnc ? dV_pnc->dV(n, v) : 0.0);
    const auto w_wf = w.shortSymbol();
    const auto n_wf = n.shortSymbol();
    const auto v_wf = v.shortSymbol();

    /*
    if (i.n() < 0) {
      std::cout << "Reduced NSI PNC matrix element  " << " <<" << w_wf
                << "||h_pnc||" << i_wf << ">> = " << pnc_wi
                << "(-Q_w/N)x10^-11 \n.";
    }

    if (i.n() < 0) {
      std::cout << "Reduced NSI PNC matrix element  " << " <<" << i_wf
                << "||h_pnc||" << v_wf << ">> = " << pnc_iv
                << "(-Q_w/N)x10^-11\n.";
    }
  */
    pnc_me.add(w, n, pnc_wn);
    pnc_me.add(n, v, pnc_nv);

    //e1 matrix elements--------------------------------------------------------------------
    const auto e1_wn = he1->reducedME(w, n) + (dV_e1 ? dV_e1->dV(w, n) : 0.0);
    const auto e1_nv = he1->reducedME(n, v) + (dV_e1 ? dV_e1->dV(n, v) : 0.0);
    //add to table
    e1_me.add(w, n, e1_wn);
    e1_me.add(n, v, e1_nv);

    //output for test
    /*
    if (i.n() < 0) {
      std::cout << "Reduced hyperfine matrix element  " << " <<" << w_wf
                << "||h_hfs||" << i_wf << ">> = " << hfs_wi << "MHz \n.";
    }
    if (i.n() < 0) {
      std::cout << "Reduced hyperfine matrix element " << " <<" << i_wf
                << "||h_hfs||" << v_wf << ">> = " << hfs_iv << "MHz \n.";
    }
    */            
  }
  //finally create vector of matrix element tables for single element output
  ME_tables.push_back(pnc_me);
  ME_tables.push_back(e1_me);
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
double pnc_WQM(std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
          const DiracSpinor &v, int I2, int Fv2, int Fw2) {

 

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

  //initialise sums
            
  double sum1=0.0;
  double sum2=0.0;



  for (const auto &n : spectrum) {
    // make sure n is not equal to v
    //get properties for state n
    double e_n = n.en();
    int kn = n.kappa();
    auto tjn = Angular::twoj_k(kn);


    //first sum
    if (n == v) {
      //do nothing
    }else{
      //sum 1
      auto phase1 = Angular::neg1pow_2(tjw - tjv);
      double sixj1 = Angular::sixj_2(tjn, tjv, 4, I2, I2, Fv2);
      double sixj2 = Angular::sixj_2(tjn, tjw, 2, Fw2, Fv2, I2);
      double angular1 = phase1 * sixj1 * sixj2;
      double e_denom = 1.0/(e_n-e_v);
      double term_n = angular1*e1_me.getv(w,n)*pnc_me.getv(n, v)*e_denom;
      sum1=sum1+term_n;
    }

    //second sum
    if(n==w){
      //do nothing
    }else{
      auto phase2 = Angular::neg1pow_2(Fw2 - Fv2);
      double sixj3 = Angular::sixj_2(tjn, tjw, 4, I2, I2, Fw2);
      double sixj4 = Angular::sixj_2(tjn, tjv, 2, Fv2, Fw2, I2);
      double angular2 = phase2 * sixj3 * sixj4;
      double e_denom = 1.0/(e_n-e_w);
      double term_n = angular2*e1_me.getv(n,v)*pnc_me.getv(w, n)*e_denom;
      sum2=sum2+term_n;
    }
    }

  double amplitude = sum1+sum2;

  return amplitude;
}


void WQM(const IO::InputBlock &input, const Wavefunction &wf) {
  // This function takes input and computes the weak perturbed hyperfine transition amplitude

  //Takes as input:
  // Input block and wavefunctions

  //check the input options
  input.check({{"transition", "List. states (e.g., 6s,6s) []"},
               {"rpa", "Include RPA? [true]"},
               {"two_I", "two times the nuclear spin (integer)"},
               {"two_Fw", "two times total angular momentum of final state w"},
               {"two_Fv", "two times total angular momentum of final state v"}});
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
  const auto r_rms =wf.get_rrms();


  // Generate operators:
  //const auto N_nuc = wf.Anuc() - wf.Znuc();
  //WQM interaction
  DiracOperator::PNC_wqm hpnc(c, t,r_rms, wf.grid());

  //E1
  DiracOperator::E1 he1(wf.grid());


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

    /*
    if(two_k==4){
      auto dVhf = ExternalField::TDHF(hfs.get(), wf.vHF());
      auto hfs_it = input.get("hfs_rpa_it", 99);
      dVpnc.solve_core(0.0, hfs_it);

    }
      */

    // Compute and store matrix elements with RPA
    std::cout << "Computing and storing all matrix elements.... ";
    std::vector<Coulomb::meTable<double>> Table = compute_me(
        &hpnc, &he1, spectrum, Fw, Fv,&dVpnc, &dVE1);
    std::cout << "Complete. \n";

    //use matrix elements to compute the amplitude
    int I2 = input.get("two_I", 0);
    const auto Fv2 = input.get("two_Fv", 0);
    const auto Fw2 = input.get("two_Fw", 0);
    //const auto two_k = input.get("two_k", 0);

    double  pnc = pnc_WQM(Table, spectrum, Fw, Fv, I2, Fv2, Fw2);

    double I = 0.5 * I2;

    //return final values 
    //first negative and positive contributiions separately

    //pre factors depend on respective operators 

    double pre_factor=sqrt((I2+3)*(I2+1)*(I+1)/(I*(I2-1)))*sqrt((Fw2+1)*(Fv2+1));

    int two_M = std::min(Fw2, Fv2); 


    double z_me_conv =Angular::neg1pow_2(Fw2 - two_M)* Angular::threej_2(Fw2, 2, Fv2, -two_M, 0, two_M);


    double pnc_final = pre_factor * pnc;
    double pnc_fullz = z_me_conv*pnc_final;

   
    std::cout << "\n The reduced matrix element is " << pnc_final << "iea_0 Q^TW e-11\n\n";
    
    std::cout << "\n The z-component of the dipole matrix element is " << pnc_fullz << "iea_0 Q^TW e-11\n\n";

  } else {
    // Compute and store matrix elements without RPA
    std::cout << "Computing and storing all matrix elements without RPA.... ";
    std::vector<Coulomb::meTable<double>> Table =
        compute_me(&hpnc, &he1, spectrum, Fw, Fv);
    std::cout << "Complete. \n";

    //use matrix elements to compute the amplitude
    int I2 = input.get("two_I", 0);
    const auto Fv2 = input.get("two_Fv", 0);
    const auto Fw2 = input.get("two_Fw", 0);

    double  pnc = pnc_WQM(Table, spectrum, Fw, Fv, I2, Fv2, Fw2);

    double I = 0.5 * I2;

    double pre_factor=sqrt((I2+3)*(I2+1)*(I+1)/(I*(I2-1)))*sqrt((Fw2+1)*(Fv2+1));

    int two_M = std::min(Fw2, Fv2); 

    double z_me_conv =Angular::neg1pow_2(Fw2 - two_M)* Angular::threej_2(Fw2, 2, Fv2, -two_M, 0, two_M);

    double pnc_final = pre_factor * pnc;
    double pnc_fullz = z_me_conv*pnc_final;

    //need additional factor to account for different definition of the reduced matrix element (phase+3j symbol)
    // double tjw = Fw.twoj();
    std::cout << "\n The reduced matrix element is " << pnc_final << "iea_0 Q^TW e-11\n\n";
    
    std::cout << "\n The z-component of the dipole matrix element is " << pnc_fullz << "iea_0 Q^TW e-11\n\n";
  
  }
}
    
  

// namespace Module
namespace {
  const Register r_WQM{
    "WQM", "SOS pnc WQM amplitude", &WQM};
  }// namespace

} // namespace Module

// namespace Module
// namespace Module