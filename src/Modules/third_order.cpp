//Module for computing the third order hyperfine induced electric dipole tranisition matrix element required for general P,T odd operators
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
#include "Amplitudes/MatrixElements.hpp"
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


//Include option for reduced matrix elements scaling factor 
// 0 (default): standard reduced matrix element definition
// 1 (Johnson unconventional... debatable whether they even do it this way) but technically this is correct

std::vector<Coulomb::meTable<double>>
compute_me_3f(const DiracOperator::TensorOperator *hpnc,
              const DiracOperator::TensorOperator *he1,
              const DiracOperator::TensorOperator *const hfs,
              const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
              const DiracSpinor &v,const ExternalField::CorePolarisation *dV_PTodd=nullptr,
              const ExternalField::CorePolarisation *dV_e1=nullptr,
              const ExternalField::DiagramRPA *dV_hf=nullptr)
             {

  //Need matrix elements between valence states in transition and intermediate states from the spectrum

  //create tables for matrix elements
  Coulomb::meTable<double> pt_me;
  Coulomb::meTable<double> e1_me; 
  Coulomb::meTable<double> hf_me;  

  // initialise vector of tables for output
  std::vector<Coulomb::meTable<double>> ME_tables;

  //fill tables: need matrix elements
  //test with matrix element table existing (parallelised)
e1_me =  Amplitudes::me_table(spectrum,spectrum,he1,dV_e1);
hf_me =  Amplitudes::me_table(spectrum,spectrum,hfs,dV_hf);
pt_me =  Amplitudes::me_table(spectrum,spectrum,hpnc,dV_PTodd);

  //summing over all states in the spectrum
  //remove now after using me table

  /*

  // get angular momentum
  double tjw = w.twoj();



  //compute all matrix elements -> including RPA if dV object is passed based on input options
  for (const auto &i : spectrum) {
  // get angular momentum
  double tji = i.twoj();


  //specify pnc matrix element definition -- 

 
  double pnc_factor1;
  double pnc_factor2;
  if(pnc_type==0){
   pnc_factor1 = 1.0;
   pnc_factor2 = 1.0;
  }else{
    pnc_factor1 = 1.0/sqrt(tjw+1);
    pnc_factor2 =  1.0/sqrt(tji+1);
  }
    //first we consider matrix elements involving the initial or final states
    //compute me
    //pnc matrix elements--------------------------------------------------------------------
    const auto pnc_wi =
        hpnc-> reducedME(w, i) + (dV_pnc ? dV_pnc->dV(w, i) : 0.0);
    const auto pnc_iv =
        hpnc-> reducedME(i, v) + (dV_pnc ? dV_pnc->dV(i, v) : 0.0);
    const auto w_wf = w.shortSymbol();
    const auto i_wf = i.shortSymbol();
    const auto v_wf = v.shortSymbol();

    
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
  
    pnc_me.add(w, i, pnc_factor1*pnc_wi);
    pnc_me.add(i, v, pnc_factor2*pnc_iv);

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
    
    if (i.n() < 0) {
      std::cout << "Reduced hyperfine matrix element  " << " <<" << w_wf
                << "||h_hfs||" << i_wf << ">> = " << hfs_wi << "MHz \n.";
    }
    if (i.n() < 0) {
      std::cout << "Reduced hyperfine matrix element " << " <<" << i_wf
                << "||h_hfs||" << v_wf << ">> = " << hfs_iv << "MHz \n.";
    }
             

    //now include additional loop over spectrum
    for (const auto &j : spectrum) {
      double pnc_factor3;
      if(pnc_type==0){
      pnc_factor3 = 1.0;
      }else{
      pnc_factor3 = 1.0/sqrt(tji+1);
  }
/*
  if(j.en()<-37577 && i.en()<-37577 ){
    const auto pnc_ij = 0.0;
      const auto e1_ij = 0.0;
      const auto hfs_ij = 0.0;

      //add to table
      //hpnc->rme3js(i.twoj(), i.twoj();
      pnc_me.add(i, j, pnc_factor3*pnc_ij);
      e1_me.add(i, j, e1_ij);
      hf_me.add(i, j, hfs_ij);
  }else{
     
    const auto pnc_ij =
          hpnc->reducedME(i, j) + (dV_pnc ? dV_pnc->dV(i, j) : 0.0);
      const auto e1_ij = he1->reducedME(i, j) + (dV_e1 ? dV_e1->dV(i, j) : 0.0);
      const auto hfs_ij =
          hfs->reducedME(i, j) + (dV_hf ? dV_hf->dV(i, j) : 0.0);

      //add to table
      //hpnc->rme3js(i.twoj(), i.twoj();
      pnc_me.add(i, j, pnc_factor3*pnc_ij);
      e1_me.add(i, j, e1_ij);
      hf_me.add(i, j, hfs_ij);
  //}
  
    }
  }
   */
  //finally create vector of matrix element tables for single element output
  ME_tables.push_back(pnc_me2);
  ME_tables.push_back(e1_me2);
  ME_tables.push_back(hf_me2);
  return ME_tables;
}
            

// To compute the reduced matrix element we split the total quantity into three summations

// <w||h_[w+hf]|v>=h_1+h_2+h_3
//Inputs
// All matrix elements for each operator
// initial and final wavefunction
//spectrum
//nuclear spin
//Total initial and final angular momenta (F=I+J)
double h1(std::vector<Coulomb::meTable<double>> ME_tables,
          const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
          const DiracSpinor &v, int I2, int Fv2, int Fw2, int twok_hf,int twok_pt) {

//double sum over i and j in spectrum 
//sum_j -> no i dependence
//sum_ij
// h1=sum_j+sum_ij

const auto sum_j =0.0;
const auto sum_ij_1 =0.0;
const auto sum_ij_2 =0.0;





  //get angular momenta of valence electronic states from wavefunctions (we actually extract 2j) from kappa
  int kw = w.kappa();
  int kv = v.kappa();
  auto tjw = Angular::twoj_k(kw);
  auto tjv = Angular::twoj_k(kv);

  //get energies for valence states
  double e_v = v.en();
  double e_w = w.en();

  //get tables
  Coulomb::meTable<double> pt_me = ME_tables.at(0);
  Coulomb::meTable<double> e1_me = ME_tables.at(1);
  Coulomb::meTable<double> hf_me = ME_tables.at(2);

  for (const auto &j : spectrum) {
    // make sure j is not equal to v
    if (j != v) {
        double sum_j_cont = 0.0;
        //get six 6 symbols
        int kj = j.kappa();
        auto tjj = Angular::twoj_k(kj);
        double sixj1=Angular::sixj_2(tjj, tjv, 2, Fv2, Fw2, I2);
        double sixj2=Angular::sixj_2(tjw, tjw, twok_hf, I2, I2, Fw2);
        double sixj3=Angular::sixj_2(tjj, tjw, twok_pt, I2, I2, Fw2);
        const auto hf = pt_me.getv(w,j);
        const auto e1 = pt_me.getv(j,v);

        double e_j = j.en();
        double e_denom2 = 1.0/(e_w-e_j)^2;
        double sum_j_cont = sixj1*sixj2*sixj3*hf*e1*e_denom2;

        sum_j += sum_j_cont;
    

        for (const auto &i : spectrum) {
            if (i != w) {
                const auto sum_ij_1_cont =0.0;
                const auto sum_ij_2_cont=0.0;

                //get info for state i
                int ki = i.kappa();
                auto tji = Angular::twoj_k(ki);
                double e_i = i.en();
                double sixj4 = Angular::sixj_2(tjj, tjv, 2, Fv2, Fw2, I2);
                double sixj5 = Angular::sixj_2(tji, tjw, twok_hf, I2, I2, Fw2);
                double sixj6 = Angular::sixj_2(tjj, tjv, twok_pt, I2, I2, Fw2);
                //get matrix elements
                double hf2 = hf.getv(w,i);
                double pt2 = hf.getv(i,j);
                double e12 = hf.getv(j,v);
                double e_denom2 = (e_w-e_i)*(e_w-e_i);


                //second part 
                double sixj7 = Angular::sixj_2(tji, tjw, twok_pt, I2, I2, Fw2);
                double sixj8 = Angular::sixj_2(tjj, tjv, twok_hf, I2, I2, Fw2);
                 //get matrix elements
                 double hf3 = hf.getv(i,j);
                 double pt3 = hf.getv(w,i);
                 double e13 = hf.getv(j,v);
                 //get phase
                 double phase = Angular::neg1pow_2(tji-tjw);
                 //evaluate corresponding contribution

                 sum_ij_1_cont = phase*sixj4*sixj5*sixj6*hf2*pt2*e12/e_denom2;
                 sum_ij_2_cont = phase*sixj4*sixj7*sixj8*hf3*pt3*e13/e_denom2;

                 sum_ij_1 +=  sum_ij_1_cont;
                 sum_ij_2 +=  sum_ij_2_cont;


        }

      }
  }
}
const auto h1 = sum_ij_1+sum_ij_2+sum_j;
return h1;

}


//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
// h2 

double h2(std::vector<Coulomb::meTable<double>> ME_tables,
  const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
  const DiracSpinor &v, int I2, int Fv2, int Fw2, int twok_hf,int twok_pt) {

//double sum over i and j in spectrum 
//sum_j -> no i dependence
//sum_ij
// h1=sum_j+sum_ij

const auto sum_ij_1 =0.0;
const auto sum_ij_2 =0.0;





//get angular momenta of valence electronic states from wavefunctions (we actually extract 2j) from kappa
int kw = w.kappa();
int kv = v.kappa();
auto tjw = Angular::twoj_k(kw);
auto tjv = Angular::twoj_k(kv);

//get energies for valence states
double e_v = v.en();
double e_w = w.en();

//get tables
Coulomb::meTable<double> pt_me = ME_tables.at(0);
Coulomb::meTable<double> e1_me = ME_tables.at(1);
Coulomb::meTable<double> hf_me = ME_tables.at(2);

for (const auto &j : spectrum) {
// make sure j is not equal to v
if (j != v) {
  for (const auto &i : spectrum) {
    if (i != w) {
        const auto sum_ij_1_cont =0.0;
        const auto sum_ij_2_cont=0.0;

        //get info for state i and j
        int ki = i.kappa();
        auto tji = Angular::twoj_k(ki);
        double e_i = i.en();
        int kj = j.kappa();
        auto tjj = Angular::twoj_k(kj);
        double e_j = j.en();
        //six j symbols
        double sixj1 = Angular::sixj_2(tji, tjj, 2, Fv2, Fw2, I2);
        double sixj2 = Angular::sixj_2(tji, tjw, twok_hf, I2, I2, Fw2);
        double sixj3 = Angular::sixj_2(tjv, tji, twok_pt, I2, I2, Fv2);
        double sixj4 = Angular::sixj_2(tji, tjw, twok_pt, I2, I2, Fw2);
        double sixj5 = Angular::sixj_2(tjv, tjj, twok_hf, I2, I2, Fv2);

        //get matrix elements
        //get matrix elements
        double hf = hf.getv(w,i);
        double pt = hf.getv(j,v);
        double e1 = hf.getv(i,j);

        double hf2 = hf.getv(j,v);
        double pt2 = hf.getv(w,i);
        
        //energy denominator
        double e_denom = (e_w-e_j)*(e_v-e_i);

         //evaluate corresponding contributions



         sum_ij_1_cont = sixj1*sixj2*sixj3*hf*pt*e1/e_denom;
         sum_ij_2_cont = sixj1*sixj4*sixj5*hf2*pt2*e12/e_denom;

         sum_ij_1 +=  sum_ij_1_cont;
         sum_ij_2 +=  sum_ij_2_cont;


}
}

}
}
}
const auto h2_phase_factor = Angular::neg1pow_2(Fv2-Fw2+tjw-tjv);
const auto h2 = h2_phase_factor*(sum_ij_1+sum_ij_2);
return h2;
}

//------------------------------------------------------------------------------------------------------
// calc h3 (third part)
std::vector<double> h1(std::vector<Coulomb::meTable<double>> ME_tables,
  const std::vector<DiracSpinor> &spectrum, const DiracSpinor &w,
  const DiracSpinor &v, int I2, int Fv2, int Fw2, int twok_hf,int twok_pt) {

//double sum over i and j in spectrum 
//sum_j -> no i dependence
//sum_ij
// h1=sum_j+sum_ij

const auto sum_j =0.0;
const auto sum_ij_1 =0.0;
const auto sum_ij_2 =0.0;





//get angular momenta of valence electronic states from wavefunctions (we actually extract 2j) from kappa
int kw = w.kappa();
int kv = v.kappa();
auto tjw = Angular::twoj_k(kw);
auto tjv = Angular::twoj_k(kv);

//get energies for valence states
double e_v = v.en();
double e_w = w.en();

//get tables
Coulomb::meTable<double> pt_me = ME_tables.at(0);
Coulomb::meTable<double> e1_me = ME_tables.at(1);
Coulomb::meTable<double> hf_me = ME_tables.at(2);

//get the diagonal valence hyperfine matrix element
const auto hf = hf_me.getv(v,v);

for (const auto &j : spectrum) {

 // make sure j is not equal to vx
 if (j != w) {
 double sum_j_cont = 0.0;
 //get six 6 symbols
 int kj = j.kappa();
 auto tjj = Angular::twoj_k(kj);
 double sixj1=Angular::sixj_2(tjw, tjj, 2, Fv2, Fw2, I2);
 double sixj2=Angular::sixj_2(tjv, tjv, twok_hf, I2, I2, Fv2);
 double sixj3=Angular::sixj_2(tjv, tjj, twok_pt, I2, I2, Fv2);
 const auto e1 = e1_me.getv(w,j);
 const auto pt = pt_me.getv(j,v);

 double e_j = j.en();
 double e_denom2 = (e_v-e_j)^2;
 double sum_j_cont = sixj1*sixj2*sixj3*hf*e1/e_denom2;

 sum_j += sum_j_cont;
}

for (const auto &i : spectrum) {
    if (i != v) {
        const auto sum_ij_1_cont =0.0;
        const auto sum_ij_2_cont=0.0;

        //get info for state i
        int ki = i.kappa();
        auto tji = Angular::twoj_k(ki);
        double e_i = i.en();
        double sixj4 = Angular::sixj_2(tjw, tjj, 2, Fv2, Fw2, I2);
        double sixj5 = Angular::sixj_2(tji, tjj, twok_hf, I2, I2, Fv2);
        double sixj6 = Angular::sixj_2(tjv, tji, twok_pt, I2, I2, Fv2);
        //get matrix elements
        double hf2 = hf.getv(j,i);
        double pt2 = hf.getv(i,v);
        double e12 = hf.getv(w,j);
        double e_denom = (e_v-e_i)*(e_v-e_j);
        //second part 
        double sixj7 = Angular::sixj_2(tji, tjj, twok_pt, I2, I2, Fv2);
        double sixj8 = Angular::sixj_2(tjv, tji, twok_hf, I2, I2, Fv2);
         //get matrix elements
         double hf3 = hf.getv(i,v);
         double pt3 = hf.getv(j,i);
         double e13 = hf.getv(w,j);
         //get phase
         double phase = Angular::neg1pow_2(tji-tjv);
         //evaluate corresponding contribution

         sum_ij_1_cont = phase*sixj4*sixj5*sixj6*hf2*pt2*e12/e_denom;
         sum_ij_2_cont = phase*sixj4*sixj7*sixj8*hf3*pt3*e13/e_denom;

         sum_ij_1 +=  sum_ij_1_cont;
         sum_ij_2 +=  sum_ij_2_cont;


}

}
}

const auto h2 = sum_ij_1+sum_ij_2+sum_j;
return h3;
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
               {"HF_k", " hyperfine operator rank"},
               {"PT_k", "PT odd operator rank"}});
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
  //for now just pnc operator --> need to generalise to take operator as input
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
    int HF_2k = 2*input.get("HF_k", 0);
    int PT_2k = 2*input.get("PT_k", 0);

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
    const auto pnc_type = input.get("pnc_type", 0);
    std::cout << "Computing and storing all matrix elements.... ";
    std::vector<Coulomb::meTable<double>> Table = compute_me_3f(
        &hpnc, &he1, hfs.get(), spectrum, Fw, Fv,pnc_type, &dVpnc, &dVE1, rpa_hf.get());
    std::cout << "Complete. \n";

    //use matrix elements to compute the amplitude
    int I2 = input.get("two_I", 0);
    const auto Fv2 = input.get("two_Fv", 0);
    const auto Fw2 = input.get("two_Fw", 0);
    //const auto two_k = input.get("two_k", 0);


    //now compute terms in total third order amplitude summation

    double h_1 = h1(Table, spectrum, Fw, Fv, I2, Fv2, Fw2,HF_2k,PT_2k );

    double h_2 = h2(Table, spectrum, Fw, Fv, I2, Fv2, Fw2, HF_2k, PT_2k);

    double h_3 = h3(Table, spectrum, Fw, Fv, I2, Fv2, Fw2, HF_2k, PT_2k);

    double I = 0.5 * I2;

    //return final values 
    //first negative and positive contributiions separately

    //pre factors depend on respective operators 

    double pre_factor_const =Angular::neg1pow_2(Fv2+I2+tjw+2)*sqrt((Fw2+1)*(Fv2+1)) ;

    if(HF_k==1){
      pre_factor_HFoperator = sqrt(I * (I + 1) * (I2 + 1));
      //magnetic dipole hyperfine

    }else{
      pre_factor_HFoperator  = sqrt((I2+3)*(I2+1)*(I+1)/(I*(I-1)));
      //electric quadrupole hyperfine

    }

    //pre factor for pnc (SI)
    pre_factor_PToperator = sqrt(2*I+1);



    double h_sum = pre_factor_const*pre_factor_HFoperator*pre_factor_PToperator * (h_1+h_2+h_3);
   

    //need additional factor to account for different definition of the reduced matrix element (phase+3j symbol)
    // double tjw = Fw.twoj();
    std::cout << "\n Total value is " << h_sum << "  MHz .with RPA\n\n";
    std::cout << "In atomic unuts value is " << h_sum / PhysConst::Hartree_MHz
              << "  au .";

  } else {
    // Compute and store matrix elements without RPA
    const auto pnc_type = input.get("pnc_type", 0);
    std::cout << "Computing and storing all matrix elements without RPA.... ";
    std::vector<Coulomb::meTable<double>> Table =
        compute_me_3f(&hpnc, &he1, hfs.get(), spectrum, Fw, Fv,pnc_type);
    std::cout << "Complete. \n";

    //use matrix elements to compute the amplitude
    int I2 = input.get("two_I", 0);
    const auto Fv2 = input.get("two_Fv", 0);
    const auto Fw2 = input.get("two_Fw", 0);
    int HF_2k = 2*input.get("HF_k", 0);
    int PT_2k = 2*input.get("PT_k", 0);

    std::vector<double> h_1 = h1(Table, spectrum, Fw, Fv, I2, Fv2, Fw2, HF_2k, PT_2k);

    std::vector<double> h_2 = h2(Table, spectrum, Fw, Fv, I2, Fv2, Fw2, HF_2k, PT_2k);

    double I = 0.5 * I2;

    double pre_factor_const =Angular::neg1pow_2(Fv2+I2+tjw+2)*sqrt((Fw2+1)*(Fv2+1)) ;

    if(HF_k==1){
      pre_factor_HFoperator = sqrt(I * (I + 1) * (I2 + 1));
      //magnetic dipole hyperfine

    }else{
      pre_factor_HFoperator  = sqrt((I2+3)*(I2+1)*(I+1)/(I*(I-1)));
      //electric quadrupole hyperfine

    }

    //pre factor for pnc (SI)
    pre_factor_PToperator = sqrt(2*I+1);



    double h_sum = pre_factor_const*pre_factor_HFoperator*pre_factor_PToperator * (h_1+h_2+h_3);
   

    //need additional factor to account for different definition of the reduced matrix element (phase+3j symbol)
    // double tjw = Fw.twoj();
    std::cout << "\n Total value is " << h_sum << "  MHz .without RPA\n\n";
    std::cout << "In atomic unuts value is " << h_sum / PhysConst::Hartree_MHz
              << "  au .";
}
    
  

// namespace Module
namespace {
  const Register r_hf_pert_weak{
    "third_order", "Third order SOS dipole amplitude incl weak x hyperfine", &third_order};
  }// namespace

} // namespace Module

// namespace Module
// namespace Module