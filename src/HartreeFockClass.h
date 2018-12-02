#pragma once
#include "ElectronOrbitals.h"
#include <vector>
//Add change update git
class HartreeFock{



public:
  HartreeFock(ElectronOrbitals &wf);

private:
  void startingApprox(ElectronOrbitals &wf);

  void form_Lambda_abk(const std::vector<int> &kappa);
  double get_Lambda_abk(int a, int b, int k) const;
  int index_from_kappa(int ka) const;
  int twoj_from_index(int i) const;
  int get_num_ks(int a, int b) const;

  std::vector<std::vector<std::vector<double> > > arr_Lambda_nmk;

  int m_max_kappa_index_core;

  int m_ngp;
  int m_num_core_states;
  std::vector<int> twoj_list;
  std::vector<int> kappa_index_list;


  void calculate_v_abk(const ElectronOrbitals &wf, int a, int b, int k,
    std::vector<double> & vabk);

  void form_vbb0(const ElectronOrbitals &wf);
  void form_vabk(const ElectronOrbitals &wf);
  std::vector<double>& get_v_abk(int a, int b, int k);

  std::vector<std::vector<double> > arr_v_bb0_r;
  std::vector<std::vector<std::vector<std::vector<double> > > > arr_v_abk_r;

  void initialise_arr_v_bb0_r();
  void initialise_arr_v_abk_r();

  void form_vdir(std::vector<double> &vdir, const ElectronOrbitals &wf,
    bool re_scale=false);
  void form_approx_vex_core(std::vector<std::vector<double> > &vex,
    const ElectronOrbitals &wf);
  void form_approx_vex_a(int a, std::vector<double> &vex_a,
    const ElectronOrbitals &wf);














// public:

  // HartreeFock(ElectronOrbitals &wf);
  //Later: can make some enums with options?



  // std::vector<double>& v_abk(int a, int b, int k);
  // std::vector<double>& v0_aa(int a) const;
  // // const std::vector<double> &HartreeFock::v_abk(int a, int b, int k) const{}
  //


// private:
//
//   // void calculate_v_abk(const ElectronOrbitals &wf,
//   //   int a, int b, int k, std::vector<double> & vabk);
//
//   // void form_Lambda_abk(const std::vector<int> &kappa);
//   // double Lambda_abk(int a, int b, int k) const;
//
//
//   void HartreeFock::form_vdir(std::vector<double> &vdir);
//
//   // int m_ngp;
//   // int m_num_core_states;
//   // std::vector<int> l_list;
//   // std::vector<int> twoj_list;
//   // std::vector<int> kappa_index_list;
//
//   int kappa_from_index(int i);
//
//
//
//   //ElectronOrbitals &wf;
//
//
//
//   void formVexA(const ElectronOrbitals &wf, int a,
//     std::vector<double> &vex_a);
//
//   //does NO checks - can cause seg-fault
//   double unsafe_Lambda_abk(int a, int b, int ik) const;
//
//
//   void form_v_abk(const ElectronOrbitals &wf);



};
