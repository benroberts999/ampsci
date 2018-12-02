#pragma once
#include "ElectronOrbitals.h"
#include <vector>

class HartreeFock{

public:

  HartreeFock(ElectronOrbitals &wf);
  //Later: can make some enums with options?



  std::vector<double>& v_abk(int a, int b, int k);
  std::vector<double>& v0_aa(int a) const;
  // const std::vector<double> &HartreeFock::v_abk(int a, int b, int k) const{}



private:

  void calculate_v_abk(const ElectronOrbitals &wf,
    int a, int b, int k, std::vector<double> & vabk);

  void form_Lambda_abk(const std::vector<int> &kappa);
  double Lambda_abk(int a, int b, int k) const;

  void startingApprox(ElectronOrbitals &wf);
  void HartreeFock::form_vdir(std::vector<double> &vdir);

  int m_ngp;
  int m_num_core_states;
  std::vector<int> l_list;
  std::vector<int> twoj_list;
  std::vector<int> kappa_index_list;
  std::vector<std::vector<std::vector<double> > > arr_Lambda_nmk;
  int kappa_from_index(int i);
  int index_from_kappa(int ka);
  int twoj_from_index(int i);

  //ElectronOrbitals &wf;

  std::vector<std::vector<double> > arr_v_bb0_r;
  std::vector<std::vector<std::vector<std::vector<double> > > > arr_v_abk_r;

  void formVexA(const ElectronOrbitals &wf, int a,
    std::vector<double> &vex_a);

  //does NO checks - can cause seg-fault
  double unsafe_Lambda_abk(int a, int b, int ik) const;


  void form_v_abk(const ElectronOrbitals &wf);



};
