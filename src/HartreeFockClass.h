#pragma once
#include "ElectronOrbitals.h"
#include <vector>

class HartreeFock{

public:

  HartreeFock(ElectronOrbitals &wf);
  //Later: can make some enums with options?

  double Lambda_abk(int a, int b, int k) const;

  std::vector<double>& v_abk(int a, int b, int k) const;
  std::vector<double>& v0_aa(int a) const;
  // const std::vector<double> &HartreeFock::v_abk(int a, int b, int k) const{}



private:
  std::vector<int> l_list;
  std::vector<int> twoj_list;
  std::vector<int> kappa_index_list;
  std::vector<std::vector<std::vector<double> > > arr_Lambda_nmk;
  std::vector<std::vector<std::vector<std::vector<double> > > > arr_v_abk_r;

  int kappa_from_index(int i);
  int index_from_kappa(int ka);
  int twoj_from_index(int i);

  //does NO checks - can cause seg-fault
  double unsafe_Lambda_abk(int a, int b, int k, int kmin) const;

  void startingApprox();
  void form_v_abk();
  void form_Lambda_abk();


};
