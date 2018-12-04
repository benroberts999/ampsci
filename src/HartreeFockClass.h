#pragma once
#include "ElectronOrbitals.h"
#include <vector>
//Add change update git
class HartreeFock{



public:
  HartreeFock(ElectronOrbitals &wf, double eps_HF=1.e-8);
  void hartree_fock_core(ElectronOrbitals &wf, double eps_HF);

  double calculate_core_energy(const ElectronOrbitals &wf);


private:

  int MAX_HART_ITS = 64;

  void starting_approx_core(ElectronOrbitals &wf);

  std::vector< std::vector<double> > vex; //into class??

  void form_core_Lambda_abk(const std::vector<int> &kappa);
  double get_Lambda_abk(int a, int b, int k) const;
  int index_from_kappa(int ka) const;
  int twoj_from_index(int i) const;

  std::vector<std::vector<std::vector<double> > > m_core_Lambda_nmk;

  int m_ngp;
  int m_num_core_states;
  std::vector<int> twoj_list;
  std::vector<int> kappa_index_list;

  int kappa_from_index(int i); //XXX
  int l_from_index(int i) const;



  void calculate_v_abk(const ElectronOrbitals &wf, int a, int b, int k,
    std::vector<double> & vabk);

  void form_vbb0(const ElectronOrbitals &wf);
  void form_vabk_core(const ElectronOrbitals &wf);
  std::vector<double>& get_v_aa0(int a);
  std::vector<std::vector<double> >& get_v_abk(int a, int b);

  // std::vector<std::vector<double> > arr_v_bb0_r;
  std::vector<std::vector<std::vector<std::vector<double> > > > m_core_v_abk_r;

  void initialise_m_core_v_abk_r();

  void form_vdir(std::vector<double> &vdir, const ElectronOrbitals &wf,
    bool re_scale=false);
  void form_approx_vex_core(std::vector<std::vector<double> > &vex,
    const ElectronOrbitals &wf);
  void form_approx_vex_a(int a, std::vector<double> &vex_a,
    const ElectronOrbitals &wf);





};
