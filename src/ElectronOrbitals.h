#pragma once
#include <vector>
#include <string>

class ElectronOrbitals{

  public:

    ElectronOrbitals(int in_z, int in_a=0, int in_ngp=2048,
        double rmin=1e-6, double rmax=250., double var_alpha=1);
    ElectronOrbitals(std::string s_in_z, int in_a=0, int in_ngp=2048,
        double rmin=1e-6, double rmax=250., double var_alpha=1);

    void clearAll();

  public:

    //orbitals:
    std::vector< std::vector<double> > f;
    std::vector< std::vector<double> > g;
    std::vector<double> en;
    //state info:
    std::vector<int> nlist;
    std::vector<int> kappa;
    //info from solveing DE
    std::vector<int> pinflist; //practical infinity
    std::vector<int> itslist; //num. iterations
    std::vector<double> epslist; //convergance
    //occupancy fraction.
    //Note: avg over non-rel for core, but rel for valence!
    std::vector<double> occ_frac;

    int num_core_states;
    int num_core_electrons; // Nc = N - M

    //grid
    std::vector<double> r;
    std::vector<double> drdt;
    std::vector<double> dror;
    int ngp;
    double h;

    //Potentials
    std::vector<double> vnuc;
    std::vector<double> vdir; //direct/local part of the electron potential

  public:

    double get_alpha() const;
    int Z() const;
    int A() const;

    std::string seTermSymbol(int ink);

    int solveLocalDirac(int n, int k, double en_a=0, int log_dele_or=0);
    int reSolveDirac(int i, double e_a=0, int log_dele_or=0);
    int reSolveDirac(int i, double e_a, const std::vector<double> &vex,
      int log_dele_or=0);

    double diracen(double z, double n, int k);
    int sphericalNucleus(double rnuc=0);
    int fermiNucleus(double t=0, double c=0);

    int getRadialIndex(double r_target);

    int solveInitialCore(std::string str_core_in, int log_dele_or=0);
    bool isInCore(int n, int k);
    int maxCore_n(void);

    void orthonormaliseOrbitals(int num_its=1);
    void orthonormaliseValence(int iv, int num_its);

    int sortedEnergyList(std::vector<int> &sort_list);

    //Single function to form grid, takes in option (w/ enums)!
    //Default: logLinear..
    int logLinearRadialGrid(int ngp_in, double r0, double rmax, double b=4.);
    int logLinearRadialGrid(double in_h, double r0, double rmax, double b=4.);

  private:

    //store internal value for alpha (allows variation)
    double alpha;

    //Atom info:
    int Z_,A_;
    //std::string atom;

    //number of electrons in each core shell (non-rel??)
    std::vector<int> num_core_shell;

  private:

    int determineCore(std::string str_core_in);

    //Grid
    int exponentialRadialGrid(int ngp_in, double r0, double rmax);
    int zeroNucleus();

    double enGuessCore(int n, int l);
    double enGuessVal(int n, int ka);

};
