#pragma once
#include <vector>

namespace ADAMS{

  const int AMO=8; //must be between 5 and 8 (for now). 7 Seems good.

  int solveDBS(std::vector<double> &f, std::vector<double> &g, double &en,
      const std::vector<double> &v, int n, int ka,
      const std::vector<double> &r, const std::vector<double> &drdt, double h,
      int &pinf, int &its, double &eps, double alpha, int log_dele_or=0);

  int findPracticalInfinity(double en,
    const std::vector<double> &v, const std::vector<double> &r, double alr);

  int findClassicalTurningPoint(double en, const std::vector<double> &v,
    int pinf);

  void trialDiracSolution(
    std::vector<double> &f, std::vector<double> &g, std::vector<double> &dg,
    double en, int ka, const std::vector<double> &v,
    const std::vector<double> &r, const std::vector<double> &drdt, double h,
    int ctp, int d_ctp, int pinf, double alpha);

  int countNodes(const std::vector<double> &f, int maxi=0);

  double largeEnergyChange(double &en,
    int &more, int &less, double &eupper, double &elower,
    double lfrac_de, bool more_nodes);

  double smallEnergyChangePT(double &anorm, double &en,
    const std::vector<double> &f, const std::vector<double> &g,
    const std::vector<double> &dg, int pinf, int ctp, int d_ctp, double alpha,
    const std::vector<double> &drdt, double h,
    int less, int more, double elower, double eupper);

  int outwardAM(std::vector<double> &p, std::vector<double> &q, double &en,
      const std::vector<double> &v, int ka,
      const std::vector<double> &r, const std::vector<double> &drdt, double h,
      int ctp, double alpha);

  int inwardAM(std::vector<double> &p, std::vector<double> &q, double &en,
      const std::vector<double> &v, int ka,
      const std::vector<double> &r, const std::vector<double> &drdt, double h,
      int ctp, int pinf, double alpha);

  int adamsMoulton(std::vector<double> &p, std::vector<double> &q, double &en,
      const std::vector<double> &v, int ka,
      const std::vector<double> &r, const std::vector<double> &drdt, double h,
      int ni, int nf, double alpha);

  void joinInOutSolutions(
    std::vector<double> &f, std::vector<double> &g, std::vector<double> &dg,
    const std::vector<double> &pin, const std::vector<double> &qin,
    const std::vector<double> &pout, const std::vector<double> &qout,
    int ctp, int d_ctp, int pinf);

  int getAdamsCoefs(std::vector<double> &mia, double &mid, double &miaa);

  int getOutwardCoefs(std::vector< std::vector<double> > &oie,
      std::vector<double> &oia, double &oid);

}
