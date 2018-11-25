#pragma once
#include <vector>

namespace ADAMS{

  const int AMO=7; //must be between 5 and 8 (for now). 7 Seems good.

  int solveDBS(std::vector<double> &p, std::vector<double> &q, double &en,
      const std::vector<double> &v, double Z, int n, int ka,
      const std::vector<double> &r, const std::vector<double> &drdt, double h, int NGP,
      int &pinf, int &its, double &eps, double alpha, int log_dele_or=0);

  int outwardAM(std::vector<double> &p, std::vector<double> &q, double &en,
      const std::vector<double> &v, double Z, int ka,
      const std::vector<double> &r, const std::vector<double> &drdt, double h, int NGP,
      int ctp, double alpha);

  int inwardAM(std::vector<double> &p, std::vector<double> &q, double &en,
      const std::vector<double> &v, int ka,
      const std::vector<double> &r, const std::vector<double> &drdt, double h, int NGP,
      int ctp, int pinf, double alpha);

  int adamsMoulton(std::vector<double> &p, std::vector<double> &q, double &en,
      const std::vector<double> &v, int ka,
      const std::vector<double> &r, const std::vector<double> &drdt, double h, int ngp,
      int ni, int nf, double alpha);

  int getAdamsCoefs(std::vector<double> &mia, double &mid, double &miaa);

  int getOutwardCoefs(std::vector< std::vector<double> > &oie,
      std::vector<double> &oia, double &oid);

}
