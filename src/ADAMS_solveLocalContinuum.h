#ifndef _SOLVECNTM_H
#define _SOLVECNTM_H
#include "ADAMS_solveLocalBS.h"
#include <fstream> //XXX remove after!
namespace ADAMS{

double fitQuadratic(double x1, double x2, double x3, double y1, double y2,
  double y3);

int solveContinuum(std::vector<double> &p, std::vector<double> &q, double en,
    std::vector<double> v, double Z, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGPb, int NGPc, int i_asym,
    double alpha);

double findSineAmplitude(std::vector<double> pc, std::vector<double> rc, int NGPc, int i_asym);

int findAsymptoticRegion(std::vector<double> pc, std::vector<double> rc,
  int NGPb, int NGPc, int i_asym);

}
#endif
