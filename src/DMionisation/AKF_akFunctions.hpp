#pragma once
#include <string>
#include <vector>
class Wavefunction;
class DiracSpinor;

namespace AKF {

double CLkk(int L, int ka, int kb);

void writeToTextFile(const std::string &fname,
                     const std::vector<std::vector<std::vector<float>>> &AK,
                     const std::vector<std::string> &nklst, double qmin,
                     double qmax, double demin, double demax);

int akReadWrite(const std::string &fname, bool write,
                std::vector<std::vector<std::vector<float>>> &AK,
                std::vector<std::string> &nklst, double &qmin, double &qmax,
                double &dEmin, double &dEmax);

int calculateK_nk(const Wavefunction &wf, const DiracSpinor &psi, int max_L,
                  double dE,
                  const std::vector<std::vector<std::vector<double>>> &jLqr_f,
                  std::vector<float> &K_nk);

int calculateKpw_nk(const Wavefunction &wf, const DiracSpinor &psi, double dE,
                    const std::vector<std::vector<double>> &jl_qr,
                    std::vector<float> &K_nk);

std::vector<std::vector<std::vector<double>>>
sphericalBesselTable(int max_L, const std::vector<double> &q_array,
                     const std::vector<double> &r);

} // namespace AKF
