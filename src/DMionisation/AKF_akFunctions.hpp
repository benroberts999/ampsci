#pragma once
#include <string>
#include <vector>
class Wavefunction;

namespace AKF {

double CLkk(int L, int ka, int kb);

void writeToTextFile(std::string fname,
                     const std::vector<std::vector<std::vector<float>>> &AK,
                     const std::vector<std::string> &nklst, double qmin,
                     double qmax, double demin, double demax);

int akReadWrite(std::string fname, bool write,
                std::vector<std::vector<std::vector<float>>> &AK,
                std::vector<std::string> &nklst, double &qmin, double &qmax,
                double &dEmin, double &dEmax);

int calculateK_nk(const Wavefunction &wf, std::size_t nk, int max_L,
                  double dE,
                  std::vector<std::vector<std::vector<double>>> &jLqr_f,
                  std::vector<float> &K_nk, double Zeff = -1);

int calculateKpw_nk(const Wavefunction &wf, std::size_t nk, double dE,
                    std::vector<std::vector<double>> &jl_qr,
                    std::vector<float> &K_nk);

void sphericalBesselTable(std::vector<std::vector<std::vector<double>>> &jLqr_f,
                          int max_L, const std::vector<double> &q_array,
                          const std::vector<double> &r);

} // namespace AKF
