#pragma once
#include <string>
#include <vector>
class Wavefunction;
class DiracSpinor;

namespace AKF {

double CLkk(int L, int ka, int kb);

void writeToTextFile_AFBE(
    const std::string &fname,
    const std::vector<std::vector<std::vector<float>>> &AK,
    const std::vector<std::string> &nklst, double qmin, double qmax,
    const std::vector<double> deion);

void writeToTextFile(const std::string &fname,
                     const std::vector<std::vector<std::vector<float>>> &AK,
                     const std::vector<std::string> &nklst, double qmin,
                     double qmax, double demin, double demax);

int akReadWrite(const std::string &fname, bool write,
                std::vector<std::vector<std::vector<float>>> &AK,
                std::vector<std::string> &nklst, double &qmin, double &qmax,
                double &dEmin, double &dEmax);

int akReadWrite_AFBE(const std::string &fname, bool write,
                     std::vector<std::vector<std::vector<float>>> &AK,
                     std::vector<std::string> &nklst, double &qmin,
                     double &qmax, std::vector<double> &deion);

int calculateK_nk(const Wavefunction &wf, const DiracSpinor &psi, int max_L,
                  double dE,
                  std::vector<std::vector<std::vector<double>>> &jLqr_f,
                  std::vector<float> &K_nk, bool alt_akf, bool force_rescale,
                  bool subtract_self, bool force_orthog, std::string dmec,
                  double Zeff = -1);

int stepK_nk(const DiracSpinor &psi, double dE,
             const std::vector<float> &AFBE_table,
             std::vector<float> &AK_nk_step);

int calculateKpw_nk(const Wavefunction &wf, const DiracSpinor &psi, double dE,
                    std::vector<std::vector<double>> &jl_qr,
                    std::vector<float> &K_nk);

void sphericalBesselTable(std::vector<std::vector<std::vector<double>>> &jLqr_f,
                          int max_L, const std::vector<double> &q_array,
                          const std::vector<double> &r);

} // namespace AKF
