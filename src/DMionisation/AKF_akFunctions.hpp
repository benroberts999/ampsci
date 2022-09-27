#pragma once
#include "DiracOperator/DiracOperator.hpp"
#include "LinAlg/Matrix.hpp"
#include <string>
#include <vector>
class Wavefunction;
class DiracSpinor;
class Grid;
namespace HF {
class HartreeFock;
}

namespace AKF {

//! Calculates K_a(dE,q) for a given core state, Fa, at specific energy deposition, dE. Sums over L and continuum states.
std::vector<double>
calculateK_nk(const Wavefunction &wf, const DiracSpinor &Fa, int max_L,
              double dE, const DiracOperator::jL &jl, bool subtract_1 = false,
              bool force_rescale = false, bool subtract_self = true,
              bool force_orthog = true, bool zeff_cont = false);

//! Heaviside step function
double Heaviside(double x);

std::vector<double> stepK_nk(const DiracSpinor &Fa, double dE,
                             const std::vector<double> &AFBE_table);

void writeToTextFile_AFBE(
    const std::string &fname,
    const std::vector<std::vector<std::vector<double>>> &AK,
    const std::vector<std::string> &nklst, const Grid &qgrid,
    const std::vector<double> deion);

/*! Writes K(dE,q) to text file (each core state sepparately)
@details
First two columns are values of dE and q.
Each subsequent  column is core state (1s, 2s, 2p ..., total) [final col is
total]. Each row is a different value of q, and each 'block' is a different dE
*/
void write_Knk_plaintext(
    const std::string &fname,
    const std::vector<std::vector<std::vector<double>>> &AK,
    const std::vector<std::string> &nklst, const Grid &qgrid,
    const Grid &Egrid);

// /*! Writes total K(dE,q) to text file in form: E/keV q/MeV K; each new E is new 'block'
void write_Ktot_plaintext(const std::string &fname,
                          const LinAlg::Matrix<double> &Keq, const Grid &Egrid,
                          const Grid &qgrid);

//! Reads/writes entire K(E,q) [and req'd details] to disk in binary form [read
//! in by dmeXsection program]
int akReadWrite(const std::string &fname, bool write,
                std::vector<std::vector<std::vector<double>>> &AK,
                std::vector<std::string> &nklst, double &qmin, double &qmax,
                double &dEmin, double &dEmax);

int akReadWrite_AFBE(const std::string &fname, bool write,
                     std::vector<std::vector<std::vector<double>>> &AK,
                     std::vector<std::string> &nklst, double &qmin,
                     double &qmax, std::vector<double> &deion);

} // namespace AKF
