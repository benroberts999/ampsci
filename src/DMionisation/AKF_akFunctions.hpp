#pragma once
#include <string>
#include <vector>
class Wavefunction;
class DiracSpinor;
class Grid;

namespace AKF {

//! Angular factor
double CLkk(int L, int ka, int kb);

/*! Writes K(dE,q) to text file (each core state sepparately)
@details
First two columns are values of dE and q.
Each subsequent  column is core state (1s, 2s, 2p ..., total) [final col is
total]. Each row is a different value of q, and each 'block' is a different dE
*/
void write_Knk_plaintext(const std::string &fname,
                         const std::vector<std::vector<std::vector<float>>> &AK,
                         const std::vector<std::string> &nklst,
                         const Grid &qgrid, const Grid &Egrid);

/*! Writes total K(dE,q) to text file
@details
Each q is a column, each dE is a row. First row/col is dE/q values
*/
void write_Ktot_plaintext(
    const std::string &fname,
    const std::vector<std::vector<std::vector<float>>> &AK, const Grid &qgrid,
    const Grid &Egrid);

//! Reads/writes entire K(E,q) [and req'd details] to disk in binary form [read
//! in by dmeXsection program]
int akReadWrite(const std::string &fname, bool write,
                std::vector<std::vector<std::vector<float>>> &AK,
                std::vector<std::string> &nklst, double &qmin, double &qmax,
                double &dEmin, double &dEmax);

//! Calculates K(q) for a given core state, psi, at specific energy deposition,
//! dE. Sums over L and continuum states.
std::vector<float>
calculateK_nk(const Wavefunction &wf, const DiracSpinor &psi, int max_L,
              double dE,
              const std::vector<std::vector<std::vector<double>>> &jLqr_f);

//! Calculates K(q) for a given core state, psi, at specific energy deposition,
//! dE, using plane-wave approximation. Note: only for testing; plane wave
//! approx. is not reasonable. The formula used here is an approximate one, that
//! applies only for large q.
std::vector<float>
calculateKpw_nk(const Wavefunction &wf, const DiracSpinor &psi, double dE,
                const std::vector<std::vector<double>> &jl_qr);

//! Given maximim L, array of q (mom. transfer), and r, returns lookup table of
//! sphericalBessel functions, J_L(qr). In form of vector<vector<vector>>,
//! j[L][q][r] - this makes it easiest to pass as function of r to radial
//! integration functions
std::vector<std::vector<std::vector<double>>>
sphericalBesselTable(int max_L, const std::vector<double> &q_array,
                     const std::vector<double> &r);

} // namespace AKF
