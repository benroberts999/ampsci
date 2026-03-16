#pragma once
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Angular/Wigner369j.hpp"

/*!
@brief
Angular provides functions and classes for calculating and storing angular
factors (3,6,9-J symbols etc.)

@details
Provides functions to:
 - Calculate wigner 3,6,9-J symbols + Clebsh-Gordon coefs etc..
 - Also provied look-up tables, which are faster than calculating symbols
on-the-fly
 - Wrapper functions to calculate wigner 3,6,9-J symbols.
 - Uses GSL:
https://www.gnu.org/software/gsl/doc/html/specfunc.html?highlight=3j#coupling-coefficients

The general equations are defined:
\f{align}{
    \frac{1}{r_{12}} &= \sum_{kq} \frac{r_<^k}{r_>^{k+1}}(-1)^q
                        C^k_{-q}(\hat{n}_1)C^k_{q}(\hat{n}_2)\\
    C^k_{q} &\equiv \sqrt{\frac{4\pi}{2k+1}} Y_{kq}(\hat{n}),
\f}
and
\f{align}{
  C^k_{ab} &\equiv \langle{\kappa_a}||C^k||{\kappa_b}\rangle 
            \equiv  (-1)^{j_a+1/2} \widetilde C^k_{ab} ,\\
            &=  (-1)^{j_a+1/2}\sqrt{[j_a][j_b]}
            \begin{pmatrix}
              {j_a} & {j_b} & {k} \\ {-1/2} & {1/2} &{0}
            \end{pmatrix}
            \pi(l_a+l_b+k). 
\f}
*/
namespace Angular {}