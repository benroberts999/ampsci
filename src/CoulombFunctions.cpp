#include "CoulombFunctions.hpp"
#include "DiracSpinor.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include <vector>

// namespace Coulomb

void Coulomb::calculate_v_abk(const DiracSpinor &phi_a,
                              const DiracSpinor &phi_b, const int k,
                              std::vector<double> &vabk)
// Calculalates v^k_ab screening function.
// Note: should only call for a>=b, and for k's with non-zero angular coefs
// (nothing bad will happen othersie, but no point!)
// Since v_ab = v_ba
//
// Stores in vabk (reference to whatever) - must already be sized corectly!
//
// r_min := min(r,r')
// rho(r') := fa(r')*fb(r') + ga(r')gb(r')
// v^k_ab(r) = Int_0^inf [r_min^k/r_max^(k+1)]*rho(f') dr'
//           = Int_0^r [r'^k/r^(k+1)]*rho(r') dr'
//             + Int_r^inf [r^k/r'^(k+1)]*rho(r') dr'
//          := A(r)/r^(k+1) + B(r)*r^k
// A(r0)  = 0
// B(r0)  = Int_0^inf [r^k/r'^(k+1)]*rho(r') dr'
// A(r_n) = A(r_{n-1}) + (rho(r_{n-1})*r_{n-1}^k)*dr
// B(r_n) = A(r_{n-1}) + (rho(r_{n-1})/r_{n-1}^(k+1))*dr
// v^k_ab(rn) = A(rn)/rn^(k+1) + B(rn)*rn^k
{
  auto irmax = std::min(phi_a.pinf, phi_b.pinf);
  auto du = phi_a.p_rgrid->du;

  double Ax = 0, Bx = 0;
  for (std::size_t i = 0; i < irmax; i++) {
    Bx += phi_a.p_rgrid->drdu[i] *
          (phi_a.f[i] * phi_b.f[i] + phi_a.g[i] * phi_b.g[i]) /
          pow(phi_a.p_rgrid->r[i], k + 1);
  }

  // For "direct" part, can't cut!
  if (phi_a == phi_b)
    irmax = phi_a.p_rgrid->ngp;

  vabk[0] = Bx * du;
  for (std::size_t i = 1; i < irmax; i++) {
    auto Fdr = phi_a.p_rgrid->drdu[i - 1] * (phi_a.f[i - 1] * phi_b.f[i - 1] +
                                             phi_a.g[i - 1] * phi_b.g[i - 1]);
    Ax = Ax + Fdr * pow(phi_a.p_rgrid->r[i - 1], k);
    Bx = Bx - Fdr / pow(phi_a.p_rgrid->r[i - 1], k + 1);
    vabk[i] = du * (Ax / pow(phi_a.p_rgrid->r[i], k + 1) +
                    Bx * pow(phi_a.p_rgrid->r[i], k));
  }
  for (std::size_t i = irmax; i < phi_a.p_rgrid->ngp; i++) {
    vabk[i] = 0; // maybe not needed?
  }
}

// } // namespace Coulomb
