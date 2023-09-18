#pragma once
#include "CorrelationPotential.hpp" // for rgrid_params
#include "Coulomb/YkTable.hpp"
#include "HF/HartreeFock.hpp"
#include "MBPT/RDMatrix.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <memory>
#include <vector>

namespace MBPT {

//------------------------------------------------------------------------------
//! Class to construct Feynman diagrams, Green's functions and polarisation op.
class Goldstone {

  // Pointer to shared radial grid (full grid)
  std::shared_ptr<const Grid> m_grid;

  using Basis = std::vector<DiracSpinor>;
  const std::pair<Basis, Basis> m_basis;
  Coulomb::YkTable m_Yeh;

  // Parameters of the sub-grid: initial/final points, stride
  std::size_t m_i0, m_imax, m_stride, m_subgrid_points;

  bool m_include_G;
  // Lowest n to polarise in polarisation operator
  // int m_min_core_n;

  // maximum multipolarity, k
  // int m_max_k;

public:
  Goldstone(const std::vector<DiracSpinor> &basis, double e_Fermi,
            MBPT::rgrid_params subgrid, int n_min_core = 1,
            bool include_G = false);

  //! Calculate Direct part of correlation potential
  GMatrix Sigma_direct(int kappa_v, double en_v,
                       const std::vector<double> &fk = {},
                       const std::vector<double> &etak = {}) const;

  // nb: can be a little faster by combining w/ direct?
  GMatrix Sigma_exchange(int kappa_v, double en_v,
                         const std::vector<double> &fk = {}) const;

private:
  inline double get_k(int k, const std::vector<double> &f) const {
    const auto sk = std::size_t(k);
    return sk < f.size() ? f[sk] : 1.0;
  }
};

} // namespace MBPT