#pragma once
#include "Coulomb/YkTable.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "MBPT/SpinorMatrix.hpp"
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
  std::pair<Basis, Basis> m_basis;
  Coulomb::YkTable m_Yeh;

  // Parameters of the sub-grid: initial/final points, stride
  std::size_t m_i0, m_stride, m_subgrid_points;

  int m_n_min_core;
  bool m_include_G;

  std::optional<HF::Breit> m_Br{};

public:
  Goldstone(const std::vector<DiracSpinor> &basis,
            const std::vector<DiracSpinor> &core, std::size_t i0,
            std::size_t stride, std::size_t size, int n_min_core = 1,
            bool include_G = false, const HF::Breit *Br = nullptr);

  //! Calculate Direct part of correlation potential
  GMatrix Sigma_direct(int kappa_v, double en_v,
                       const std::vector<double> &fk = {},
                       const std::vector<double> &etak = {},
                       int n_max_core = 99) const;

  // Calculate Exchange part of correlation potential
  GMatrix Sigma_exchange(int kappa_v, double en_v,
                         const std::vector<double> &fk = {}) const;

  // Calculate both parts of correlation potential.
  GMatrix Sigma_both(int kappa_v, double en_v,
                     const std::vector<double> &fk = {},
                     const std::vector<double> &etak = {},
                     int n_max_core = 99) const;

  // Calculates 2-body Breit correction to correlation potential. Must have Breit
  GMatrix dSigma_Breit2(int kappa_v, double en_v,
                        const std::vector<double> &fk = {},
                        const std::vector<double> &etak = {},
                        int n_max_core = 99, int m_max_n_breit = -1) const;

  const std::pair<Basis, Basis> &basis() const { return m_basis; }
  const Coulomb::YkTable &Yeh() const { return m_Yeh; }

  std::size_t stride() const { return m_stride; }
  int n_min() const { return m_n_min_core; }
  int lmax() const { return DiracSpinor::max_l(m_basis.second); }

private:
  inline double get_k(int k, const std::vector<double> &f) const {
    const auto sk = std::size_t(k);
    return sk < f.size() ? f[sk] : 1.0;
  }
};

} // namespace MBPT