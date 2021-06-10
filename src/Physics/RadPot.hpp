#pragma once
#include "FGRadPot.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Maths/Interpolator.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <vector>

namespace QED {

//! Class holds Flambaum-Ginges QED Radiative Potential
class RadPot {

public:
  //****************************************************************************
  //! Scale factors for Uehling, high, low, magnetic, Wickman-Kroll
  struct Scale {
    double u, h, l, m, wk;
  };

  //****************************************************************************
  //! Extra fitting for s,p,d etc. states.
  /*! @details Will assume factor for higher l's same as last one.
      e.g., {1,1} => 1,1,1,1,...
      e.g., {1,0} => 1,0,0,0,...
  */
  struct Xl {

  private:
    std::vector<double> m_x;

  public:
    Xl();
    Xl(std::vector<double> x);

    double operator()(int l) const;
  };

  //****************************************************************************
private:
  double m_Z, m_rN, m_rcut;
  Scale m_f;
  Xl m_xl;
  bool print;

  std::vector<double> mVu{};  // Uehling
  std::vector<double> mVh{};  // High-freq electric SE (without A)
  std::vector<double> mVl{};  // Low-freq electric SE (without B)
  std::vector<double> mHm{};  // Magnetic FF
  std::vector<double> mVwk{}; // Approx Wickman-Kroll

public:
  //! Empty constructor
  RadPot();

  //! Constructor: will build potential
  /*! @details
    rcut is maxum radius (atomic units) to calc potential for.
  */
  RadPot(const std::vector<double> &r, double Z, double rN = 0.0,
         double rcut = 0.0, Scale f = {1.0, 1.0, 1.0, 1.0, 0.0}, Xl xl = {},
         bool tprint = true, bool do_readwrite = true);

  bool read_write(const std::vector<double> &r, IO::FRW::RoW rw);

  void form_potentials(const std::vector<double> &r);

  //! Returns entire electric part of potential
  std::vector<double> Vel(int l = 0) const;
  //! Returns H_mag
  std::vector<double> Hmag(int) const;
  std::vector<double> Vu(int l = 0) const;
  std::vector<double> Vl(int l = 0) const;
  std::vector<double> Vh(int l = 0) const;

  template <typename Func>
  std::vector<double> fill(Func f, const std::vector<double> &r);

  template <typename Func>
  std::vector<double> fill(Func f, const std::vector<double> &r,
                           std::size_t stride);
};

//****************************************************************************
//******************************************************************************
template <typename Func>
std::vector<double> RadPot::fill(Func f, const std::vector<double> &r) {
  std::vector<double> v;
  v.resize(r.size());

  const auto rcut = m_rcut == 0.0 ? r.back() : m_rcut;

  // index for r cut-off
  const auto icut = std::size_t(std::distance(
      begin(r),
      std::find_if(begin(r), end(r), [rcut](auto ri) { return ri > rcut; })));

#pragma omp parallel for
  for (auto i = 0ul; i < icut; ++i) {
    // nb: Use H -> H+V (instead of H-> H-V), so change sign!
    v[i] = -f(m_Z, r[i], m_rN);
  }
  return v;
}

//******************************************************************************
template <typename Func>
std::vector<double> RadPot::fill(Func f, const std::vector<double> &r,
                                 std::size_t stride) {

  const auto rcut = m_rcut == 0.0 ? r.back() : m_rcut;

  // index for r cut-off
  const auto icut =
      std::size_t(std::distance(
          begin(r), std::find_if(begin(r), end(r),
                                 [rcut](auto ri) { return ri > rcut; }))) /
      stride;

  std::vector<double> tv, tr;
  tv.resize(icut);
  tr.resize(icut);

#pragma omp parallel for
  for (auto i = 0ul; i < icut; ++i) {
    // nb: Use H -> H+V (instead of H-> H-V), so change sign!
    tv[i] = -f(m_Z, r[i * stride], m_rN);
    tr[i] = r[i * stride];
  }
  return stride == 1 ? tv : Interpolator::interpolate(tr, tv, r);
}

} // namespace QED
