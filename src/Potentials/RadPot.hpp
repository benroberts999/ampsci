#pragma once
#include "FGRadPot.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/InputBlock.hpp"
#include "Maths/Interpolator.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <vector>

//! Radiative QED corrections (Flambaum-Ginges Radiative Potenti)
namespace QED {

//==============================================================================
//! Constructs and stores the Flambaum-Ginges QED Radiative Potential
class RadPot {

public:
  //============================================================================
  //! Scale factors for Uehling, high, low, magnetic, Wickman-Kroll
  struct Scale {
    double u, h, l, m, wk;
  };

  //============================================================================
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
    //! The factors as given (one per l; the last applies to higher l)
    const std::vector<double> &values() const { return m_x; }

    double operator()(int l) const;
  };

  //============================================================================
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
         bool tprint = true, bool do_readwrite = true,
         const std::string &label = "");

  //! The constructor arguments other than the grid: enough to rebuild the
  //! potential (on the same grid) with the Params constructor
  struct Params {
    double Z;
    double rN;
    double rcut;
    Scale f;
    std::vector<double> xl;
  };
  //! As params() of an existing potential
  Params params() const { return {m_Z, m_rN, m_rcut, m_f, m_xl.values()}; }
  //! Constructs from params(); see the main constructor
  RadPot(const std::vector<double> &r, const Params &params, bool tprint = true,
         bool do_readwrite = true, const std::string &label = "")
    : RadPot(r, params.Z, params.rN, params.rcut, params.f, Xl(params.xl),
             tprint, do_readwrite, label) {}

  bool read_write(const std::vector<double> &r, IO::FRW::RoW rw,
                  const std::string &label_x = "");

  void form_potentials(const std::vector<double> &r);

  //! Returns entire electric part of potential
  std::vector<double> Vel(int l = 0) const;
  //! Returns H_mag (magnetic self-energy form vactor)
  std::vector<double> Hmag(int) const;
  //! Uehling potential
  std::vector<double> Vu(int l = 0) const;
  //! Low-frequency electric self-energy potential
  std::vector<double> Vl(int l = 0) const;
  //! High-frequency electric self-energy potential
  std::vector<double> Vh(int l = 0) const;

  template <typename Func>
  std::vector<double> fill(Func f, const std::vector<double> &r);

  template <typename Func>
  std::vector<double> fill(Func f, const std::vector<double> &r,
                           std::size_t stride);
};

//============================================================================
//==============================================================================
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

//==============================================================================
template <typename Func>
std::vector<double> RadPot::fill(Func f, const std::vector<double> &r,
                                 std::size_t stride) {

  const auto rcut = m_rcut == 0.0 ? r.back() : m_rcut;

  // index for r cut-off
  const auto icut =
    std::size_t(std::distance(
      begin(r),
      std::find_if(begin(r), end(r), [rcut](auto ri) { return ri > rcut; }))) /
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

//==============================================================================
//! Function constructs a Radiative potential with given input parameters; rN_au is nuclear radius (not rms radius), in atomic units
RadPot ConstructRadPot(const std::vector<double> &r, double Z_eff, double rN_au,
                       const IO::InputBlock &input = {}, bool print = true,
                       bool do_readwrite = true, const std::string &label = "");

} // namespace QED
