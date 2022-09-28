#include "RadPot.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Maths/Interpolator.hpp"
#include "Physics/FGRadPot.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <vector>

namespace QED {

//==============================================================================
RadPot::Xl::Xl() : m_x{1.0} {}
RadPot::Xl::Xl(std::vector<double> x) : m_x{std::move(x)} {}

double RadPot::Xl::operator()(int l) const {
  if (m_x.empty())
    return 1.0;
  if (l >= (int)m_x.size())
    return m_x.back();
  return m_x.at(std::size_t(l));
}

//==============================================================================
RadPot::RadPot()
    : m_Z(0),
      m_rN(0.0),
      m_rcut(0.0),
      m_f({0.0, 0.0, 0.0, 0.0, 0.0}),
      m_xl{},
      print(true) {}

//==============================================================================
RadPot::RadPot(const std::vector<double> &r, double Z, double rN, double rcut,
               Scale f, Xl xl, bool tprint, bool do_readwrite)
    : m_Z(Z), m_rN(rN), m_rcut(rcut), m_f(f), m_xl(xl), print(tprint) {

  bool read_ok = false;
  if (do_readwrite)
    read_ok = read_write(r, IO::FRW::RoW::read);
  if (!read_ok) {
    form_potentials(r);
    if (do_readwrite)
      read_write(r, IO::FRW::RoW::write);
  }
}

//==============================================================================
bool RadPot::read_write(const std::vector<double> &r, IO::FRW::RoW rw) {

  std::string label = "_";
  if (m_f.u != 0.0)
    label += 'u';
  if (m_f.h != 0.0)
    label += 'h';
  if (m_f.l != 0.0)
    label += 'l';
  if (m_f.m != 0.0)
    label += 'm';
  if (m_f.wk != 0.0)
    label += "w";
  if (m_rN == 0.0)
    label += "_pt";
  const auto fname = std::to_string(int(m_Z + 0.1)) + label + ".qed";

  const auto readQ = rw == IO::FRW::read;

  if (readQ && !IO::FRW::file_exists(fname))
    return false;

  const auto rw_str = !readQ ? "Writing to " : "Reading from ";
  std::cout << rw_str << "QED rad. pot. file: " << fname << "\n" << std::flush;

  std::fstream iofs;
  IO::FRW::open_binary(iofs, fname, rw);

  // read/write info + Grid:
  std::vector<double> t_r;
  if (!readQ) {
    t_r = r; // copy; cannot write const ref to file?
  }
  rw_binary(iofs, rw, m_rN, m_rcut, t_r);
  if (readQ) {
    std::cout << "(rcut = " << m_rcut << "au, rN=" << m_rN * PhysConst::aB_fm
              << " fm)\n";
  }
  // read/write potentials::
  rw_binary(iofs, rw, mVu, mVh, mVl, mHm, mVwk);

  const auto grid_same =
      (r.size() == t_r.size() && (qip::compare_eps(r, t_r).first < 1.0e-3));

  if (readQ && !grid_same) {
    std::cout << "Interpolating QED rad-pot onto current grid.\n";

    if (!mVu.empty())
      mVu = Interpolator::interpolate(t_r, mVu, r);
    if (!mVh.empty())
      mVh = Interpolator::interpolate(t_r, mVh, r);
    if (!mVl.empty())
      mVl = Interpolator::interpolate(t_r, mVl, r);
    if (!mHm.empty())
      mHm = Interpolator::interpolate(t_r, mHm, r);
    if (!mVwk.empty())
      mVwk = Interpolator::interpolate(t_r, mVwk, r);
  }

  return true;
}

//==============================================================================
void RadPot::form_potentials(const std::vector<double> &r) {
  if (print) {
    std::cout << "Forming QED radiative potential ";
    if (m_rcut > 0.0) {
      std::cout << "(rcut = " << m_rcut << ") ";
    }
    std::cout << "with Rn = " << m_rN * PhysConst::aB_fm << "fm\n";
  }

  if (m_f.u != 0.0) {
    if (print)
      std::cout << "Uehling; scale=" << m_f.u << "\n";
    mVu = fill(FGRP::V_Uehling, r);
  }
  if (m_f.h != 0.0) {
    if (print)
      std::cout << "Self-energy (high freq); scale=" << m_f.h << "\n";
    const auto stride = std::max(r.size() / 400, 2ul);
    // Stride of 5 with 2000 points
    mVh = fill(FGRP::V_SEh, r, stride);
  }
  if (m_f.l != 0.0) {
    if (print)
      std::cout << "Self-energy (low freq); scale=" << m_f.l << "\n";
    mVl = fill(FGRP::V_SEl, r);
  }
  if (m_f.m != 0.0) {
    if (print)
      std::cout << "Self-energy (magnetic); scale=" << m_f.m << "\n";
    mHm = fill(FGRP::H_Magnetic, r);
  }
  if (m_f.wk != 0.0) {
    if (print)
      std::cout << "Wickman-Kroll; scale=" << m_f.wk << "\n";
    mVwk = fill(FGRP::V_WK, r);
  }
}

//==============================================================================
std::vector<double> RadPot::Vel(int l) const {
  using namespace qip::overloads;
  const auto a = FGRP::Fit::A(m_Z, l);
  const auto b = FGRP::Fit::B(m_Z, l);
  const auto xl = m_xl(l);
  return xl * ((m_f.u * mVu) + (a * m_f.h * mVh) + (b * m_f.l * mVl) +
               (m_f.wk * mVwk));
}
//------------------------------------------------------------------------------
std::vector<double> RadPot::Hmag(int l) const {
  using namespace qip::overloads;
  const auto xl = m_xl(l);
  return xl * m_f.m * mHm;
}
//------------------------------------------------------------------------------
std::vector<double> RadPot::Vu(int l) const {
  using namespace qip::overloads;
  const auto xl = m_xl(l);
  return xl * m_f.u * mVu;
}
//------------------------------------------------------------------------------
std::vector<double> RadPot::Vl(int l) const {
  using namespace qip::overloads;
  const auto b = FGRP::Fit::B(m_Z, l);
  const auto xl = m_xl(l);
  return xl * b * m_f.l * mVl;
}
//------------------------------------------------------------------------------
std::vector<double> RadPot::Vh(int l) const {
  using namespace qip::overloads;
  const auto a = FGRP::Fit::A(m_Z, l);
  const auto xl = m_xl(l);
  return xl * a * m_f.h * mVh;
}

//==============================================================================
RadPot ConstructRadPot(const std::vector<double> &r, double Z_eff, double rN_au,
                       const IO::InputBlock &input, bool print,
                       bool do_readwrite) {

  input.check(
      {{"",
        "QED Radiative potential will be included if this block is present"},
       {"", "The following 5 are all doubles. Scale to include * potential; "
            "usually either 0.0 or 1.0, but can take any value:"},
       {"Ueh", "  Uehling (vacuum pol). [1.0]"},
       {"SE_h", "  self-energy high-freq electric. [1.0]"},
       {"SE_l", "  self-energy low-freq electric. [1.0]"},
       {"SE_m", "  self-energy magnetic. [1.0]"},
       {"WK", "  Wickman-Kroll. [0.0]"},
       {"rcut", "Maximum radius (au) to calculate Rad Pot for [5.0]"},
       {"scale_rN", "Scale factor for Nuclear size. 0 for pointlike, 1 for "
                    "typical [1.0]"},
       {"scale_l", "List of doubles. Extra scaling factor for each l e.g., "
                   "1,0,1 => include for s and d, but not for p [1.0]"}});

  const auto x_Ueh = input.get({"RadPot"}, "Ueh", 1.0);
  const auto x_SEe_h = input.get({"RadPot"}, "SE_h", 1.0);
  const auto x_SEe_l = input.get({"RadPot"}, "SE_l", 1.0);
  const auto x_SEm = input.get({"RadPot"}, "SE_m", 1.0);
  const auto x_wk = input.get({"RadPot"}, "WK", 0.0);
  const auto rcut = input.get({"RadPot"}, "rcut", 5.0);
  const auto scale_rN = input.get({"RadPot"}, "scale_rN", 1.0);
  const auto x_spd = input.get({"RadPot"}, "scale_l", std::vector{1.0});

  return RadPot(r, Z_eff, scale_rN * rN_au, rcut,
                {x_Ueh, x_SEe_h, x_SEe_l, x_SEm, x_wk}, x_spd, print,
                do_readwrite);
}

} // namespace QED
