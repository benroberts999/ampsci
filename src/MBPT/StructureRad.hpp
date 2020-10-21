#pragma once
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>
class Wavefunction;
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}
namespace HF {
class HartreeFock;
}

namespace MBPT {

class StructureRad {

  /*
Note: the way this is written, SPLINE states, not HF should be used for <w,v}
(legs of diagram). HF states prob fine, so long as valence/basis is orthog. But
definitely, cannot use Brueckner orbs. The reason is that {w,v} "mix" in with
the basis states, in 'Z' integral. Note: It could be re-written to avoid this
"issue", but it's not clear if that is actually better.

For now, does not include RPA, but that is a trivial addition.

NOT TESTED

XX ALSO: add new Yk.Qk, Xk, Wk, Zk etc. tests to Coulomb Tests!

nb: With ~30 splines up to 'g', takes ~ 10 mins to calc each term.
Can probably be made much faster..

-> Can put energy cut-off in denominator
-> can have a min and max 'n'

1. Transform from X,Z -> Q,W/P
2. Try to make 'v' and 'w' appear _first_ (so always uses actual states)

  */

public:
  StructureRad(const std::vector<DiracSpinor> &basis, double en_core,
               std::pair<int, int> nminmax = {0, 999});

private:
  std::vector<DiracSpinor> mBasis;
  Coulomb::YkTable mY;

  // XXX Energy threshold; when denominators larger than this, don't calc.
  // nb: might not be good... just a test for now
  // double e_thresh = 60.0; //??

  // Fa.en < en_core ==> core state!
  // nb: take as 0.5*[maxE(core)-minE(val)]

  // nb: it seems conter-intuative, but this copy makes it FASTER!
  std::vector<DiracSpinor> mCore{}, mExcited{};

public:
  std::pair<double, double> srTB(const DiracOperator::TensorOperator *const h,
                                 const DiracSpinor &w, const DiracSpinor &v,
                                 double omega) const;

  double srC(const DiracOperator::TensorOperator *const h, const DiracSpinor &w,
             const DiracSpinor &v) const;

public:
  double t1(int K, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double t2(int K, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double t3(int K, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double t4(int K, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &c) const;

  //--

  double b1(int K, const DiracSpinor &w, const DiracSpinor &c,
            const DiracSpinor &v, const DiracSpinor &r) const;
  double b2(int K, const DiracSpinor &w, const DiracSpinor &c,
            const DiracSpinor &v, const DiracSpinor &r) const;
  double b3(int K, const DiracSpinor &w, const DiracSpinor &c,
            const DiracSpinor &v, const DiracSpinor &r) const;
  double b4(int K, const DiracSpinor &w, const DiracSpinor &c,
            const DiracSpinor &v, const DiracSpinor &r) const;

  //--

  double c1(int K, const DiracSpinor &w, const DiracSpinor &a,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double c2(int K, const DiracSpinor &w, const DiracSpinor &a,
            const DiracSpinor &v, const DiracSpinor &c) const;
  double c3(int K, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &m) const;
  double c4(int K, const DiracSpinor &w, const DiracSpinor &r,
            const DiracSpinor &v, const DiracSpinor &m) const;
};

} // namespace MBPT
