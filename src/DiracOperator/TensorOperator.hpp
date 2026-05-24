#pragma once
#include "Angular/Wigner369j.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include <cmath>
#include <memory>
#include <string>
#include <vector>

/*!
  @brief Dirac operators: TensorOperator base class and derived implementations for single-particle (one-body) spherical tensor operators.

  @details
  This namespace contains the TensorOperator base class and all derived
  operator classes. Operators act on DiracSpinor objects and compute reduced
  matrix elements (RMEs) of the form:

  \f[
    \langle a \| \hat{h} \| b \rangle = A_{ab} \cdot R_{ab}
  \f]

  where \f$ A_{ab} \f$ is a purely angular factor and \f$ R_{ab} \f$ is a radial
  integral. New operators are implemented by deriving from TensorOperator (or
  ScalarOperator for rank-0 operators) and overriding the relevant virtual
  functions; see TensorOperator for details.

  Operators are typically constructed via @ref generate(), which takes a name
  string and an InputBlockLegacy; see GenerateOperator.hpp for the full list of
  available operators.

  See @ref TensorOperator class documentation for main descriptions.
  All usable tensor operators dervive from that virtual class.

  From the command line, use
  ```shell
    ampsci -o
  ```
  to get a list of available operators.
  For a specific operator 'OperatorName', use:
  ```shell
    ampsci -o OperatorName
  ```
  to see available run-time options for that operator.
  These options are passed to @ref generate().

  @note New operator classes must be registered in the @ref operator_list in
  GenerateOperator.hpp to be accessible at runtime.
*/
namespace DiracOperator {

//! Parity of operator
enum class Parity { even, odd, Error };

/*!
  @brief Specifies whether an operator's matrix elements are real or imaginary.
  @details
  Operators must have purely real or purely imaginary reduced matrix elements.
  For imaginary operators, the imaginary part is computed and stored; the real
  part is identically zero. This distinction affects the Hermitian symmetry
  relation 
  \f[ 
    \langle b \| \hat{h} \| a \rangle = \pm \langle a \| \hat{h} \| b \rangle 
  \f]
*/
enum class Realness { real, imaginary, Error };

//! Type of matrix element returned
/*! @details
 Specifies which form of matrix element is used or returned.

  @value Reduced    : Reduced matrix element
  @value Stretched  : Matrix element for stretched states (m = j)
  @value HFConstant : Hyperfine (A,B,C...) constant

  For off-diagonal, j:=min(ja,jb)
*/
enum class MatrixElementType { Reduced, Stretched, HFConstant, Error };

//! Convert string to Parity
Parity parse_Parity(const std::string &s);

//! Convert Parity to string
std::string parse_Parity(Parity p);

//! Convert string to Realness
Realness parse_Realness(const std::string &s);

//! Convert Realness to string
std::string parse_Realness(Realness r);

//! Convert string to MatrixElementType
MatrixElementType parse_MatrixElementType(const std::string &s);

//! Convert MatrixElementType to string
std::string parse_MatrixElementType(MatrixElementType t);

//==============================================================================
/*!
  @brief General tensor operator (virtual base class); 
  all single-particle (one-body) tenosor operators derive from this.

  @details

  Represents a single-particle (one-body) tensor operator 
  \f$ \hat{h} \f$ of rank-\f$ k \f$ and parity \f$ \pi \f$, 
  thats acts on Dirac spinors. 
  The core quantity is the _reduced matrix element_ (RME):

  \f[
    \langle a \| \hat{h} \| b \rangle \equiv A_{ab} \cdot R_{ab}
  \f]

  where \f$ A_{ab} \f$ is the angular factor (see angularF()) and \f$ R_{ab} \f$ 
  is the radial integral (see radialIntegral()).

  The default radial integral is:

  \f[
    R_{ab} = c \int_0^\infty v(r)\left(
      C_{ff}\,f_a f_b + C_{fg}\,f_a g_b +
      C_{gf}\,g_a f_b + C_{gg}\,g_a g_b
    \right)\,{\rm d}r
  \f]

  where \f$ c \f$ is an overall multiplicative constant (1 by default),
  \f$ v(r) \f$ is an optional radial function (stored in m_vec), 
  and the \f$ C_{xy} \f$ are angular coefficients returned by @ref angularCff() etc.

  - Operators must be spherical tensor operators of defined rank and parity
  - Operators must have pure real, or pure imaginary matrix elements
  - If matrix elements are imaginary, class returns the imaginary part
  - For operators with imaginary matrix elements, the @ref Realness member variable (m_Realness) must be set to Realness::imaginary via the virtual base class constructor @ref TensorOperator() (impacts symmetry)
  - (For brevity, we may refer to operators whose matrix elements are imaginary as 'imaginary operators'. This is sometimes confusing.)

  ## Implementing a derived operator

  TensorOperator is a virtual base class and cannot be constructed directly.

  @note All derived operators must implement angularF()
  
  Beyond that, there are two levels of customisation:

  ### Standard case -- override angular factors only:

  Implement angularF() (required), and optionally 
  @ref angularCff(), angularCgg(), angularCfg(), angularCgf().
  The radial integral is then handled
  automatically by the default @ref radial_rhs() and @ref radialIntegral().
  The \f$ v(r) \f$ radial function, and overall constant \f$ c \f$ are set
  within the constructor @ref TensorOperator::TensorOperator(), which should
  be called by any deriving class.
  Then, the default radial integral as above will be used.

  ### Non-standard case -- override the radial functions:

  If the radial integral cannot be expressed in the above standard/default form 
  (e.g., it involves derivatives, non-local terms, or more complicated radial structure), 
  you shuold override radial_rhs() **and** radialIntegral() directly. 
  In that case both should be consistent with each other:

  - radial_rhs(ka, Fb) returns the spinor \f$ \delta F_b \f$ such that

  \f[ F_a \cdot \delta F_b = R_{ab} \f] 

  - radialIntegral(Fa, Fb) returns \f$ R_{ab} \f$ directly
  - Generally, this should be directly checked with a unit test

  @note 
  - If radialIntegral() is overridden then radial_rhs() _must_ be overridden
  (and vice versa), or results will be inconsistent. 
  It's OK (but inefficient) to define radialIntegral(a,b) = a*radial_rhs(b).
  But that is not the default.

  @warning 
  - Frequency-dependent operators: if the operator depends on the
  transition frequency, or energy/momentum transfer,
  (e.g., dynamic multipole operators), pass
  freq_dep=true to the constructor and override updateFrequency(). This
  function must be called with the current frequency before computing
  any matrix elements. Failing to override it will abort at runtime.

  @note
  - You may not construct a TensorOperator directly. Construct one of the
  derived classes; see [Operators/include.hpp](Operators/include.hpp) for the list.

  --

  @note
  - To expose a new operator at runtime (e.g., via `ampsci -o`), add it to
  the @ref operator_list in GenerateOperator.hpp and provide a corresponding
  `generate_XYZ()` factory function.
*/
class TensorOperator {
protected:
  /*! 
    @brief Constructs a specific tensor operator. Called by derived classes.
    @details
    Initialises the basic properties of a tensor operator.

    @param rank_k     Tensor rank k of the operator.
    @param pi         Parity of the operator (Parity::even or Parity::odd).
    @param constant   Multiplicative constant c, included in the radial integral.
    @param vec        Overall v=f(r) radial function, as std::vector. 
                      May be impty, in which case it is not used (equivilant to vector of 1)
    @param RorI       Specifies whether the matrix element is real or imaginary
                      (Realness::real or Realness::imaginary).
    @param freq_dep   Indicates whether the operator depends on frequency
                      (e.g., dynamic external fields).

    @note TensorOperator is a virtual base class and may not be constructed
    directly. It is intended to be initialised by derived operator classes.
  */
  TensorOperator(int rank_k, Parity pi, double constant = 1.0,
                 const std::vector<double> &vec = {},
                 Realness RorI = Realness::real, bool freq_dep = false)
    : m_rank(rank_k),
      m_parity(pi),
      m_Realness(RorI),
      m_freqDependantQ(freq_dep),
      m_constant(constant),
      m_vec(vec) {};

public:
  //! Used to pass generic parameters to update() function
  /*! @details
   Override in derived class, then use as:

   auto const* params = dynamic_cast<const Params*>(&p);

   Rarely used.
  */
  struct Params {
    virtual ~Params() = default;
  };

public:
  virtual ~TensorOperator() = default;
  TensorOperator(const TensorOperator &) = default;
  TensorOperator &operator=(const TensorOperator &) = default;
  TensorOperator(TensorOperator &&) = default;
  TensorOperator &operator=(TensorOperator &&) = default;

protected:
  int m_rank;
  Parity m_parity;
  Realness m_Realness;
  bool m_freqDependantQ{false};

protected:
  // these may be updated for frequency-dependant operators
  double m_constant; // included in radial integral
  std::vector<double> m_vec;

public:
  //! Returns true if the operator is frequency-dependent (requires updateFrequency() calls).
  bool freqDependantQ() const { return m_freqDependantQ; }

public:
  //! Returns true if <a|h|b> = 0 by rank/parity selection rules.
  bool isZero(int ka, int kb) const;
  //! Overload taking DiracSpinors; forwards to isZero(ka, kb).
  bool isZero(const DiracSpinor &Fa, const DiracSpinor &Fb) const;

  //! Returns true if the matrix element is non-zero by angular momentum and parity selection rules (arguments are 2j and pi as integers).
  bool selectrion_rule(int twoJA, int piA, int twoJB, int piB) const {
    if (twoJA == twoJB && twoJA == 0.0)
      return false;

    if (Angular::triangle(twoJA, twoJB, 2 * m_rank) == 0)
      return false;

    return (m_parity == Parity::even) == (piA == piB);
  }

  //! @brief Updates the operator for a new frequency omega.
  /*!
    @details
    Must be overridden by any frequency-dependent operator (i.e., where
    freqDependantQ() returns true). Called before computing matrix elements
    whenever the frequency changes.

    The base class implementation aborts -- if a frequency-dependent operator
    is constructed but this function is not overridden, it will abort at
    runtime when called.

    @param omega   Frequency in atomic units.

    @warning Must be implemented in any derived class that sets freq_dep=true.
    Calling this on a non-frequency-dependent operator is a logic error.
  */
  virtual void updateFrequency(const double) {
    std::cout << "Must reimplement updateFrequency()\n";
    std::cout << this->name() << "\n";
    std::abort();
  };

  /*!
    Updates the rank of operator (rarely used). Generally also updates parity

    @details 
    Use with caution.
    @warning Will usually have to call updateFrequency() after this! 
    Updating rank often breaks frequency-dependence of radial functions 
    (e.g., spherical Bessel functions)
  */
  virtual void updateRank(int) {
    std::cout << "Must reimplement updateRank() is needed\n";
    std::cout << this->name() << "\n";
    std::abort();
  }

  //! Returns a const ref to the stored vector v
  const std::vector<double> &getv() const { return m_vec; }

  //! Returns the "overall" constant c
  double getc() const { return m_constant; }

  //! returns true if operator is imaginary (has imag MEs)
  bool imaginaryQ() const { return (m_Realness == Realness::imaginary); }

  //! Rank k of operator
  int rank() const { return m_rank; }

  //! returns parity, as integer (+1 or -1)
  int parity() const { return (m_parity == Parity::even) ? 1 : -1; }

  //! returns relative sign between <a||x||b> and <b||x||a>
  int symm_sign(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
    const auto sra_i = imaginaryQ() ? -1 : 1;
    const auto sra = Angular::neg1pow_2(Fa.twoj() - Fb.twoj());
    return sra_i * sra;
  }

  //! Returns "name" of operator (e.g., 'E1')
  virtual std::string name() const { return "Operator"; };
  //! Returns units of operator as a string (usually au, may be MHz, etc.)
  virtual std::string units() const { return "au"; };

public:
  // These are needed for radial integrals
  // Usually just constants, but can also be functions of kappa

  /*!
    @brief Angular coefficient C_ff for the f_a*f_b term of the radial integral.
    
    @details
    The default radial integral is structured as:

    \f[
      R_{ab} = c\int_0^\infty v(r)\left(
        C_{ff}\,f_a f_b + C_{fg}\,f_a g_b +
        C_{gf}\,g_a f_b + C_{gg}\,g_a g_b
      \right)\,{\rm d}r
    \f]

    These coefficients are often constants, but may depend on
    \f$ \kappa_a, \kappa_b \f$ for operators with angular-momentum-dependent
    coupling between large and small components (e.g., spin-dependent
    operators). Override in derived classes as needed.

    @param kappa_a  kappa \f$ \kappa_a \f$ for left-hand-side (bra)
    @param kappa_b  kappa \f$ \kappa_b \f$ for right-hand-side (ket)

    @note Only relevant when using the default radial_rhs()/radialIntegral().
    If those are overridden, these are not called.
  */
  virtual double angularCff(int kappa_a, int kappa_b) const {
    (void)kappa_a, (void)kappa_b;
    return 1.0;
  }
  //! Angular coefficient C_gg for the g_a*g_b term of the radial integral.
  virtual double angularCgg(int, int) const { return 1.0; }
  //! Angular coefficient C_fg for the f_a*g_b term of the radial integral.
  virtual double angularCfg(int, int) const { return 0.0; }
  //! Angular coefficient C_gf for the g_a*f_b term of the radial integral.
  virtual double angularCgf(int, int) const { return 0.0; }

  /*!
    @brief Dispatches to angularCff/fg/gf/gg based on component indices x, y.
    @details
    Convenience dispatcher: x and y index the spinor component (0 = large/f,
    1 = small/g) of the bra and ket respectively, and returns the
    corresponding angular coefficient C_xy(ka, kb).

    See @ref angularCff()

    | x | y | returns        |
    |---|---|----------------|
    | 0 | 0 | angularCff()   |
    | 0 | 1 | angularCfg()   |
    | 1 | 0 | angularCgf()   |
    | 1 | 1 | angularCgg()   |

    @param x       Bra component index: 0=f (large), 1=g (small).
    @param y       Ket component index: 0=f (large), 1=g (small).
    @param kappa_a  kappa \f$ \kappa_a \f$ for left-hand-side (bra)
    @param kappa_b  kappa \f$ \kappa_b \f$ for right-hand-side (ket)

    @note Only relevant when using the default radial_rhs()/radialIntegral().
    If those are overridden, these are not called.
  */
  double angularCxy(uint8_t x, uint8_t y, int kappa_a, int kappa_b) const {
    assert((x <= 1 && y <= 1) && "x and y must be 0 or 1 in angularCxy");
    return x == 0 ? (y == 0 ? angularCff(kappa_a, kappa_b) :
                              angularCfg(kappa_a, kappa_b)) :
                    (y == 0 ? angularCgf(kappa_a, kappa_b) :
                              angularCgg(kappa_a, kappa_b));
  }

public:
  /*!
    @brief Angular factor A_ab linking the radial integral to the RME.
    
    @details
    All derived operators must implement this. It gives the purely angular
    part of the reduced matrix element:

    \f[
      \langle a \| \hat{h} \| b \rangle \equiv A_{ab} \cdot R_{ab}
    \f]

    where \f$ R_{ab} \f$ is returned by radialIntegral(). For most operators,
    \f$ A_{ab} \f$ is a product of Clebsch-Gordan / 3j coefficients and
    depends only on \f$ \kappa_a, \kappa_b \f$ 
    (and the rank \f$ k \f$ and parity \f$ \pi \f$ of the operator).

    @note This is a pure virtual function -- every derived operator must
    provide an implementation.
  */
  virtual double angularF(const int, const int) const = 0;

  //! Returns a polymorphic copy of the operator at its current state,
  //! or nullptr if cloning is not supported by the derived class.
  virtual std::unique_ptr<TensorOperator> clone() const { return nullptr; }

  //! @brief Computes the right-hand spinor dF_b for the radial integral.
  /*!
    @details
    Returns \f$ \delta F_b \f$ such that the radial integral satisfies:

    \f[
      R_{ab} = F_a \cdot \delta F_b
             = \int_0^\infty \left(f_a\,\delta f_b + g_a\,\delta g_b\right)\,{\rm d}r
    \f]

    The default implementation constructs \f$ \delta F_b \f$ using the stored
    radial function \f$ v(r) \f$ and the angular coefficients:

    \f[
      \delta F_b(r) = c\,v(r)
      \begin{pmatrix}
        C_{ff}\,f_b(r) + C_{fg}\,g_b(r) \\
        C_{gf}\,f_b(r) + C_{gg}\,g_b(r)
      \end{pmatrix}
    \f]

    This is used by reduced_rhs() to build \f$ \langle a \| \hat{h} \| b \rangle \f$ 
    as a spinor-valued quantity, enabling perturbation theory and TDHF.
    Override this for operators whose radial structure cannot be expressed in
    this standard form.

    @param kappa_a   Relativistic quantum number \f$ \kappa_a \f$ of the bra state
                     (needed to evaluate the angular coefficients).
    @param Fb        Ket DiracSpinor \f$ F_b \f$ .

    @return DiracSpinor \f$ \delta F_b \f$ .

    @warning If this is overridden, radialIntegral() should also be overridden
    consistently (and vice versa), so that reducedME() and reduced_rhs() remain
    consistent.
  */
  virtual DiracSpinor radial_rhs(const int kappa_a,
                                 const DiracSpinor &Fb) const;

  //! @brief Radial integral R_ab, defined by RME = angularF(a,b) * radialIntegral(a,b).
  /*!
    @details
    Returns the radial part \f$ R_{ab} \f$ of the reduced matrix element:

    \f[
      \langle a \| \hat{h} \| b \rangle = A_{ab} \cdot R_{ab}
    \f]

    where \f$ A_{ab} \f$ is angularF().

    The default implementation evaluates \f$ R_{ab} = F_a \cdot \delta F_b \f$ ,
    using the default radial structure:

    \f[
      R_{ab} = c\int_0^\infty v(r)\left(
        C_{ff}\,f_a f_b + C_{fg}\,f_a g_b +
        C_{gf}\,g_a f_b + C_{gg}\,g_a g_b
      \right)\,{\rm d}r
    \f]

    Override this for operators that do not fit this standard form. If
    radial_rhs() is also overridden, both must remain mutually consistent.

    @warning If radial_rhs() is overridden but radialIntegral() is not (or
    vice versa), reducedME() and reduced_rhs() will give inconsistent results.
  */
  virtual double radialIntegral(const DiracSpinor &Fa,
                                const DiracSpinor &Fb) const;

  /*!
    @brief 3j-symbol factor linking the full ME to the RME.
    @details
    \f[ (-1)^{j_a - m_a} \begin{pmatrix} j_a & k & j_b \\ -m_a & q & m_b \end{pmatrix} \f]
    such that \f$\langle a | \hat{h} | b \rangle = {\rm rme3js} \times \langle a \| \hat{h} \| b \rangle\f$.
    All arguments are twice the actual value (2*j, 2*m, 2*q).
  */
  double rme3js(int twoja, int twojb, int two_mb = 1, int two_q = 0) const;

  //! Overload of rme3js taking DiracSpinors.
  double rme3js(const DiracSpinor &Fa, const DiracSpinor &Fb, int two_mb = 1,
                int two_q = 0) const {
    return rme3js(Fa.twoj(), Fb.twoj(), two_mb, two_q);
  }

  //! Returns angularF(ka,kb) * radial_rhs(ka,Fb); spinor-valued RME action on Fb, used in perturbation theory/TDHF.
  DiracSpinor reduced_rhs(const int ka, const DiracSpinor &Fb) const;

  //! As reduced_rhs but for the conjugate direction; Fb * reduced_lhs(ka, Fb) = <b||h||a>.
  DiracSpinor reduced_lhs(const int ka, const DiracSpinor &Fb) const;

  //! Returns the reduced matrix element <a||h||b> = A_ab * R_ab.
  double reducedME(const DiracSpinor &Fa, const DiracSpinor &Fb) const;

  //! Returns "full" matrix element, for optional (ma, mb, q) [taken as int 2*].
  //! If not specified, returns z-component (q=0), with ma=mb=min(ja,jb)
  double fullME(const DiracSpinor &Fa, const DiracSpinor &Fb,
                std::optional<int> two_ma = std::nullopt,
                std::optional<int> two_mb = std::nullopt,
                std::optional<int> two_q = std::nullopt) const;

  //! Returns the factor to convert a reduced ME to a different form (Reduced, Stretched, or HFConstant); see MatrixElementType.
  double matel_factor(MatrixElementType type, int twoJa, int twoJb) const;

  //! Overload of matel_factor taking DiracSpinors.
  double matel_factor(MatrixElementType type, const DiracSpinor &Fa,
                      const DiracSpinor &Fb) const {
    return matel_factor(type, Fa.twoj(), Fb.twoj());
  }
};

//============================================================================
//============================================================================
/*!
  @brief Rank-0 (scalar) tensor operator; derives from TensorOperator with k=0.
  @details
  Convenience base for scalar operators. angularF() returns \f$ \sqrt{2|\kappa|} \f$ 
  when \f$ |\kappa_a|=|\kappa_b| \f$ , and zero otherwise. The four C_xy coefficients
  default to (1,0,0,1) (diagonal ff+gg), but can be set via in_g for operators
  with large/small component mixing (e.g., PNC operators).
*/
class ScalarOperator : public TensorOperator {
public:
  //! General scalar operator constructor. in_g = {C_ff, C_fg, C_gf, C_gg}.
  ScalarOperator(Parity pi, double in_coef,
                 const std::vector<double> &in_v = {},
                 const std::array<int, 4> &in_g = {1, 0, 0, 1},
                 Realness rori = Realness::real)
    : TensorOperator(0, pi, in_coef, in_v, rori),
      c_ff(in_g[0]),
      c_fg(in_g[1]),
      c_gf(in_g[2]),
      c_gg(in_g[3]) {}

  //! Convenience constructor: even-parity, real, diagonal (C_ff=C_gg=1, C_fg=C_gf=0) scalar operator.
  ScalarOperator(const std::vector<double> &in_v = {}, double in_coef = 1.0)
    : TensorOperator(0, Parity::even, in_coef, in_v),
      c_ff(1.0),
      c_fg(0.0),
      c_gf(0.0),
      c_gg(1.0) {}

public:
  virtual double angularF(const int ka, const int kb) const override {
    // For scalar operators, <a||h||b> = RadInt / 3js
    // 3js:= 1/(Sqrt[2j+1]) ... depends on m???
    // |k| = j+1/2
    return (std::abs(ka) == std::abs(kb)) ? std::sqrt(2.0 * std::abs(ka)) : 0.0;
  }

private:
  const double c_ff, c_fg, c_gf, c_gg;

protected:
  double virtual angularCff(int, int) const override { return c_ff; }
  double virtual angularCgg(int, int) const override { return c_gg; }
  double virtual angularCfg(int, int) const override { return c_fg; }
  double virtual angularCgf(int, int) const override { return c_gf; }
};

//------------------------------------------------------------------------------
//! Speacial operator: 0
class NullOperator final : public ScalarOperator {
public:
  NullOperator() : ScalarOperator(Parity::even, 0, {}) {}

protected:
  double angularCff(int, int) const override final { return 0.0; }
  double angularCgg(int, int) const override final { return 0.0; }
  double angularCfg(int, int) const override final { return 0.0; }
  double angularCgf(int, int) const override final { return 0.0; }
};

//******************************************************************************
// Helper functions: Useful for several operators
//******************************************************************************

//! Pab function: Int[ (fa*gb + pm*ga*fb) * t(r) , dr]. pm = +/-1 (usually)
/*! @details
Note: does not know selection rules, so only evaluate if non-zero
*/
double Pab(double pm, const std::vector<double> &t, const DiracSpinor &Fa,
           const DiracSpinor &Fb);

//! Rab function: Int[ (fa*fb + pm*ga*gb) * t(r) , dr]. pm = +/-1 (usually)
/*! @details
Note: does not know selection rules, so only evaluate if non-zero
*/
double Rab(double pm, const std::vector<double> &t, const DiracSpinor &Fa,
           const DiracSpinor &Fb);

//! Pab_rhs function: dF_ab += A * t(r) * (g, pm*f) , pm=+/-1 (usually).
//! NOTE: uses +=, so can combine. Ensure empty to begin.
/*! @details
Note: does not know selection rules, so only evaluate if non-zero.
Should have Fa*Pab_rhs = A * Pab.
*/
void Pab_rhs(double pm, const std::vector<double> &t, DiracSpinor *dF,
             const DiracSpinor &Fb, double A = 1.0);

//! Rab_rhs function: dF_ab += A * t(r) * (f, pm*g) , pm=+/-1 (usually).
//! NOTE: uses +=, so can combine. Ensure empty to begin.
/*! @details
Note: does not know selection rules, so only evaluate if non-zero.
Should have Fa*Rab_rhs = A * Rab.
*/
void Rab_rhs(double pm, const std::vector<double> &t, DiracSpinor *dF,
             const DiracSpinor &Fb, double A = 1.0);

//! Vab function: Int[ (fa*gb ) * t(r) , dr].
double Vab(const std::vector<double> &t, const DiracSpinor &Fa,
           const DiracSpinor &Fb);
//! Wab function: Int[ (ga*fb ) * t(r) , dr].
double Wab(const std::vector<double> &t, const DiracSpinor &Fa,
           const DiracSpinor &Fb);

//! Gab function: Int[ (ga*gb ) * t(r) , dr].
double Gab(const std::vector<double> &t, const DiracSpinor &Fa,
           const DiracSpinor &Fb);

//! Gab_rhs function: dF += a * t(r) * (0, g_b). NOTE: uses +=, so can combine.
void Gab_rhs(const std::vector<double> &t, DiracSpinor *dF,
             const DiracSpinor &Fb, double a);

// Same - for constant t(r)=c

//! Pab[1] function: Int[ (fa*gb + pm*ga*fb) , dr]. pm = +/-1 (usually)
double Pab(double pm, const DiracSpinor &Fa, const DiracSpinor &Fb);

//! Rab[1] function: Int[ (fa*fb + pm*ga*gb) , dr] = Int[ (pm-1)*ga*gb) , dr].
//! NOTE: assumes NOT diagonal, using orthogonality condition.
double Rab(double pm, const DiracSpinor &Fa, const DiracSpinor &Fb);

//! Pab_rhs[1] function: dF_ab += A  * (g, pm*f) , pm=+/-1 (usually).
void Pab_rhs(double pm, DiracSpinor *dF, const DiracSpinor &Fb, double A = 1.0);

//! Rab_rhs[1] function: dF_ab += A * (f, pm*g)  = dF_ab += A * (0, (pm-1)*g).
//! NOTE: assumes NOT diagonal, using orthogonality condition.
void Rab_rhs(double pm, DiracSpinor *dF, const DiracSpinor &Fb, double A = 1.0);

//! Gab = Int[ ga*gb , dr] - (just relativistic correction part of integral)
double Gab(const DiracSpinor &Fa, const DiracSpinor &Fb);

//! Gab_rhs(r) += a*g_b(r). Note: uses += so may be sumulative
void Gab_rhs(DiracSpinor *dF, const DiracSpinor &Fb, double a);

} // namespace DiracOperator
