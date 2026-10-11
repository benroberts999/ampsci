#pragma once
#include <vector>

//! Set of two-parameter parametric electron potentials
namespace Parametric {

/*!
  @brief Parametric potential type.

  @details
  - Green : Green potential, parameters H, d (see green)
  - Tietz : Tietz potential, parameters g, t (see tietz)
*/
enum class Type { Green, Tietz };

/*!
  @brief Green parametric potential at r; electron part only (no nuclear).

  @details
  \f[
    V(r) = \frac{Z}{r} - \frac{1}{r}
      \left[ 1 + \frac{Z-1}{H\,(e^{r/d} - 1) + 1} \right].
  \f]
  The total potential is V(r) - Z/r. Returns 0 at r = 0.
*/
double green(int Z, double r, double H, double d);

/*!
  @brief Tietz parametric potential at r; electron part only (no nuclear).

  @details
  As implemented, in terms of the parameter names here:
  \f[
    V(r) = \frac{Z}{r} - \frac{1}{r}
      \left[ 1 + \frac{(Z-1)\,e^{-t r}}{(1 + g r)^2} \right].
  \f]
  The total potential is V(r) - Z/r. Returns 0 at r = 0.

  @warning The definition names the 3rd/4th parameters (t, g), opposite to
           this declaration (g, t). As called, g and t therefore take the
           opposite roles to the standard Tietz form,
           \f$ e^{-g r}/(1 + t r)^2 \f$.
*/
double tietz(int Z, double r, double g, double t);

/*!
  @brief Default Green parameters H and d for given Z, tuned to match HF core.

  @details
  Tabulated for selected Z, linearly interpolated in between. For Z > 113,
  H = 5.0, d = 0.75.

  @return Always 0; H and d are set by reference.

  @warning For Z = 1, H = d = 0, for which green() is nan.
*/
int defaultGreenCore(int z, double &H, double &d);

/*!
  @brief Default Green parameters H and d for given Z, tuned to experiment.

  @details
  Fitted to the lowest valence J=1/2 states for selected Z; crude quadratic
  fit in Z otherwise.

  @return Always 0; H and d are set by reference.
*/
int defaultGreen(int z, double &H, double &d);

/*!
  @brief Default Tietz parameters g and t for given Z, tuned to experiment.

  @details
  Fitted to the lowest valence J=1/2 states for selected Z; crude quadratic
  fit in Z otherwise.

  @return Always 0; g and t are set by reference.
*/
int defaultTietz(int z, double &g, double &t);

/*!
  @brief Green parametric potential (electron part) on given radial points.

  @details
  See green(). If H*d is (near) zero, the core-tuned defaults are used
  (defaultGreenCore).

  @param z        Nuclear charge
  @param r_array  Radial points
  @param H        Green parameter H (0: default)
  @param d        Green parameter d (0: default)

  @return V(r) at each point of @p r_array; excludes nuclear potential.
*/
std::vector<double> GreenPotential(int z, const std::vector<double> &r_array,
                                   double H = 0, double d = 0);

/*!
  @brief Tietz parametric potential (electron part) on given radial points.

  @details
  See tietz(). If g*t is (near) zero, the experiment-tuned defaults are used
  (defaultTietz).

  @param z        Nuclear charge
  @param r_array  Radial points
  @param g        Tietz parameter g (0: default)
  @param t        Tietz parameter t (0: default)

  @return V(r) at each point of @p r_array; excludes nuclear potential.
*/
std::vector<double> TietzPotential(int z, const std::vector<double> &r_array,
                                   double g = 0, double t = 0);
} // namespace Parametric
