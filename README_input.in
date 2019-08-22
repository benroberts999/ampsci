// This is an anotated input file, explaining what each available option is.
// Most options have a default, and can be left blank (or removed entirely).
// A few inputs (e.g., Atom/Z, and HartreeFock/core) are required
// Use c++ style comments (comments will not be read in)
// This file should work as an input file (it will run every available module).
// For now, no spell-checker, so mispelt options will likely just be ignored.

Atom {
  Z = Cs; // required
  A;
  varAlpha2;
}
// Atom block:
// Z: Can enter as string (Z=Cs) or integer (Z=55). Required (no default)
// A: nuclear mass number. Leave blank to look up default value
// varAlpha2: scale factor for (inverse) speed of light^2 [alpha in atomic un.]
// default=1. alpha^2 = varAlpha2 * alpha_real^2.
// put very small number to achieve non-relativistic limit (useful for tests)

HartreeFock {
  core = [Xe]; // required
  valence = 7sp5df;
  sortOutput;
  method;
  convergence;
  orthonormaliseValence;
}
// core: Core configuration. Requied (no default)
// Format: [Noble gas],extra (comma separated, no spaces)
// can enter entire sting also, e.g., 1s2,2s2,2p6,....
// (As well as Noble gas, can use Zn,Cd,Hg,Cn). Can also add negative values
// E.g. (V^N-1):
//    * For Cs: '[Xe]'
//    * For Au: '[Xe],4f14,5d10' or '[Hg],6s-2'
//    * For Tl: '[Xe],4f14,5d10,6s2' or '[Hg]'
//    * For I (V^N): '[Cd],5p5' or '[Xe],5p-1'
// enter like: []  (or 1s0) to do H-like system
// valence: which valence states to calculate
// e.g., "7sp5df" will do s and p states up to n=7, and d and f up to n=5
// sortOutput: Sort output by energy. true(dflt) or false.
// method: which method to use. can be:
// HartreeFock(default), ApproxHF, Hartree, GreenPRM, TietzPRM
// convergence: level we try to converge to. Can be blank
// orthonormaliseValence: orthogonalise valence states? false by default

Nucleus {
  rrms;
  type;
  skin_t;
}
// All of these can be left blank, in which case default values will be looked
// up (using the A value given above. A can also be given in this block)
// rrms: nuclear charge radius (in femptometres = 10^-15m)
// type: Which distribution to use for nucleus? Options are:
// Fermi (default), spherical, zero
// skin_t: skin thickness [only used by Fermi distro] 2.3 by default

Grid {
  r0 = 1e-6;
  rmax = 120.0;
  ngp = 1600;
  type;
  b;
}
// r0: grid starting point (in atomic units)
// rmax: Final grid point (in atomic units)
// ngp: number of points in the grid
// type: What type of grid to use? options are:
// loglinear (default), logarithmic, linear
// Note: 'linear' grid requires a _very_ large number of points to work,
// and should essentially never be used.
// b: only used for loglinear grid; the grid is roughly logarithmic below
// this value, and linear after it. Default is 4.0 (atomic units)

//******************************************************************************
//******************************************************************************
/*
Matrix Elements:
Each MatrixElements block will be run in order, calculate (reduced) matrix
elements of given operator.
You can comment-out just the block name, and the block will be skiped.
There are some options that apply for any operator; and then there are some
options specific to each operator
*/

MatrixElements::E1 {
  printBoth;    // applies to all operators:
  onlyDiagonal; // applies to all operators:
}
// This block details common options, that apply to all operators:
// printBoth: Print <a|h|b> and <b|h|a> ? false by default.
// For _some_ operators (e.g., involving derivatives), this is a good test of
// numerical error. For most operators though, values will be trivially the
// same (note: reduced matrix elements, sign may be different).
// onlyDiagonal: False by default. If true, will only print diagonal MEs <a|h|a>

// Following show options for each available operator:

MatrixElements::E1 {
  gauge; // lform (default), or vform
}

MatrixElements::r {
  power = 1;
  // with power = n, will calculate ME's of operator |r|^n (any real number)
}

MatrixElements::pnc {
  c;
  t;
}
// spin-independent (Qw) PNC operator
// Output given in units of i(-Q/N)e-11
// c: half-density radius, and t: skin thickness. Used for Nuclear distribution
// Will use default values (from _charge_ distribution) if none given

MatrixElements::hfs {
  F(r); // = ball, shell, pointlike, VolotkaBW, doublyOddBW
  mu;
  I;
  rrms;
  printF;
  units = MHz; // will be in atomic units by default.
  // -----
  // the following are only used in the "VolotkaBW" case
  // both are optional (will be deduced otherwise)
  parity; // parity of unpaired valence nucleon
  gl;     // =1 for valence proton, =0 for valence neutron
  // -----
  // The following are only read in if F(r) = doublyOddBW
  // but are _required_  in that case
  // only used for doubly-odd (current values are for 212-Fr)
  mu1 = 4.00;
  I1 = 4.5;
  l1 = 5.;
  gl1 = 1;
  I2 = 0.5;
  l2 = 1.;
}
// mu, I, rrms are mu (in nuc. magnetons), nuclear spin and nuclear rms radius,
// respectively. If none given, looks up default values for given Z,A (Atom)
// F(r): which model to use for nuclear magnetisation distro. default is ball,
// VolotkaBW uses single-particle formula from:
// Volotka et al, Phys. Rev. A 78, 062507 (2008).
// printF: true or false. Will print value of F(r) to a file
// Note: gives reduced matrix elements, not A constants (see below for those)
// will add option to this for A constants soon (?)

//******************************************************************************
//******************************************************************************

// Modules:
// (work in same way as matrix elements)

Module::Tests {
  // tests of numerical errors:
  orthonormal = true; // check orthonormality
  Hamiltonian = true; // check eigenvalues of Hamiltonian
}

Module::BohrWeisskopf {
  mu;
  I;
  rrms;
  // only used for doubly-odd (current values for 212-Fr)
  mu1 = 4.00;
  I1 = 4.5;
  l1 = 5.;
  gl1 = 1;
  I2 = 0.5;
  l2 = 1.;
  printF = false;
}
// Calculates Bohr-Weisskopf effect. Takes same input as MatrixElements::hfs
// (except doen't take in F(r), since it runs for each F(r))

Module::WriteOrbitals { label = outputLabel; }
// Writes the core + valence orbitals (and the total electron density) to a
// file, in GNUplot friendly format.
// The (optional) label will be appended to
// the output file name Plot in GNUPLOT. For example, try this:
// plot for [i=2:20] "file.txt" u 1:i every :::0::0  w l t columnheader(i)
// plot for [i=2:20] "file.txt" u 1:i every :::1::1  w l t columnheader(i)
// plot "file.txt" u 1:2 every :::2::2  w l t "Core Density"

Module::FitParametric {
  method = Green;     // Green, Tietz
  statesToFit = core; // core, valence, both
  fitWorst;           // false (default), true;
}
// Performs a 2D fit to determine the best-fit values for the given
// two-parameter parametric potential (Green, or Tietz potentials)
// Does fit to Hartree-Fock energies.
// (easy to update so that it takes a list of experimental energies instead, but
// haven't done that yet.)
// fitWorst: if true, will optimise fit for the worst
// state. If false, uses least squares for the fit. False is defualt

Module::pnc { transition = 6, -1, 5, 2; }
// Calculates pnc amplitude. This is just a basic example.
// (uses valence states for sum over intermediate states
// so.. need a large rmax (Grid/rmax), and many p valence states..)
// input order is: n_a, kappa_a, n_b, kappa_b for a |a> -> |b> transition

Module::AtomicKernal {
  Emin = 0.01;
  Emax = 4.0;
  Esteps = 25;
  qmin = 0.001;
  qmax = 4.0;
  qsteps = 100;
  max_l_bound = 1;
  max_L = 2;
  output_text = true;
  output_binary = true;
  label = test_new;
  use_plane_waves = false;
}
// Calculates the "Atomic Kernal" (for scattering/ionisation) for each core
// orbital, as a function of momentum transfer (q), and energy deposition (dE).
// Writes result to human-readable (and gnuplot-friendly) file, and/or binary.
// * For definitions/details, see:
// B.M. Roberts, V.V. Flambaum, arXiv:1904.07127 (2019).
// B.M.Roberts, V.A.Dzuba, V.V.Flambaum, M.Pospelov, Y.V.Stadnik, Phys.Rev.D 93,
// 115037 (2016) [arXiv:1604.04559]
// * Uses self-consistent Hartree Fock method
// (optionally, can use parametric potential, which is faster but less accurate)
// * Note: need quite a dense grid [large number of points] for
//   * a) highly oscillating J_L function at low r, and
//   * b) to solve equation for high-energy continuum states.
// * Sums over 'all' continuum angular momentum states (and multipolarities)
//   * Maximum values for l are input parameters
// Binary output from this program is read in by dmeXSection program
// Note: tested only for neutral atoms (V^N potential).
// Also: tested mainly for high values of q
