\page tutorial_exotic Exotic atoms

\brief Calculations for exotic atoms (e.g., muonic)


This assumes you already have ampsci compiled and have a basic understanding of how to run and use it.

* See [Compilation](\ref compilation) for compilation instructions
* See [Basic Tutorial](\ref tutorial_basic) for getting started if unfamiliar

## Starting point:

At a minimum, you must set the `Atom`, `Grid` and `Exotic` blocks:

For example:

```java
Atom {
  Z = Cs;
  A = 133;
}

Grid {
  r0 = 1e-8;
  rmax = 1.0;
  num_points = 4000;
}

Exotic {
  states = 3sp;
}
```

The optional `Nucleus` block allows you to set nuclear charge distribution parameters. If not set, defaults (based on Z,A) will be used. Uses Fermi distribution with charge radii from Angeli tables by default.

Note:

- Use `ampsci -a` to see a list of available input blocks
- And, for example, `ampsci -a Exotic` to see list of options for 'Exotic' block

With the above options in text file `input_file.in`, run as:
`./ampsci input_file.in`

Output should be:

```text
Running for Cs, Z=55 A=133
Fermi nucleus;  r_rms = 4.8041, c_hdr = 5.67073, t = 2.3
Loglinear (b=0.333333) grid: 1e-08 -> 1.0, N=4000, du=0.00179
========================================================

Cs-133
Core: [] V^N-55 - H-like
E_c = 0.000000

---------------------------------------------
Exotic Cs
M = 206.76828270 m_e = 105.658375496267 MeV

Energies - without screening:
nk    Rinf  eps    R_rms (a0)    E (au)            E (keV)
1s+   0.05  1e-16  2.08947e-04  -2.205529115e+05  -6.001550463e+03
2s+   0.09  9e-15  6.65035e-04  -6.648031490e+04  -1.809021527e+03
3s+   0.13  4e-15  1.39951e-03  -3.131039258e+04  -8.519991861e+02
2p-   0.08  1e-15  4.65662e-04  -8.038004543e+04  -2.187252463e+03
3p-   0.13  3e-15  1.15009e-03  -3.554126473e+04  -9.671270822e+02
2p+   0.08  2e-15  4.81379e-04  -7.813279599e+04  -2.126101690e+03
3p+   0.13  2e-16  1.17161e-03  -3.492090173e+04  -9.502461450e+02
```

- nk is state (+/- means j+/-1/2),  p- = p_1/2, p+=p_3/2 etc.
- Rinf radius where wavefunction effectively goes to zero. Should be smaller than rmax from the Grid.
- eps - convergence parameter for solcing eigenvalue problem. Check carefully - should be very small number <~10^-11
- R_rms - the rms radius for this state. Helpful for determining reasonable grid parameters
- E - binding energy, in atomic units (Hartree) and keV

Hint: use `ampsci -c` from command line to get useful unit conversions.
Even better: try these out for actual conversions:
`ampsci -c 5 keV'
`ampsci -c 200 nm'

## Important: Making sure the Grid settings are OK

Very important to choose a Grid that actually works.
In the end, you should always check that making small changes to the grid does not impact the result.

- The first way to make sure Grid is OK is to check 'Rinf' (that your wavefunctions all fit inside the grid radius).
- More importantly, check 'eps':

For example, with this Grid:

```java
Grid {
  r0 = 1e-7;
  rmax = 10.0;
  num_points = 3000;
}
```

```text
Energies - without screening:
nk    Rinf  eps    R_rms (a0)    E (au)            E (keV)
1s+   0.05  7e-15  2.08947e-04  -2.205529117e+05  -6.001550468e+03
2s+   0.09  4e-15  6.65035e-04  -6.648031492e+04  -1.809021527e+03
3s+   0.13  1e-04  1.50559e-03  -3.131393670e+04  -8.520956264e+02
2p-   0.08  4e-16  4.65662e-04  -8.038004540e+04  -2.187252462e+03
3p-   0.12  1e-15  1.15009e-03  -3.554126472e+04  -9.671270819e+02
2p+   0.08  1e-15  4.81379e-04  -7.813279599e+04  -2.126101690e+03
3p+   0.13  0e+00  1.17161e-03  -3.492090173e+04  -9.502461450e+02
```

Everything is OK except for '3s+', which only converged to '1e-04'.

Using a more reasonable Grid (you can make educated guess, but best bet is usually trial+error):

```java
Grid {
  r0 = 1e-8;
  rmax = 1.0;
  num_points = 4000;
}
```

```text
Energies - without screening:
nk    Rinf  eps    R_rms (a0)    E (au)            E (keV)
1s+   0.05  1e-16  2.08947e-04  -2.205529115e+05  -6.001550463e+03
2s+   0.09  9e-15  6.65035e-04  -6.648031490e+04  -1.809021527e+03
3s+   0.13  4e-15  1.39951e-03  -3.131039258e+04  -8.519991861e+02
2p-   0.08  1e-15  4.65662e-04  -8.038004543e+04  -2.187252463e+03
3p-   0.13  3e-15  1.15009e-03  -3.554126473e+04  -9.671270822e+02
2p+   0.08  2e-15  4.81379e-04  -7.813279599e+04  -2.126101690e+03
3p+   0.13  2e-16  1.17161e-03  -3.492090173e+04  -9.502461450e+02
```

## Include effect of electron screening

Can add a Hartree-Fock block to include electrons alongside the muon

```java
HartreeFock {
  core = [Xe];
}
```

This will solve Hartree-Fock for the electrons in core configuration 'Xe', for example. Must have fewer electrons than Z to converge (since there will also be a muon). In reality, atom is usually very ionised, so probably `[He]` or someting more reasonable.
Both effects are accounted for: impact of muon on electrons, and impact of electrons on the muon.
Due to radius of their respective orbits, the effect of screening on the muon is very small.
(Effect of muon on electrons is not small - it essentially subtracts 1 from Z!)

Note: The same 'Grid' is used for both the electron and muon parts of the problem.
Therefore, you will need quite a large grid for this to work well. For example:

```java
Grid {
  r0 = 1e-8;
  rmax = 50.0; // accomodate electron _and_ muons
  num_points = 20000;
}
```

```text
Energies - without screening:
nk    Rinf  eps    R_rms (a0)    E (au)            E (keV)
1s+   0.05  2e-15  2.08947e-04  -2.205529115e+05  -6.001550463e+03
2s+   0.09  2e-16  6.65035e-04  -6.648031490e+04  -1.809021527e+03
3s+   0.13  5e-16  1.39951e-03  -3.131039258e+04  -8.519991861e+02
2p-   0.08  4e-16  4.65662e-04  -8.038004543e+04  -2.187252463e+03
3p-   0.13  6e-16  1.15009e-03  -3.554126473e+04  -9.671270822e+02
2p+   0.08  1e-15  4.81379e-04  -7.813279599e+04  -2.126101690e+03
3p+   0.13  2e-15  1.17161e-03  -3.492090173e+04  -9.502461450e+02

....

Exotic energies - with screening:
nk    Rinf  eps    R_rms (a0)    E (au)            E (keV)
1s+   0.05  7e-16  2.08947e-04  -2.202175480e+05  -5.992424758e+03
2s+   0.09  1e-15  6.65037e-04  -6.614512626e+04  -1.799900579e+03
3s+   0.13  8e-16  1.39953e-03  -3.097576144e+04  -8.428934089e+02
2p-   0.08  2e-15  4.65663e-04  -8.004476127e+04  -2.178128916e+03
3p-   0.13  9e-15  1.15011e-03  -3.520640973e+04  -9.580152136e+02
2p+   0.08  2e-16  4.81380e-04  -7.779751836e+04  -2.116978321e+03
3p+   0.13  4e-16  1.17163e-03  -3.458606467e+04  -9.411347644e+02
```

So, small, but noticable impact.

## Matrix elements, hyperfine structure, anomaly

- Use `ampsci -o` to get the list of available operators
- And, for example, `ampsci -o hfs` to see all options for the hfs operator

Can simple calculate matrix elements of any operator:

```java
Module::matrixElements{
  operator = hfs;
  options{
    // magnetisation distribution model:
    F = ball;
    // magnetic rms radius:
    rrms = 7.0; // by default, will be same as charge rms
  }

  // Note: rpa (many-body corrections)
  // will run for muonic atoms, but is completely meaningless
  rpa = false;

  off-diagonal = false;
  diagonal = true;
}
```

**NOTE** always set rpa=false (or ignore the rpa output) for muonic atoms.
RPA is a many-body correction. The code calculates this for muons, even though it shouldn't. It's meaningless for muons, and will not converge. Results mean nothing.

There's als a module specificly for hyperfine anomaly.
It just calculates matrix elements twice, one for pointlike, and takes ratio.

```java
Module::hfAnomaly {
  // same options as hfs operator
  hfs_options{
    F = ball;
    rrms = 7.0;
  }
}
```
