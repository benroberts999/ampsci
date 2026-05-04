# f-dependent Breit

## Fix interface:
* First, user options
* Then interface through HF

## HF::Breit
* Make one Breit integral function
* Then, it will calculate static or f-dep breit based on option
* That way: only need to make change at top level of code, then HF::Breit takes care of the rest!

## Everywhere
* omega and w are used as variable names, when alpha*omega (q, momentum exchange) are passed.
* This is begging for mistakes
* In most similar places, ampsci passes the frequency, not the momentum, so:
  * Ideally, omega (not alpha*omega) should be passed to all Breit functions
  * then when alpha*omega is used inside functions, rename to 'qw' or similar
  * `q*r` passed to Spherical bessel, `q = alpha*w`
  * Except, that would require bein *very* careful everywhere - maybe better to change to qw
* Which frequency should use used for exchange integrals P, W ?

## Failing unit tests

* For the raw Breit/Coulomb integrals, some current tests fail
* This might not actually be wrong, but should be checked
* Odd behaviour for large w
  * "turn-around"

```
gk_ab:
k = 0 : gk_ab
    r      w=1e+01      w=5e+00      w=5e-01      w=5e-02      w=1e-03      w=1e-09       static
1e-06  -6.7901e-02  -7.8807e-02  -6.6484e-02  -6.6266e-02  -6.6264e-02  -6.6264e-02  -6.6264e-02 +1.3e-08
1e-02  -6.7901e-02  -7.8807e-02  -6.6484e-02  -6.6266e-02  -6.6264e-02  -6.6264e-02  -6.6264e-02 +1.3e-08
1e-01  -6.7901e-02  -7.8807e-02  -6.6484e-02  -6.6266e-02  -6.6264e-02  -6.6264e-02  -6.6264e-02 +1.3e-08
5e+00  -6.7900e-02  -7.8807e-02  -6.6484e-02  -6.6266e-02  -6.6264e-02  -6.6264e-02  -6.6264e-02 +1.3e-08
1e+02  -6.7900e-02  -7.8807e-02  -6.6484e-02  -6.6266e-02  -6.6264e-02  -6.6264e-02  -6.6264e-02 +1.3e-08
```

* Even odder for the 'vk_ab' case, e.g.,
* Again, might not be wrong, but should be understood

```
k = 3 [vk_ab]
    r      w=1e+01      w=5e+00      w=5e-01      w=5e-02      w=1e-03      w=1e-09
1e-04  -1.3844e-06  -1.3768e-06  -1.6965e-06  -3.3945e-05  -1.3708e-06  -1.3708e-06
1e-02  -3.3940e-03  -3.2926e-03  -3.2647e-03  -3.2725e-03  -3.2643e-03  -3.2643e-03
1e-01  -1.8971e-02  -1.4843e-02  -1.3801e-02  -1.3792e-02  -1.3792e-02  -1.3792e-02
5e+00   3.6032e-04   6.3125e-04   2.3847e-05   1.5970e-05   1.5917e-05   1.5917e-05
1e+02   1.0074e-05  -4.3840e-05   7.2995e-07   5.5782e-09   2.1890e-09   2.1890e-09
```