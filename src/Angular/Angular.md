# Angular

## Angular::Ck_ab (class)

 * Lookup table for C^k and 3j symbols (special m=1/2 case)
 * Builds 3j symbol lookup table for given maximum k and maximum j (2j)
 * 3j symbols, special case: (ja jb k \\ -1/2, 1/2, 0)
 * Ck_ab      := <ka||C^k||kb>        [symmetric up to +/- sign]
 * TildeCk_ab := (-1)^{ja+1/2} * ckab [symmetric]
 * Slightly faster than calculating on-the-fly

```cpp
Ck_ab(const int in_max_K = 0, const int in_max_twoj = 0);
void fill_maxK_twojmax(const int in_max_K, const int in_max_twoj);
```
 * Constructor: fills 3j table with all non-zero 3j-symbols up to and including the maximum specified k and j values (entered as integer 2j)
 * 'fill' function will extend the lookup table to new maximum values

```cpp
double get_tildeCkab_mutable(int k, int ka, int kb);
double get_Ckab_mutable(int k, int ka, int kb);
double get_3jkab_mutable(int k, int ka, int kb);
double get_tildeCkab(int k, int ka, int kb) const;
double get_Ckab(int k, int ka, int kb) const;
double get_3jkab(int k, int ka, int kb) const;
```
 * The _getters_ -- lookup values from table + return them
 * Note: input is _kappa_ (not j or 2j) [Ck includes parity]
 * The _mutable_ versions will calculate (+store) the required 3j symbols if user asks for one that doesn't exist (up to new max_2j and max_K)
   * This is 'safer' (in that it won't segfault), but is slower and not thread-safe [calls std::vector::push_back]
 * not-mutable versions are faster, and thread-safe (read-only) - but perform no bounds checking, so will seg-fault if you ask for a symbol that isn't calculated.

## Angular::SixJ (class)

 * Lookup table for 6j symbols: { ja, jb, k \\ jc, jd, l}
 * j's half-integer (called using integer 2j). Integer k and l
 * Much faster than calculating on the fly.
 * Stores 6j symbols up to + including given max_K and max_2j (stores all allowed l)

```cpp
SixJ(int in_max_k, int in_max_twoj);
void fill(int in_max_k, int in_max_twoj);
```
 * Constructor: fills 6j table with all non-zero 6j-symbols up to and including the maximum specified k and j values (entered as integer 2j)
 * 'fill' function will extend the lookup table to new maximum values

```cpp
double get_6j_mutable(int tja, int tjb, int tjc, int tjd, int k, int l);
double get_6j(int tja, int tjb, int tjc, int tjd, int k, int l) const;
```
 * The _getters_ -- lookup values from table + return them
 * Note: input is integer _2j_ (not j or kappa)
 * The _mutable_ versions will calculate (+store) the required 6j symbols if user asks for one that doesn't exist (up to new max_2j and max_K)
   * This is 'safer' (in that it won't segfault), but is slower and not thread-safe [calls std::vector::push_back]
 * not-mutable versions are faster, and thread-safe (read-only) - but perform no bounds checking, so will seg-fault if you ask for a symbol that isn't calculated.
