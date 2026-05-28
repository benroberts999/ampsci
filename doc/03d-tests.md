\page contributing_tests Testing

\brief Writing and running tests for ampsci using Catch2.

* Where possible, add new unit/integration tests for any new functionality that you add.
  * Tests are strongly encouraged for core ampsci changes.
  * Tests are not required for standalone modules.
  * Code uses `Catch2` style tests - easiest way is to copy existing tests, or see example below
  * All tests must be in a cpp file that ends with
    * `.tests.cpp`
    * This ensures they are compiled correctly into the right executable
* More importantly, ensure all existing tests pass before submitting a pull request.
  * In dev mode, `tests` will be compiled by default
  * Otherwise, compile the tests executable

## Running tests

* Running tests will produce junk output files. These all have `deleteme` in their filenames, so can be cleared easily. The Makefile has a target for this: `make remove_junk`
* The full set of tests can take a very long time; though running just the unit tests should be very quick.
 * The set of tests that are run can be selected by passing "filters" to the commandline `tests` executable (see examples below)


In dev mode, `tests` will be compiled by default.
  * Otherwise, compile the tests executable:
<div class="shell-block">
```bash
make tests
```
</div>

* The following will run, respectively:
  * Just the 'unit' tests (very quick, but doesn't check numerical accuracy)
  * All test tests, except the slowest ones (will run for several minutes)
  * All the tests (may take >30 minutes)

<div class="shell-block">
```bash
./tests [unit]
```
</div>
<div class="shell-block">
```bash
./tests ~[slow]
```
</div>
<div class="shell-block">
```bash
./tests
```
</div>

Then, remove the junk output files (all have `deleteme` in the filename):

<div class="shell-block">
```bash
make remove_junk
```
</div>

@note
On mac, you may need to escape the filters:
<div class="shell-block">
```bash
./tests "~[slow]"
```
</div>


### Writing tests using Catch2

* Try to write **small, focused tests** that check a single behaviour.
* Give tests **clear descriptive names** so failures are easy to understand.
* Use `REQUIRE` for conditions that must hold; use `CHECK` when later checks should still run.
* Try to test **edge cases** as well as normal inputs (e.g. zero, limits, invalid values).
* Keep tests **deterministic and fast**, avoiding randomness unless explicitly controlled.
* **NOTE** ampsci often produces output files. This is annoying when running many tests; it also can "damage" the integrity of the tests, since sometimes these files are read in too!
  * As a general rule, unless you are testing the read/write, we should attempt to not write data files during tests
  * If they must be written: (a) including a random string and the phrase `deleteme` in the output filename.
  * This can usually be achieved by adding this string to the "run label"
    * The random string stops this file from being read in by a subsequent test if that's not what we want, and `deleteme` makes it easy to clear all files using a sript

```cpp
const auto label = "deleteme_" + qip::random_string(4);
```

Examples:

```cpp
TEST_CASE("a name for the test", "[tag]") {
    // test body
}
```

* The **test name** should clearly describe the behaviour being tested.
  * Write names as **short sentences describing expected behaviour**, not implementation details.
* **Tags** (the strings in `[...]`) group related tests so they can be run selectively.
  * A test may have **multiple tags**, e.g. `[math][fast]`.

Common tags used in this project:

* `[unit]` -- tests a small component and runs very quickly.
* `[integration]` -- checks numerical accuracy and interaction between parts of the code.
* `[slow]` -- long-running tests (≈1 minute or more).

* If a test is checking numerical accuracy and how the code interfaces with the rest of ampsci, mark it '[integration]'
* If a test is for one small component and runs very quickly, mark it `[unit]`
* If a test is slow (>~1 min), mark it as `[slow]`
* Then, give at least one tag for the catagory being tested, e.g. `[MBPT]`
  * usually, this will match the directory the file is in
  * May often have multiple tags, e.g., `[MBPT][RPA][Basis]`

Example:

```cpp
#include "catch2/catch.hpp"
#include "math/square.hpp" // example, not real code

TEST_CASE("square computes correct values", "[math][unit]") {
    REQUIRE(square(2) == 4);
    REQUIRE(square(-2) == 4);
    REQUIRE(square(3) == 9);
}

TEST_CASE("example showing Approx for floating point", "[math][unit]") {

    const double y = square(2.0);

    // Default tolerance (~machine precision)
    REQUIRE(y == Catch::Approx(4.0));

    // Relative tolerance (epsilon)
    REQUIRE(y == Catch::Approx(4.0).epsilon(1e-12));

    // Absolute tolerance (delta)
    REQUIRE(y == Catch::Approx(4.0).margin(1e-10));
}

TEST_CASE("square handles zero input", "[math][unit]") {
    REQUIRE(square(0) == 0);
    REQUIRE(square(0ul) == 0ul);
    REQUIRE(square(0.0f) == 0.0f);
    REQUIRE(square(0.0) == 0.0);
}
```
