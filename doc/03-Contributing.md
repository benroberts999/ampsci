# Contributing

Updates, additions, changes are welcome.
Please follow the guidelines below to help keep the codebase clean, consistent, and maintainable, and to help make merging your changes in easy.

## Development Setup

* Ensure the build mode is set to **development** (`dev`):
  * Set `MODE=dev` in the `Makefile`.
* This turns on all the warnings.
* While these can be very annoying at first, they are all there for a good reason.
* Code should not flag any warnings to be merged into main.

## General Principles

* The most common case is to simply add a new module (see [Modules](\ref modules))
  * Since all modules are independent, nothing you do here can impact anything else, so feel free to add as much as you like to your new module file!
* You may also need to add new functionality, or extend existing methods.
* Prefer minimal and targeted changes; prefer _adding_ rather than making changes to parts of the code
  * Avoid modifying existing code unless necessary, and keep any changes you do have to make small and local if possible to avoid conflicting changes with other users
  * If you make "temporary" changes to other parts of the code (e.g., for testing), ensure this is clearly marked and reverted before committing your changes
    * It's very easy to accidentally change functionality that might not be intended or noticed by others.
* With ease of use as a high priority, we prefer not to rely on too many external libraries, which can be painful to install/link, particularly on older systems

## Git Workflow

* Branch from the appropriate parent branch (usually `dev`).
* Create a new branch with a clear, descriptive name.
* It's a good idea to regularly incorporate updates from the parent branch, which will likely also be undergoing continuing changes
  * Rebasing (`git rebase dev`) is preferred.
  * If you are new to rebasing, consider doing on a temporary branch first.
  * Ideally, if your changes are well targetted, there will not be conflicts, however, they do arise
  * Resolve conflicts carefully and test after rebasing.
  * It is OK to submit a pull request without rebasing from dev if it's too difficuly; but that just means I'll have to do it. If it's too hard to resolve all conflicts, the pull request may be rejected

<div class="shell-block">
```bash
git rebase dev
```
</div>

* In theory, it's a good idea to have a nice clear commit history
* At the same time, too many commits is better than too few!
* Try to commit small changes at a time, rather than one giant commit
  * Each commit should represent one logical change, making it easier to review, revert, or bisect if needed.
  * Commits can always be squashed together at the end
* Write clear, descriptive commit messages, use a concise imperative title (e.g. "Adds such and such function"), explaining what the change does.
  * Use detailed description if required, putting as much information as required

## Code Formatting Style: clang-format

* Consistent formatting minimises unnecessary git diffs and improves readability.
* Rather than doing this by hand, it is enforced by a tool [**clang-format**](https://clang.llvm.org/docs/ClangFormat.html), which automatically formats the code
* Install the tool (for example):
  * Linux: `sudo apt install clang-format`
  * macOS: `brew install clang-format`
* A style file, which defines the style rules, is provided in the repository `src/.clang-format`.
  * Unfortunately, different versions of clang-format make slightly different choices. v16+ is ideal, but the differences are usually small
* The makefile has a target to directly run clang-format. It also normalises a few common differences between different clang-format version. It will only format files you have changed since last git commit

<div class="shell-block">
```bash
make clang_format
```
</div>

### VSCode with clang-format

The VSCode C/C++ extension ships its own bundled formatting tool, which may differ in version from the system-installed one and produce inconsistent formatting.
To avoid this, configure VSCode to use the system clang-format via `.vscode/settings.json` (this file is not tracked by git -- VSCode may have already created it, otherwise create it in the project directoty).

First, find the path to your clang-format binary:

<div class="shell-block">
```bash
which clang-format
```
</div>

Common locations:

* Linux: `/usr/bin/clang-format`
* macOS (Homebrew): `/opt/homebrew/bin/clang-format`

Then add the following to `.vscode/settings.json`:

```json
{
  "C_Cpp.clang_format_path": "/PATH_TO_CLANG_FORMAT/clang-format",
  "editor.formatOnSave": true
}
```

Replace `PATH_TO_CLANG_FORMAT` with the path returned from `which clang-format`

--------------------

## Documenting

* All new functions/classes should be documented using Doxygen-style comments
* Usually, functions and parameters should be names so that their meaning is obvious.
  * If additional explanation is needed, add it in the `@details` section.
* The `@brief` line should be about one sentence. Avoid special symbols and
  succinctly describe what the function or class does.  
  * This text is what appears in IDE/editor tooltips when hovering over the symbol.
* LaTeX is supported, but uses Doxygen's formula delimiters instead of standard LaTeX:
  display equations use `\f[` ... `\f]` (not `\[` ... `\]`), and inline equations use `\f$` ... `\f$` (not `$` ... `$`)

Example:

@include comment_example.cpp

### Common Doxygen commands

| Command | Purpose |
|-------|-------|
| `@brief` | Short summary of the function/class |
| `@details` | Longer description, often including equations |
| `@param name` | Describe a function parameter (often not needed if named well) |
| `@tparam name` | Describe a template parameter (often not needed if named well) |
| `@return` | Describe the return value (often not needed if named well) |
| `@note` | Additional information |
| `@warning` | Important warning for users |
| `@see` | Reference related functions or classes |
| `@ref` | Link to another documented symbol |

## Testing

* Where possible, add new unit/integration tests for new functionality.
  * Tests are strongly encouraged for core library changes.
  * Tests are not required for standalone modules.
  * Code uses `Catch2` style tests - easiest way is to copy existing tests, or see example below
  * All tests must be in a cpp file that ends with
    * `.tests.cpp`
    * This ensures they are compiled correctly into the right executable
* More importantly, ensure existing tests pass before submitting a pull request.
  * In dev mode, `tests` will be compiled by default
  * Otherwise, compile the tests executable
* Running tests will produce junk output files. These all have `deleteme` in their filenames, so can be cleared easily. The Makefile has a target for this: `make remove_junk`

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

* `[unit]` — tests a small component and runs very quickly.
* `[integration]` — checks numerical accuracy and interaction between parts of the code.
* `[slow]` — long-running tests (≈1 minute or more).

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
