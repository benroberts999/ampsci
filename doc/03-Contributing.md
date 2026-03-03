# Contributing

Updates, additions, changes are welcome.
Please follow the guidelines below to help keep the codebase clean, consistent, and maintainable, and to help make merging your changes in easy.

## Development Setup

* Ensure the build mode is set to **development** (`dev`):
  * Set `MODE=dev` in the `Makefile`.
* This turns on all the warnings.
* While these can be very annoying at first, they are all there for a good reason.
* Code shuold not flag any warnings to be merged into main.

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
  * Rebasing (`git rebase`) is preferred.
  * If you are new to rebasing, consider doing on a temporary branch first.
  * Ideally, if your changes are will targetted, there will not be conflicts, however, they do arise
  * Resolve conflicts carefully and test after rebasing.

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

## Code Formatting Style

* Consistent formatting minimises unnecessary git diffs and improves readability.
* Rather than doing this by hand, it is enforced by a tool [**clang-format**](https://clang.llvm.org/docs/ClangFormat.html), which automatically formats the code
* Install the tool:
  * `sudo apt install clang-format` and/or `brew install clang-format`
  * Or the extension in VsCode
* It's recommended to set up the VSCode extension to format the code automatically.
* A style file, which defines the style rules, is provided in the repository `src/.clang-format`.
* The makefile also has a target to directly run clang-format, which is a good idea before committing:

<div class="shell-block">
```bash
make clang_format
```
</div>

## Documenting

* All new functions/classes should be documented using Doxygen-style comments
* Usually, functions and parameters should be names so that their meaning is obvious. If not, add the description to the `details`
* The `breif` first line should be ~1 line, use no special symbols, and succinctly state what the function/class does; this is what gets shown to users in IDE/editor on hover
* The
  * e.g.,

```cpp

//! Takes an x and calculates something with it
/*!
  @details
  A more detailed description, that may contain LaTeX (using '\f' that doxygen expectes)
  \f[ 
    y = mx + b 
  \f]
  or \f$ y = mx + b \f$ for inline equations.
   - detailed description of input x, what it is, what it must be
   - detailed description of what output expected
*/
double my_new_function(double x);
```

## Testing

* Where possible, add new unit/integration tests for new functionality.
  * Tests are strongly encouraged for core library changes.
  * Tests are not required for standalone modules.
  * Code uses `Catch2` style tests - easiest way is to copy existing tests
  * All tests must be in file that ends with
    * `.tests.cpp`
    * This ensures they are compiled correctly into the right executable
* More importantly, ensure existing tests pass before submitting a pull request.
  * In dev mode, `tests` will be compiled by default
  * Otherwise, compile the tests executable

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
