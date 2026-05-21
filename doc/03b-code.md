\page contributing_code Code Contributions

\brief Pull requests, development workflow, git, and code formatting.

If you plan for your code to be merged into the main ampsci program, the following guidelines will be helpful.

* See also: [Documenting](\ref contributing_docs)
* See also: [Testing](\ref contributing_tests)

## Pull requests

* Open a pull request on [GitHub](https://github.com/benroberts999/ampsci/pulls), ideally against the `dev` branch.
* Ensure all existing tests pass and clang-format has been run before submitting.
* A short description of what the PR does and why is sufficient; link to a related issue if one exists.
* Try to avoid conflicts with the rest of the code -- see below.

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
