\page contributing_docs Documenting

\brief Writing Doxygen documentation for ampsci.

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
