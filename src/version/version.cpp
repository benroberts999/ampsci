#include "version.hpp"
#include <string>

// Macro translates constants to "strings"
#define XSTRING(s) STRING(s)
#define STRING(s) #s

// Constants refer to git revision info.
// These are passed in at compile time via -D flag
// If not set, the code will still work (these will just be blank)
// These are in the .cpp file, so we only need to re-build this file (and
// re-link) whenever we want updated git version info [most applicable for the
// GITMODIFIED option, which is not easy to implement otherwise]
#ifndef GITBRANCH
#define GITBRANCH
#endif
#ifndef GITREVISION
#define GITREVISION
#endif
#ifndef GITMODIFIED
#define GITMODIFIED
#endif
#ifndef CXXVERSION
#define CXXVERSION
#endif
#ifndef COMPTIME
#define COMPTIME
#endif

//==============================================================================
namespace version {

static const std::string git_branch = XSTRING(GITBRANCH);
static const std::string git_revision = XSTRING(GITREVISION);
static const std::string git_modified = XSTRING(GITMODIFIED);
static const std::string cxx_version = XSTRING(CXXVERSION);
static const std::string compiled_time = XSTRING(COMPTIME);
static const std::string ampsci_version = XSTRING(AMPSCI_VERSION);

std::string version() {
  return git_revision.empty() ?
             ampsci_version :
         git_modified.empty() ?
             ampsci_version + " [" + git_branch + "/" + git_revision + "]" :
             ampsci_version + " [" + git_branch + "/" + git_revision + "]*\n" +
                 " *(Modified: " + git_modified + ")";
}

std::string compiled() { return cxx_version + " " + compiled_time; }

} // namespace version
