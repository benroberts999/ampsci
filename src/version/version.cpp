#include "version.hpp"
#include <dlfcn.h>
#include <gsl/gsl_version.h>
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
#ifndef EXTERNAL_MODULES
#define EXTERNAL_MODULES
#endif
// GSL_VERSION defined by GSL library (if it exists)
#ifndef GSL_VERSION
#define GSL_VERSION ""
#endif
// _OPENMP defined by openMP library (if it exists)
#ifdef _OPENMP
#define OMP_VERSION _OPENMP
#else
#define OMP_VERSION 0
#endif

//==============================================================================
namespace version {

// git branch for current compilation (if available)
static const std::string git_branch = XSTRING(GITBRANCH);
// git hash for current compilation (if available)
static const std::string git_revision = XSTRING(GITREVISION);
// List of files that have been modified since last git commit
static const std::string git_modified = XSTRING(GITMODIFIED);
// Compiler version information
static const std::string cxx_version = XSTRING(CXXVERSION);
// Date and time of compilation
static const std::string compiled_time = XSTRING(COMPTIME);
// External modules compiled in (space-separated filenames, may be empty)
static const std::string external_modules = XSTRING(EXTERNAL_MODULES);
// ampsci version
static const std::string ampsci_version = std::string(AMPSCI_VERSION);
// gsl library version
static const std::string gsl_version = std::string(GSL_VERSION);
// OpenMP library version
static const std::string omp_version = XSTRING(OMP_VERSION);

std::string version() {
  std::string v = git_revision.empty() ? ampsci_version :
                  git_modified.empty() ? ampsci_version + " [" + git_branch +
                                           "/" + git_revision + "]" :
                                         ampsci_version + " [" + git_branch +
                                           "/" + git_revision + "]*\n" +
                                           " *(Modified: " + git_modified + ")";
  if (!external_modules.empty())
    v += "\n External modules: " + external_modules;
  return v;
}

std::string compiled() { return cxx_version + " " + compiled_time; }

static std::string blas_info() {
  using fn_str = const char *(*)();
  using fn_mkl = void (*)(char *, int);
  auto oblas_cfg =
    reinterpret_cast<fn_str>(dlsym(RTLD_DEFAULT, "openblas_get_config"));
  if (oblas_cfg)
    return std::string("OpenBLAS: ") + oblas_cfg();
  auto mkl_ver =
    reinterpret_cast<fn_mkl>(dlsym(RTLD_DEFAULT, "MKL_Get_Version_String"));
  if (mkl_ver) {
    char buf[256] = {};
    mkl_ver(buf, 256);
    return std::string("Intel MKL: ") + buf;
  }
  return "Reference LAPACK/BLAS";
}

std::string blas_threads() {
  using fn_int = int (*)();
  auto oblas_thr =
    reinterpret_cast<fn_int>(dlsym(RTLD_DEFAULT, "openblas_get_num_threads"));
  if (oblas_thr)
    return "OpenBLAS: " + std::to_string(oblas_thr()) + " threads.";
  auto mkl_thr =
    reinterpret_cast<fn_int>(dlsym(RTLD_DEFAULT, "mkl_get_max_threads"));
  if (mkl_thr)
    return "MKL: " + std::to_string(mkl_thr()) + " threads.";
  return "";
}

std::string libraries() {
  return "  GSL (GNU Scientific Libraries): " + gsl_version + '\n' +
         "  OpenMP: " + omp_version + '\n' + "  LAPACK/BLAS: " + blas_info();
}

} // namespace version
