#pragma once
#include <iostream>
#include <string>

#define XSTRING(s) STRING(s)
#define STRING(s) #s

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
namespace GitInfo {

const std::string git_branch = XSTRING(GITBRANCH);
const std::string git_revision = XSTRING(GITREVISION);
const std::string git_modified = XSTRING(GITMODIFIED);
const std::string cxx_version = XSTRING(CXXVERSION);
const std::string compiled_time = XSTRING(COMPTIME);

void print_git_info() {
  std::cout << "AMPSCI v: " << GitInfo::git_branch << "/"
            << GitInfo::git_revision << "\n";
  if (!GitInfo::git_modified.empty())
    std::cout << "Modified: " << GitInfo::git_modified << "\n";
  std::cout << "Compiled: " << GitInfo::cxx_version << "\n          "
            << GitInfo::compiled_time << "\n";
}

} // namespace GitInfo
