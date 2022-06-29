#pragma once
#include <string>

// Manually define major/minor ampsci versions
#define AMPSCI_VERSION 0.0
#define AMPSCI_MAJOR_VERSION 0
#define AMPSCI_MINOR_VERSION 0

//==============================================================================
//! Information about the ampsci code (version, compiler etc.)
namespace version {

//! String with version info, including git branch/revision, and if any files
//! have been modified since the last commit
std::string version();

//! String with compilation info, including which compiler was used and the time
//! of compilation
std::string compiled();

} // namespace version
