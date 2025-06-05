################################################################################
# Options + settings for ampsci Makefile system.
# Included by the main Makefile. Defines flags, libraries, and build logic.
################################################################################

# Parallel build settings
ParallelBuild ?= 1
ifeq ($(Build),debug)
  ParallelBuild=1
endif
ifneq ($(ParallelBuild),1)
  MAKEFLAGS += -j$(ParallelBuild)
  $(info Compiling in parallel with N=$(ParallelBuild) threads)
endif

# Compiler warning flags (set base, then append per-compiler)
BASE_WARN = -Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Wshadow -Weffc++ -Wsign-conversion -Wno-psabi

GCC_WARN = -Wsuggest-override -Wnon-virtual-dtor -Wcast-align \
  -Woverloaded-virtual -Wduplicated-cond -Wduplicated-branches \
  -Wuseless-cast -Wformat

# -Wnull-dereference --> have to remove with g++-13

CLANG_WARN = -Wheader-hygiene -Wno-unused-function

# Special case for g++-7
ifeq ($(CXX),g++-7)
  WARN = -Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Weffc++ \
    -Wsign-conversion -Wno-unused-variable -Wno-shadow
else
  WARN = $(BASE_WARN)
  ifneq (,$(findstring clang++,$(CXX)))
    WARN += $(CLANG_WARN)
  else ifneq (,$(findstring g++,$(CXX)))
    WARN += $(GCC_WARN)
    ifeq ($(CXX),g++-13)
      WARN += -Wno-null-dereference
    endif
  endif
endif


# Optimisation flags
OPT ?= -O3

ifeq ($(Build),release)
  WARN = -w -Wno-psabi
  OPT += -g0 -DNDEBUG -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
endif
ifeq ($(Build),debug)
  UseOpenMP = no
  WARN += -Wno-unknown-pragmas
  OPT = -O0 -g3 -fno-omit-frame-pointer
endif
# -fno-exceptions -fno-rtti not needed
# -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF are from GSL

# Normalize UseOpenMP to lowercase for consistent checks
USEOPENMP_LC := $(shell echo "$(UseOpenMP)" | tr '[:upper:]' '[:lower:]')

OMPLIB ?= -fopenmp
OMP ?= $(OMPLIB)
ifneq ($(USEOPENMP_LC),true)
ifneq ($(USEOPENMP_LC),yes)
  OMP =
  WARN += -Wno-unknown-pragmas
endif
endif

################################################################################
# Library settings
LIBS ?= -lgsl -lgslcblas -llapack -lblas

################################################################################
# Include and linker flags
INCLUDES = -I$(SD)
ifneq ($(PathForGSL),)
  INCLUDES += -I$(PathForGSL)/include/
  LIBS += -L$(PathForGSL)/lib/
endif

CXXFLAGS = $(CXXSTD) $(OPT) $(OMP) $(WARN) $(INCLUDES)
LDFLAGS = $(CXXSTD) $(OPT) $(OMP)

# Optional: user can add sanitizer flags via SANITIZERS variable
# SANITIZERS ?=
# CXXFLAGS += $(SANITIZERS)
# LDFLAGS += $(SANITIZERS)

# These should be used with clang in debug mode only
# MSAN = -fsanitize=memory
# ASAN = -fsanitize=address -fsanitize=leak
# TSAN = -fsanitize=thread
# USAN = -fsanitize=undefined -fsanitize=unsigned-integer-overflow
# ifeq ($(Build)+$(CXX),debug+clang++)
# 	CXXFLAGS += -g $(ASAN) -fno-omit-frame-pointer
# endif
# echo $(Build)+$(CXX)
# MSAN_SYMBOLIZER_PATH=/usr/lib/llvm-6.0/bin/llvm-symbolizer ./ampsci
# ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer ./ampsci

# Compile and link commands
COMP = $(CXX) $(CARGS) -MMD -MP -c -o $@ $< $(CXXFLAGS)
LINK = $(CXX) $(LARGS) -o $@ $^ $(LDFLAGS) $(LIBS)
