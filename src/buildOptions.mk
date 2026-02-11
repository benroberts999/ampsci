# BUILD OPTIONS:

## Set directories for source files (SRC), and output built object files (BUILD)
SRC=./src
BUILD=./build

################################################################################
# Optimisation flags
OPT?=-O3

################################################################################
# Library settings
LIBS ?= -lgsl -lgslcblas -llapack -lblas

################################################################################
# Clean variables for later use

OMPLIB  := $(strip $(OMPLIB))
MODE    := $(strip $(MODE))
SRC     := $(strip $(SRC))
BUILD   := $(strip $(BUILD))
OPT     := $(strip $(OPT))
CXXNAME := $(notdir $(word 1,$(CXX)))

################################################################################
# Compiler warning flags (set base, then append per-compiler)

BASE_WARN = -Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Wshadow -Weffc++ -Wsign-conversion -Wno-psabi

GCC_WARN = -Wsuggest-override -Wnon-virtual-dtor -Wcast-align -Woverloaded-virtual -Wduplicated-cond -Wduplicated-branches -Wuseless-cast -Wformat

CLANG_WARN = -Wheader-hygiene -Wno-unused-function

# nb: must check for clang++ first, since clang++ contains g++!
WARN = $(BASE_WARN)
ifneq (,$(findstring clang++,$(CXXNAME)))
  WARN += $(CLANG_WARN)
else ifneq (,$(findstring g++,$(CXXNAME)))
  WARN += $(GCC_WARN)
endif

################################################################################
# Build-mode switches:

ifeq ($(MODE),release)
  WARN = -w -Wno-psabi
  CARGS += -g0 -DNDEBUG -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
endif
ifeq ($(MODE),debug)
  OMPLIB=
  OPT = -O0
  CARGS += -g3 -fno-omit-frame-pointer
endif

# suppress unknown OpenMP pragmas when OpenMP is not enabled.
ifeq ($(OMPLIB),)
  WARN += -Wno-unknown-pragmas
endif


################################################################################
# Include and linker flags
INCLUDES = -I$(SRC)
ifneq ($(PathForGSL),)
  INCLUDES += -I$(PathForGSL)/include/
  LIBS += -L$(PathForGSL)/lib/
endif

CXXFLAGS = $(CXXSTD) $(OPT) $(OMPLIB) $(WARN) $(INCLUDES)
LDFLAGS = $(OMPLIB)

# Compile and link commands
COMPILE = $(CXX) $(CARGS) -MMD -MP -c -o $@ $< $(CXXFLAGS)
LINK = $(CXX) $(LARGS) -o $@ $^ $(LDFLAGS) $(LIBS)
