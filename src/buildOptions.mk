# BUILD OPTIONS:

## c++ standard. must be at least c++17
CXXSTD ?= -std=c++17

## Set directories for source files (SRC), and output built object files (BUILD)
SRC = ./src
BUILD ?= ./build

################################################################################

## Parallelism for build. Default to JOBS if -j not passed
## Directly -j from command line seems to override this, but maybe not gaurenteed
MAKEFLAGS += $(DEFAULT_MFLAGS)

################################################################################
## Clean variables (strip spaces) for later use

OMPLIB  := $(strip $(OMPLIB))
MODE    := $(strip $(MODE))
SRC     := $(strip $(SRC))
BUILD   := $(strip $(BUILD))
OPT     := $(strip $(OPT))
CXXNAME := $(notdir $(word 1,$(CXX)))
GSL_PATH:= $(strip $(GSL_PATH))

################################################################################
## Compiler warning flags (set base, then append per-compiler)

BASE_WARN = -Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Wshadow -Weffc++ -Wsign-conversion -Wno-psabi

GCC_WARN = -Wsuggest-override -Wnon-virtual-dtor -Wcast-align -Woverloaded-virtual -Wduplicated-cond -Wduplicated-branches -Wuseless-cast -Wformat

CLANG_WARN = -Wheader-hygiene -Wno-unused-function

## nb: must check for clang++ first, since 'clang++' contains 'g++'!
WARN = $(BASE_WARN)
ifneq (,$(findstring clang++,$(CXXNAME)))
  WARN += $(CLANG_WARN)
else ifneq (,$(findstring g++,$(CXXNAME)))
  WARN += $(GCC_WARN)
endif

################################################################################
## Build-mode switches:

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
## Include, compiler, and linker flags

INCLUDES = -I$(SRC)
ifneq ($(GSL_PATH),)
  INCLUDES += -I$(GSL_PATH)/include/
  LIBS += -L$(GSL_PATH)/lib/
endif

CXXFLAGS = $(CXXSTD) $(OPT) $(OMPLIB) $(WARN) $(INCLUDES)
LDFLAGS = $(OMPLIB)

# Compile and link commands
COMPILE = $(CXX) $(CARGS) -MMD -MP -c -o $@ $< $(CXXFLAGS)
LINK = $(CXX) $(LARGS) -o $@ $^ $(LDFLAGS) $(LIBS)

################################################################################
## Git info for version details

## Build / Git metadata (collected once)
NOW := $(shell date +%Y-%m-%d' '%H:%M' '%Z 2>/dev/null)
GITREVISION := $(shell git rev-parse --short HEAD 2>/dev/null)
GITBRANCH   := $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)
GITMODIFIED := $(shell git status -s 2>/dev/null)
CXXVERSION  := $(shell $(CXX) --version 2>/dev/null | sed -n '1p' | sed 's/(/[/; s/)/]/')
CXXPATH := $(shell which $(word 1,$(CXX)) 2>/dev/null)

## Preprocessor flags (passed to C++ compiler)
GITFLAGS =
GITFLAGS += -D GITREVISION="$(GITREVISION)"
GITFLAGS += -D GITBRANCH="$(GITBRANCH)"
GITFLAGS += -D GITMODIFIED="$(GITMODIFIED)"
GITFLAGS += -D CXXVERSION="$(CXXVERSION)"
GITFLAGS += -D COMPTIME="$(NOW)"
