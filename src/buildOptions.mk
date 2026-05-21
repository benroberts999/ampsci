# BUILD OPTIONS:

## c++ standard. must be at least c++17
CXXSTD ?= -std=c++17

## Set directories for source files (SRC), and output built object files (BUILD)
SRC = ./src
BUILD ?= ./build

################################################################################
## External modules: list extra .cpp files to compile and link into ampsci/tests.
## Paths may be absolute or relative to the directory where make is invoked.
## Globs are supported (e.g. ext/*.cpp). Examples:
##   EXTERNAL_MODULES = /abs/path/MyModule.cpp ../shared/OtherModule.cpp
##   EXTERNAL_MODULES = ext/*.cpp
## Objects land in $(BUILD_DIR)/ExternalModules/. Each file must have a unique
## basename (two files named MyMod.cpp from different directories will collide).
## Paths with spaces are not supported. Recompile after changing this list.
EXTERNAL_MODULES ?=

## External operators: same as EXTERNAL_MODULES but for custom DiracOperator files.
## Usage is identical -- list .cpp files containing your operator class + Register<T>.
##   EXTERNAL_OPERATORS = /abs/path/MyOperator.cpp
##   EXTERNAL_OPERATORS = ext/operators/*.cpp
EXTERNAL_OPERATORS ?=

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

BASE_WARN = \
  -Wall -Wextra -Wpedantic \
  -Wdouble-promotion -Wconversion -Wsign-conversion -Wsign-compare \
  -Wshadow -Wunused-parameter -Weffc++ \
  -Wzero-as-null-pointer-constant -Wnonnull -Wdeprecated -Wno-psabi \
  -Wimplicit-fallthrough -Wvla

# add soon:  -Wold-style-cast -Wswitch-enum  -Wswitch-default

GCC_WARN = \
  -Wsuggest-override -Wnon-virtual-dtor -Woverloaded-virtual \
  -Wduplicated-cond -Wduplicated-branches \
  -Wcast-align -Wuseless-cast -Wformat \
  -Wpessimizing-move -Wstrict-aliasing -Wmaybe-uninitialized -Wlogical-op

CLANG_WARN = \
  -Wheader-hygiene -Wno-unused-function \
  -Wno-unknown-warning-option \
  -Wno-invalid-utf8 -Wno-c2x-extensions -Wno-c2y-extensions

## nb: must check for clang++ first, since 'clang++' contains 'g++'!
WARN = $(BASE_WARN)
ifneq (,$(findstring clang++,$(CXXNAME)))
  WARN += $(CLANG_WARN)
else ifneq (,$(findstring g++,$(CXXNAME)))
  WARN += $(GCC_WARN)
endif

# suppress unknown OpenMP pragmas when OpenMP is not enabled.
ifeq ($(OMPLIB),)
  WARN += -Wno-unknown-pragmas
endif

################################################################################
## Build-mode switches:

ifeq ($(MODE),release)
  WARN = -w -Wno-psabi -Wno-unknown-pragmas
  EXTRA_CXXFLAGS += -g0 -DNDEBUG -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
endif

ifeq ($(MODE),debug)
  OMPLIB=
  OPT = -O0
  EXTRA_CXXFLAGS += -g3 -fno-omit-frame-pointer
endif

## clang sanitisers:

ifeq ($(MODE),asan)
  # AddressSanitizer + UndefinedBehaviourSanitizer. Use clang++ for best output.
  # OMP disabled to avoid false positives from OMP runtime.
  OMPLIB=
  OPT = -O1
  EXTRA_CXXFLAGS += -g3 -fno-omit-frame-pointer -fsanitize=address,undefined
  LDFLAGS += -fsanitize=address,undefined
endif

ifeq ($(MODE),tsan)
  # ThreadSanitizer. Use clang++ for best output.
  OPT = -O1
  EXTRA_CXXFLAGS += -g3 -fno-omit-frame-pointer -fsanitize=thread
  LDFLAGS += -fsanitize=thread
endif

ifeq ($(MODE),ubsan)
  # UndefinedBehaviourSanitizer only: lower overhead than asan, OMP kept.
  OPT = -O1
  EXTRA_CXXFLAGS += -g3 -fno-omit-frame-pointer -fsanitize=undefined
  LDFLAGS += -fsanitize=undefined
endif

################################################################################
## Include, compiler, and linker flags

INCLUDES = -I$(SRC)

# Common one: GSL (other library includes may be added to CXXFLAGS and LDLIBS)
ifneq ($(GSL_PATH),)
  INCLUDES += -I$(GSL_PATH)/include/
  LDFLAGS += -L$(GSL_PATH)/lib/
endif

CXXFLAGS += $(CXXSTD) $(OPT) $(OMPLIB) $(WARN) $(INCLUDES) $(EXTRA_CXXFLAGS)
LDFLAGS += $(OMPLIB)

# Compile and link commands
COMPILE = $(CXX) $(CXXFLAGS) -MMD -MP -c -o $@ $<
LINK = $(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

################################################################################
## Git info for version details

## Build / Git metadata (collected once)
NOW := $(shell date +%Y-%m-%d' '%H:%M' '%Z 2>/dev/null)
GITREVISION := $(shell git rev-parse --short HEAD 2>/dev/null)
GITBRANCH   := $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)
GITMODIFIED := $(shell git diff --name-only HEAD 2>/dev/null)
CXXVERSION  := $(shell $(CXX) --version 2>/dev/null | sed -n '1p' | sed 's/(/[/; s/)/]/')
CXXPATH := $(shell which $(word 1,$(CXX)) 2>/dev/null)

## Preprocessor flags (passed to C++ compiler)
GITFLAGS =
GITFLAGS += -D GITREVISION="$(GITREVISION)"
GITFLAGS += -D GITBRANCH="$(GITBRANCH)"
GITFLAGS += -D GITMODIFIED="$(GITMODIFIED)"
GITFLAGS += -D CXXVERSION="$(CXXVERSION)"
GITFLAGS += -D COMPTIME="$(NOW)"
GITFLAGS += -D EXTERNAL_MODULES="$(foreach src,$(filter %.cpp,$(wildcard $(EXTERNAL_MODULES))),$(notdir $(src)))"
GITFLAGS += -D EXTERNAL_OPERATORS="$(foreach src,$(filter %.cpp,$(wildcard $(EXTERNAL_OPERATORS))),$(notdir $(src)))"
