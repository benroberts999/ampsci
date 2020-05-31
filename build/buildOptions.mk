# Options + settings for makefile for each 'Main'

## will return the Operating system name
detected_OS := $(shell uname -s)
$(info )
$(info Detected operating system: $(detected_OS))


# runs make in //
ifeq ($(Build),debug)
  ParallelBuild=0
endif
ifneq ($(ParallelBuild),)
ifneq ($(ParallelBuild),0)
ifneq ($(ParallelBuild),1)
  MAKEFLAGS += -j$(ParallelBuild)
  $(info Compiling in parallel with N=$(ParallelBuild) cores)
endif
endif
endif
$(info )

#Warnings:
WARN=-Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Wshadow -Weffc++ -Wsign-conversion
# -Wuseless-cast
# -Wfloat-equal

# Changes to warning based on compiler:
ifeq ($(CXX),clang++)
  WARN += -Wheader-hygiene
endif
ifeq ($(CXX),g++)
  WARN += -Wsuggest-override -Wnon-virtual-dtor -Wcast-align -Woverloaded-virtual -Wduplicated-cond -Wduplicated-branches -Wnull-dereference
#-Wsuggest-final-methods  -Wsuggest-final-types
endif

# Changes to optimisation based on build setting:
OPT=-O3
ifeq ($(Build),release)
  WARN=-w
endif
ifeq ($(Build),debug)
  UseOpenMP=no
  WARN+=-Wno-unknown-pragmas
  OPT=-O0 -g
endif

# If not using openMP, turn off 'unkown pragmas' warning.
OMP=-fopenmp
ifneq ($(UseOpenMP),yes)
  OMP=
  WARN+=-Wno-unknown-pragmas
endif

################################################################################
# Linking + Compiling:

CXXFLAGS= $(CXXSTD) $(OPT) $(OMP) $(WARN) -I$(SD)
LIBS=-lgsl -lgslcblas

# GSL location (if different from assumed default)
ifneq ($(PathForGSL),)
  CXXFLAGS+=-I$(PathForGSL)/include/
  LIBS+=-L$(PathForGSL)/lib/
endif

# Any other "extra" includes:
ifneq ($(ExtraInclude),)
  tmpInc = $(addprefix -I,$(ExtraInclude))
  CXXFLAGS+= $(tmpInc)
endif
ifneq ($(ExtraLink),)
  tmpLink = $(addprefix -L,$(ExtraLink))
  LIBS+= $(tmpLink)
endif
ifneq ($(ExtraFlags),)
  CXXFLAGS+= $(ExtraFlags)
endif

#These should be used with clang in debug mode only
MSAN = -fsanitize=memory
ASAN = -fsanitize=address -fsanitize=leak
TSAN = -fsanitize=thread
USAN = -fsanitize=undefined -fsanitize=unsigned-integer-overflow
ifeq ($(Build)+$(CXX),debug+clang++)
	CXXFLAGS += -g $(ASAN) -fno-omit-frame-pointer
endif
# echo $(Build)+$(CXX)
# MSAN_SYMBOLIZER_PATH=/usr/lib/llvm-6.0/bin/llvm-symbolizer ./diracSCAS
# ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer ./diracSCAS

#Command to compile objects and link them
COMP=$(CXX) -c -o $@ $< $(CXXFLAGS)
LINK=$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)
