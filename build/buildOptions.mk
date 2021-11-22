# Options + settings for makefile for each 'Main'

## will return the Operating system name
detected_OS := $(shell uname -s)
$(info )
$(info Detected operating system: $(detected_OS))

# runs make in //
ifeq ($(Build),debug)
  ParallelBuild=1
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
WARN=-Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Wshadow \
-Weffc++ -Wsign-conversion
# -Wfloat-equal

# Changes to warning based on compiler:
ifeq ($(findstring clang++,$(CXX)),clang++)
  WARN += -Wheader-hygiene -Wno-unused-function
else
ifeq ($(findstring g++,$(CXX)),g++)
  WARN += -Wsuggest-override -Wnon-virtual-dtor -Wcast-align \
  -Woverloaded-virtual -Wduplicated-cond -Wduplicated-branches \
  -Wnull-dereference -Wuseless-cast -Wformat
endif
endif
# nb: have to put g++ in 'else' block, since clanG++ contains g++
# can't use equal, since compiler might be e.g., g++-11

# Changes to optimisation based on build setting:
OPT ?= -O3

ifeq ($(Build),release)
  WARN=-w
  OPT+=-g0 -DNDEBUG -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
endif
ifeq ($(Build),debug)
  UseOpenMP=no
  WARN+=-Wno-unknown-pragmas
  OPT=-O0 -g3 -fno-omit-frame-pointer
endif
# -fno-exceptions -fno-rtti not needed
# -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF are from GSL

# If not using openMP, turn off 'unkown pragmas' warning.
OMP?=-fopenmp
ifneq ($(UseOpenMP),yes)
  OMP=
  WARN+=-Wno-unknown-pragmas
endif

################################################################################
# Linking + Compiling:

CXXFLAGS= $(CXXSTD) $(OPT) $(OMP) $(WARN) -I$(SD)
LIBS=-lgsl -lgslcblas

ifeq ($(RunProfiler),yes)
  CXXFLAGS+=-DIOPROFILER
endif

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
# MSAN_SYMBOLIZER_PATH=/usr/lib/llvm-6.0/bin/llvm-symbolizer ./ampsci
# ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer ./ampsci

#Command to compile objects and link them
COMP=$(CXX) -MMD -c -o $@ $< $(CXXFLAGS)
ifeq ($(ParallelBuild),1)
	COMP=time $(CXX) -MMD -MP -c -o $@ $< $(CXXFLAGS)
endif
LINK=$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)
