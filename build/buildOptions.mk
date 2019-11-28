# Options + settings for makefile for each 'Main'

# runs make in //
ifneq ($(Build),debug)
  MAKEFLAGS += -j12
endif

#Warnings:
WARN=-Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Wshadow
# -Weffc++
# -Wfloat-equal
# -Wsign-conversion

# Changes to warning based on compiler:
ifeq ($(CXX),clang++)
  WARN += -Wno-sign-conversion -Wheader-hygiene
endif
ifeq ($(CXX),g++)
  WARN += -Wsuggest-override -Wsuggest-final-types #-Wsuggest-final-methods
endif

# Changes to optimisation based on build setting:
OPT=-O3
ifeq ($(Build),release)
  WARN=-w
endif
ifeq ($(Build),debug)
#  UseOpenMP=no
#	WARN+=-Wno-unknown-pragmas
	OPT=-O0 -g
endif

OMP=-fopenmp
ifneq ($(UseOpenMP),yes)
  OMP=
	WARN+=-Wno-unknown-pragmas
endif

CXXFLAGS= $(CXXSTD) $(OPT) $(OMP) $(WARN) -I$(SD)
LIBS=-lgsl -lgslcblas

#These should be used with clang in debug mode only
MSAN = -fsanitize=memory
ASAN = -fsanitize=address -fsanitize=leak
TSAN = -fsanitize=thread
USAN = -fsanitize=undefined -fsanitize=unsigned-integer-overflow
ifeq ($(Build)+$(CXX),debug+clang++)
	CXXFLAGS += -g $(ASAN) -fno-omit-frame-pointer
endif
# echo $(Build)+$(CXX)
# MSAN_SYMBOLIZER_PATH=/usr/lib/llvm-6.0/bin/llvm-symbolizer ./hartreeFock
# ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer ./hartreeFock

#Command to compile objects and link them
COMP=$(CXX) -c -o $@ $< $(CXXFLAGS)
LINK=$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)
