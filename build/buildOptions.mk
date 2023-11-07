# Options + settings for makefile for each 'Main'

# runs make in //
ParallelBuild ?= 1
ifeq ($(Build),debug)
  ParallelBuild=1
endif
ifneq ($(ParallelBuild),1)
  MAKEFLAGS += -j$(ParallelBuild)
  $(info Compiling in parallel with N=$(ParallelBuild) threads)
endif

# Warnings:
WARN=-Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Wshadow \
-Weffc++ -Wsign-conversion -Wno-psabi

#-Wno-psabi: #https://stackoverflow.com/questions/48149323/what-does-the-gcc-warning-project-parameter-passing-for-x-changed-in-gcc-7-1-m

# Changes to warning based on compiler:
GCC_Compilers:=g++ g++13 g++12 g++-11 g++-10 g++-9 g++-8 g++-7
CLANG_Compilers:=clang++ clang++-11 clang++-10 clang++-9 clang++-8 clang++-7 clang++-6
ifneq ($(filter $(CXX),$(GCC_Compilers)),)
  WARN += -Wsuggest-override -Wnon-virtual-dtor -Wcast-align \
  -Woverloaded-virtual -Wduplicated-cond -Wduplicated-branches \
  -Wnull-dereference -Wuseless-cast -Wformat
endif
ifneq ($(filter $(CXX),$(CLANG_Compilers)),)
  WARN += -Wheader-hygiene -Wno-unused-function
endif

#g++-7 gives some faulty warnings with c++17 (structured bindings, constexpr)
ifneq ($(filter $(CXX),g++-7),)
  WARN=-Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Weffc++ \
  -Wsign-conversion -Wno-unused-variable -Wno-shadow
endif


# Changes to optimisation based on build setting:
OPT?=-O3

ifeq ($(Build),release)
  WARN=-w -Wno-psabi
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
ifneq ($(UseOpenMP),true)
ifneq ($(UseOpenMP),yes)
  OMP=
  WARN+=-Wno-unknown-pragmas
endif
endif

################################################################################
# Allows different blas standard (e.g., openblas) by setting BLAS=openblas
BLASLIB?=blas

################################################################################
# Linking + Compiling:

CXXFLAGS= $(CXXSTD) $(OPT) $(OMP) $(WARN) -I$(SD)
LIBS=-lgsl -lgslcblas -llapack -l$(BLASLIB)

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
COMP=$(CXX) $(CARGS) -MMD -MP -c -o $@ $< $(CXXFLAGS)
ifeq ($(ParallelBuild),1)
	COMP=time $(CXX) -MMD -MP -c -o $@ $< $(CXXFLAGS)
endif
LINK=$(CXX) $(LARGS) -o $@ $^ $(CXXFLAGS) $(LIBS)
