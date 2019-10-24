# Which compiler: (g++, clang++) [no spaces]
CXX=g++
#CXX=clang++

# Use OpenMP (parallelisation) yes/no:
UseOpenMP=yes

#release, dev, debug (changes warnings + optimisation level)
#Build=release
Build=dev

#optional: set directory for executables (by default: current directory)
XD=.
################################################################################
################################################################################
#Set directories for input files/source code (ID),
# output object files (OD), and executables (XD)
ID=./src
OD=./obj

# runs make in //
ifneq ($(Build),debug)
  MAKEFLAGS += -j12
endif

WARN=-Wall -Wpedantic -Wextra -Wdouble-promotion -Wconversion -Wshadow
# -Weffc++
# -Wfloat-equal
# -Wsign-conversion

ifeq ($(CXX),clang++)
  WARN += -Wno-sign-conversion -Wheader-hygiene
endif
ifeq ($(CXX),g++)
  WARN += -Wsuggest-override -Wsuggest-final-types -Wsuggest-final-methods
endif

OPT=-O3
ifeq ($(Build),release)
  WARN=-w
endif
ifeq ($(Build),debug)
  UseOpenMP=no
	WARN+=-Wno-unknown-pragmas
	OPT=-O0 -g
endif

OMP=-fopenmp
ifneq ($(UseOpenMP),yes)
  OMP=
	WARN+=-Wno-unknown-pragmas
endif

CXXFLAGS= -std=c++17 $(OPT) $(OMP) $(WARN) -I$(ID)
LIBS=-lgsl -lgslcblas

#These should be used with clang in debug mode only
MSAN = -fsanitize=memory
ASAN = -fsanitize=address
TSAN = -fsanitize=thread
USAN = -fsanitize=undefined -fsanitize=unsigned-integer-overflow
#CXXFLAGS += -g $(MSAN) -fno-omit-frame-pointer
# MSAN_SYMBOLIZER_PATH=/usr/lib/llvm-6.0/bin/llvm-symbolizer ./hartreeFock

#Command to compile objects and link them
COMP=$(CXX) -c -o $@ $< $(CXXFLAGS)
LINK=$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

################################################################################
#Allow exectuables to be placed in another directory:
ALLEXES = $(addprefix $(XD)/, \
 hartreeFock wigner dmeXSection periodicTable \
)

DEFAULTEXES = $(addprefix $(XD)/, \
 hartreeFock wigner periodicTable \
)
#dmeXSection

#Default make rule:
all: checkObj checkXdir $(DEFAULTEXES) doneMessage

################################################################################
## Dependencies: ... this is getting dumb... CMAKE?

$(OD)/Adams_bound.o: $(ID)/Adams/Adams_bound.cpp $(ID)/Adams/Adams_bound.hpp \
$(ID)/Adams/DiracODE.hpp $(ID)/Adams/Adams_coefs.hpp \
$(ID)/Maths/Matrix_linalg.hpp $(ID)/Maths/NumCalc_quadIntegrate.hpp
	$(COMP)

$(OD)/Adams_continuum.o: $(ID)/Adams/Adams_continuum.cpp \
$(ID)/Adams/Adams_continuum.hpp $(ID)/Adams/DiracODE.hpp \
$(ID)/Adams/Adams_bound.hpp
	$(COMP)

$(OD)/Adams_Greens.o: $(ID)/Adams/Adams_Greens.cpp $(ID)/Adams/DiracODE.hpp \
$(ID)/Adams/Adams_Greens.hpp $(ID)/Adams/Adams_bound.hpp
	$(COMP)

$(OD)/AtomInfo.o: $(ID)/Physics/AtomInfo.cpp $(ID)/Physics/AtomInfo.hpp \
$(ID)/Physics/AtomInfo_PeriodicTable.hpp
	$(COMP)

$(OD)/AKF_akFunctions.o: $(ID)/DMionisation/AKF_akFunctions.cpp \
$(ID)/DMionisation/AKF_akFunctions.hpp $(ID)/Physics/Wigner_369j.hpp \
$(ID)/IO/FileIO_fileReadWrite.hpp $(ID)/Maths/NumCalc_quadIntegrate.hpp \
$(ID)/Physics/PhysConst_constants.hpp $(ID)/Maths/SphericalBessel.hpp
	$(COMP)

$(OD)/ContinuumOrbitals.o: $(ID)/Dirac/ContinuumOrbitals.cpp \
$(ID)/Dirac/ContinuumOrbitals.hpp $(ID)/Adams/DiracODE.hpp \
$(ID)/Physics/PhysConst_constants.hpp
	$(COMP)

$(OD)/CoulombIntegrals.o: $(ID)/HF/CoulombIntegrals.cpp \
$(ID)/HF/CoulombIntegrals.hpp $(ID)/Maths/NumCalc_quadIntegrate.hpp
	$(COMP)

$(OD)/DiracSpinor.o: $(ID)/Dirac/DiracSpinor.cpp $(ID)/Dirac/DiracSpinor.hpp \
$(ID)/Maths/NumCalc_quadIntegrate.hpp
	$(COMP)

$(OD)/dmeXSection.o: $(ID)/DMionisation/dmeXSection.cpp \
$(ID)/IO/FileIO_fileReadWrite.hpp $(ID)/Maths/NumCalc_quadIntegrate.hpp \
$(ID)/Physics/PhysConst_constants.hpp
	$(COMP)

$(OD)/Grid.o: $(ID)/Maths/Grid.cpp $(ID)/Maths/Grid.hpp
	$(COMP)

$(OD)/hartreeFock.o: $(ID)/hartreeFock.cpp $(ID)/Physics/AtomInfo.hpp \
$(ID)/Physics/Nuclear.hpp $(ID)/Maths/Grid.hpp
	$(COMP)

$(OD)/HartreeFockClass.o: $(ID)/HF/HartreeFockClass.cpp \
$(ID)/HF/HartreeFockClass.hpp $(ID)/Physics/Wigner_369j.hpp \
$(ID)/Dirac/Operators.hpp $(ID)/Dirac/DiracOperator.hpp \
$(ID)/Maths/NumCalc_quadIntegrate.hpp $(ID)/Adams/DiracODE.hpp
	$(COMP)

$(OD)/Module_runModules.o: $(ID)/Modules/Module_runModules.cpp \
$(ID)/Modules/Module_runModules.hpp $(ID)/Dirac/DiracOperator.hpp \
$(ID)/Dirac/Operators.hpp
	$(COMP)

$(OD)/Module_atomicKernal.o: $(ID)/DMionisation/Module_atomicKernal.cpp \
$(ID)/DMionisation/Module_atomicKernal.hpp \
$(ID)/DMionisation/AKF_akFunctions.hpp $(ID)/Physics/PhysConst_constants.hpp
	$(COMP)

$(OD)/Module_fitParametric.o: $(ID)/Modules/Module_fitParametric.cpp \
$(ID)/Modules/Module_fitParametric.hpp $(ID)/Modules/Module_fitParametric.hpp \
$(ID)/Physics/Nuclear.hpp
	$(COMP)

$(OD)/Module_matrixElements.o: $(ID)/Modules/Module_matrixElements.cpp \
$(ID)/Modules/Module_matrixElements.hpp $(ID)/Physics/PhysConst_constants.hpp \
$(ID)/Physics/Nuclear.hpp $(ID)/Dirac/Operators.hpp
	$(COMP)

$(OD)/Nuclear.o: $(ID)/Physics/Nuclear.cpp $(ID)/Physics/Nuclear.hpp \
$(ID)/Physics/Nuclear_DataTable.hpp
	$(COMP)

$(OD)/Parametric_potentials.o: $(ID)/Physics/Parametric_potentials.cpp \
$(ID)/Physics/Parametric_potentials.hpp
	$(COMP)

$(OD)/periodicTable.o: $(ID)/periodicTable.cpp $(ID)/Physics/Nuclear.hpp \
$(ID)/Physics/Nuclear_DataTable.hpp $(ID)/Physics/AtomInfo.hpp \
$(ID)/Physics/AtomInfo_PeriodicTable.hpp
	$(COMP)

$(OD)/StandardHaloModel.o: $(ID)/DMionisation/StandardHaloModel.cpp \
$(ID)/DMionisation/StandardHaloModel.hpp
	$(COMP)

$(OD)/UserInput.o: $(ID)/IO/UserInput.cpp $(ID)/IO/UserInput.hpp \
$(ID)/IO/FileIO_fileReadWrite.hpp
	$(COMP)

$(OD)/Wavefunction.o: $(ID)/Dirac/Wavefunction.cpp $(ID)/Dirac/Wavefunction.hpp\
$(ID)/Adams/DiracODE.hpp \
$(ID)/Physics/Nuclear.hpp $(ID)/Physics/PhysConst_constants.hpp
	$(COMP)

$(OD)/wigner.o: $(ID)/wigner.cpp $(ID)/IO/FileIO_fileReadWrite.hpp \
$(ID)/Physics/Wigner_369j.hpp
	$(COMP)

################################################################################
# Just to save typing: Many programs depend on these combos:

BASE = $(addprefix $(OD)/, \
 Adams_bound.o Wavefunction.o DiracSpinor.o AtomInfo.o Nuclear.o Grid.o \
)

HF = $(addprefix $(OD)/, \
 HartreeFockClass.o CoulombIntegrals.o Parametric_potentials.o \
 Adams_Greens.o \
)

CNTM = $(addprefix $(OD)/, \
 Adams_continuum.o ContinuumOrbitals.o \
)

MODS = $(addprefix $(OD)/, \
 Module_runModules.o Module_atomicKernal.o AKF_akFunctions.o \
 Module_matrixElements.o Module_fitParametric.o \
)

################################################################################
# Link + build all final programs

$(XD)/hartreeFock: $(BASE) $(HF) $(CNTM) $(OD)/hartreeFock.o \
$(OD)/UserInput.o $(MODS)
	$(LINK)

$(XD)/dmeXSection: $(BASE) $(CNTM) $(HF) $(OD)/dmeXSection.o \
$(OD)/AKF_akFunctions.o $(OD)/StandardHaloModel.o
	$(LINK)

$(XD)/wigner: $(OD)/wigner.o
	$(LINK)

$(XD)/periodicTable: $(OD)/periodicTable.o $(OD)/AtomInfo.o $(OD)/Nuclear.o
	$(LINK)

################################################################################

checkObj:
	@if [ ! -d $(OD) ]; then \
	  echo '\n ERROR: Directory: '$(OD)' doesnt exist - please create it!\n'; \
	  false; \
	fi

checkXdir:
	@if [ ! -d $(XD) ]; then \
		echo '\n ERROR: Directory: '$(XD)' doesnt exist - please create it!\n'; \
		false; \
	fi

doneMessage:
		@echo 'done'

.PHONY: clean do_the_chicken_dance checkObj checkXdir
clean:
	rm -f $(ALLEXES) $(OD)/*.o
do_the_chicken_dance:
	@echo 'Why would I do that?'
