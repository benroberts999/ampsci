#Set directories for input files/source code (ID),
# output object files (OD), and executables (XD)
ID =./src
OD =./obj
XD =.

CXX=g++
#CXX=clang++

OPT=-O3
OMP=-fopenmp

WARN=-Wpedantic -Wall -Wextra -Wdouble-promotion -Wconversion
# -Weffc++
# -Wshadow
# -Wfloat-equal
# -Wsign-conversion
# -fmax-errors=n # useful when changing a lot

ifeq ($(CXX),clang++)
  WARN += -Wno-sign-conversion -Wheader-hygiene
endif
ifeq ($(CXX),g++)
  WARN += -Wlogical-op
endif

CXXFLAGS= -std=c++14 $(OPT) $(OMP) $(WARN)
LIBS=-lgsl -lgslcblas

MSAN = -fsanitize=memory
ASAN = -fsanitize=address
TSAN = -fsanitize=thread
USAN = -fsanitize=undefined -fsanitize=unsigned-integer-overflow
#CXXFLAGS += -g $(TSAN) -fno-omit-frame-pointer

#Command to compile objects and link them
COMP=$(CXX) -c -o $@ $< $(CXXFLAGS)
LINK=$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

################################################################################
#Allow exectuables to be placed in another directory:
ALLEXES = $(addprefix $(XD)/, \
 fitParametric hartreeFock wigner dmeXSection nuclearData \
)

DEFAULTEXES = $(addprefix $(XD)/, \
 hartreeFock wigner nuclearData dmeXSection \
)

#Default make rule:
all: checkObj checkXdir $(DEFAULTEXES)

################################################################################
## Dependencies:

$(OD)/Adams_bound.o: $(ID)/Adams/Adams_bound.cpp $(ID)/Adams/Adams_bound.hpp \
$(ID)/DiracSpinor.hpp $(ID)/Grid.hpp $(ID)/Matrix_linalg.hpp \
$(ID)/NumCalc_quadIntegrate.hpp
	$(COMP)

$(OD)/AtomInfo.o: $(ID)/Physics/AtomInfo.cpp $(ID)/Physics/AtomInfo.hpp
	$(COMP)

$(OD)/Adams_continuum.o: $(ID)/Adams/Adams_continuum.cpp \
$(ID)/Adams/Adams_bound.hpp $(ID)/Adams/Adams_continuum.hpp\
$(ID)/DiracSpinor.hpp $(ID)/Grid.hpp $(ID)/NumCalc_quadIntegrate.hpp
	$(COMP)

$(OD)/AKF_akFunctions.o: $(ID)/DMionisation/AKF_akFunctions.cpp \
$(ID)/DMionisation/AKF_akFunctions.hpp \
$(ID)/Physics/AtomInfo.hpp $(ID)/ContinuumOrbitals.hpp $(ID)/Wavefunction.hpp \
$(ID)/FileIO_fileReadWrite.hpp $(ID)/NumCalc_quadIntegrate.hpp \
$(ID)/Physics/PhysConst_constants.hpp \
$(ID)/DMionisation/SBF_sphericalBessel.hpp $(ID)/Physics/Wigner_369j.hpp
	$(COMP)

$(OD)/Module_runModules.o: $(ID)/Module_runModules.cpp \
$(ID)/Module_runModules.hpp $(ID)/DiracOperator.hpp $(ID)/HartreeFockClass.hpp \
$(ID)/DMionisation/Module_atomicKernal.hpp \
$(ID)/Operators.hpp $(ID)/UserInput.hpp $(ID)/Wavefunction.hpp
	$(COMP)

$(OD)/Module_atomicKernal.o: $(ID)/DMionisation/Module_atomicKernal.cpp \
$(ID)/DMionisation/Module_atomicKernal.hpp \
$(ID)/DMionisation/AKF_akFunctions.hpp \
$(ID)/Physics/AtomInfo.hpp $(ID)/ChronoTimer.hpp $(ID)/ContinuumOrbitals.hpp \
$(ID)/Grid.hpp $(ID)/Physics/PhysConst_constants.hpp \
$(ID)/Wavefunction.hpp
	$(COMP)

$(OD)/Module_matrixElements.o: $(ID)/Module_matrixElements.cpp \
$(ID)/Module_matrixElements.hpp $(ID)/Physics/PhysConst_constants.hpp \
$(ID)/Physics/Nuclear.hpp $(ID)/Operators.hpp $(ID)/UserInput.hpp  \
$(ID)/HartreeFockClass.hpp $(ID)/Wavefunction.hpp
	$(COMP)

$(OD)/ContinuumOrbitals.o: $(ID)/ContinuumOrbitals.cpp \
$(ID)/ContinuumOrbitals.hpp $(ID)/Adams/Adams_bound.hpp \
$(ID)/Adams/Adams_continuum.hpp $(ID)/Physics/AtomInfo.hpp \
$(ID)/Wavefunction.hpp $(ID)/Grid.hpp $(ID)/Physics/PhysConst_constants.hpp
	$(COMP)

$(OD)/CoulombIntegrals.o: $(ID)/CoulombIntegrals.cpp $(ID)/CoulombIntegrals.hpp\
$(ID)/DiracSpinor.hpp $(ID)/NumCalc_quadIntegrate.hpp
	$(COMP)

$(OD)/dmeXSection.o: $(ID)/DMionisation/dmeXSection.cpp \
$(ID)/DMionisation/AKF_akFunctions.hpp $(ID)/ChronoTimer.hpp \
$(ID)/FileIO_fileReadWrite.hpp $(ID)/Grid.hpp $(ID)/NumCalc_quadIntegrate.hpp \
$(ID)/Physics/PhysConst_constants.hpp $(ID)/DMionisation/StandardHaloModel.hpp
	$(COMP)

$(OD)/fitParametric.o: $(ID)/fitParametric.cpp \
$(ID)/Physics/AtomInfo.hpp $(ID)/ChronoTimer.hpp $(ID)/FileIO_fileReadWrite.hpp\
$(ID)/HartreeFockClass.hpp $(ID)/NumCalc_quadIntegrate.hpp \
$(ID)/Parametric_potentials.hpp $(ID)/Physics/PhysConst_constants.hpp \
$(ID)/Wavefunction.hpp
	$(COMP)

$(OD)/hartreeFock.o: $(ID)/hartreeFock.cpp $(ID)/Physics/AtomInfo.hpp \
$(ID)/ChronoTimer.hpp $(ID)/DiracOperator.hpp $(ID)/UserInput.hpp \
$(ID)/Physics/Nuclear.hpp $(ID)/Operators.hpp $(ID)/Wavefunction.hpp
	$(COMP)

$(OD)/HartreeFockClass.o: $(ID)/HartreeFockClass.cpp $(ID)/HartreeFockClass.hpp\
$(ID)/Physics/AtomInfo.hpp $(ID)/CoulombIntegrals.hpp $(ID)/DiracSpinor.hpp \
$(ID)/Grid.hpp $(ID)/NumCalc_quadIntegrate.hpp $(ID)/Wavefunction.hpp \
$(ID)/Parametric_potentials.hpp $(ID)/Physics/Wigner_369j.hpp
	$(COMP)

$(OD)/Parametric_potentials.o: $(ID)/Parametric_potentials.cpp \
$(ID)/Parametric_potentials.hpp
	$(COMP)

$(OD)/StandardHaloModel.o: $(ID)/DMionisation/StandardHaloModel.cpp \
$(ID)/DMionisation/StandardHaloModel.hpp $(ID)/Grid.hpp
	$(COMP)

$(OD)/UserInput.o: $(ID)/UserInput.cpp $(ID)/UserInput.hpp \
$(ID)/FileIO_fileReadWrite.hpp
	$(COMP)

$(OD)/Wavefunction.o: $(ID)/Wavefunction.cpp $(ID)/Wavefunction.hpp \
$(ID)/Adams/Adams_bound.hpp $(ID)/Physics/AtomInfo.hpp $(ID)/DiracSpinor.hpp \
$(ID)/Grid.hpp $(ID)/Physics/Nuclear.hpp $(ID)/Physics/PhysConst_constants.hpp
	$(COMP)

$(OD)/wigner.o: $(ID)/wigner.cpp $(ID)/FileIO_fileReadWrite.hpp \
$(ID)/Physics/Wigner_369j.hpp
	$(COMP)

$(OD)/nuclearData.o: $(ID)/nuclearData.cpp $(ID)/Physics/Nuclear.hpp \
$(ID)/Physics/Nuclear_DataTable.hpp $(ID)/Physics/AtomInfo.hpp
	$(COMP)

################################################################################
# Hust to save typing: Many programs depend on these combos:

BASE = $(addprefix $(OD)/, \
 Adams_bound.o Wavefunction.o AtomInfo.o\
)

HF = $(addprefix $(OD)/, \
 HartreeFockClass.o CoulombIntegrals.o Parametric_potentials.o \
)

CNTM = $(addprefix $(OD)/, \
 Adams_continuum.o ContinuumOrbitals.o \
)

MODS = $(addprefix $(OD)/, \
 Module_runModules.o Module_atomicKernal.o AKF_akFunctions.o \
 Module_matrixElements.o \
)

################################################################################
# Link + build all final programs

$(XD)/fitParametric: $(BASE) $(HF) $(OD)/fitParametric.o \
$(OD)/Parametric_potentials.o
	$(LINK)

$(XD)/hartreeFock: $(BASE) $(HF) $(CNTM) $(OD)/hartreeFock.o $(OD)/UserInput.o \
$(MODS)
	$(LINK)

$(XD)/dmeXSection: $(BASE) $(CNTM) $(HF) $(OD)/dmeXSection.o \
$(OD)/AKF_akFunctions.o $(OD)/StandardHaloModel.o
	$(LINK)

$(XD)/wigner: $(OD)/wigner.o
	$(LINK)

$(XD)/nuclearData: $(OD)/nuclearData.o $(OD)/AtomInfo.o
	$(LINK)

################################################################################

checkObj:
	@if [ ! -d $(OD) ]; then \
	  echo '\n ERROR: Directory: '$(OD)' doesnt exist - please create it!\n'; \
	  false; \
	else \
	  echo 'OK'; \
	fi

checkXdir:
	@if [ ! -d $(XD) ]; then \
		echo '\n ERROR: Directory: '$(XD)' doesnt exist - please create it!\n'; \
		false; \
	fi

.PHONY: clean do_the_chicken_dance checkObj checkXdir
clean:
	rm -f $(ALLEXES) $(OD)/*.o
do_the_chicken_dance:
	@echo 'Why would I do that?'
