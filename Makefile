#Set directories for input files/source code (ID),
# output object files (OD), and executables (XD)
ID =./src
OD =./obj
XD =.

CXX=g++
# CXX=clang++

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

# Useful for tests:
#clang++ -c -std=c++11 -O3 -Wpedantic -Wall -Wextra -Wdouble-promotion -Wconversion src/FILE -o ./junk.jnk

CXXFLAGS= -std=c++11 $(OPT) $(OMP) $(WARN)
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

# # All programs depend on these generic common headers:
# COMH = $(addprefix $(ID)/, \
#  DiracSpinor.hpp DiracOperator.hpp Wavefunction.hpp Operators.hpp \
#  Nuclear.hpp Nuclear_DataTable.hpp \
#  AtomInfo.hpp PhysConst_constants.hpp FileIO_fileReadWrite.hpp \
#  NumCalc_quadIntegrate.hpp Grid.hpp Matrix_linalg.hpp ChronoTimer.hpp \
#  Wigner_369j.hpp \
# )

# # Rule for files that have .cpp AND a .h file
# # They depend 'only' on their own header, + generic common headers
# $(OD)/%.o: $(ID)/%.cpp $(ID)/%.hpp $(COMH)
# 	$(COMP)
#
# # Rule for files that _don't_ have a .h header. (mains)
# # These also depend on the common headers
# $(OD)/%.o: $(ID)/%.cpp $(COMH)
# 	$(COMP)
#
# 	# Here: List rules for any other progs that don't fit above rules?
# $(OD)/dummy.o: $(ID)/dummy.cpp $(COMH) $(ID)/otherHeader.hpp
# 	$(COMP)

$(OD)/ADAMS_bound.o: $(ID)/ADAMS_bound.cpp $(ID)/ADAMS_bound.hpp \
$(ID)/DiracSpinor.hpp $(ID)/Grid.hpp $(ID)/Matrix_linalg.hpp \
$(ID)/NumCalc_quadIntegrate.hpp
	$(COMP)

$(OD)/ADAMS_continuum.o: $(ID)/ADAMS_continuum.cpp $(ID)/ADAMS_bound.hpp \
$(ID)/ADAMS_continuum.hpp $(ID)/DiracSpinor.hpp $(ID)/Grid.hpp \
$(ID)/NumCalc_quadIntegrate.hpp
	$(COMP)

$(OD)/AKF_akFunctions.o: $(ID)/DMionisation/AKF_akFunctions.cpp \
$(ID)/DMionisation/AKF_akFunctions.hpp \
$(ID)/AtomInfo.hpp $(ID)/ContinuumOrbitals.hpp $(ID)/Wavefunction.hpp \
$(ID)/FileIO_fileReadWrite.hpp $(ID)/NumCalc_quadIntegrate.hpp \
$(ID)/PhysConst_constants.hpp $(ID)/DMionisation/SBF_sphericalBessel.hpp \
$(ID)/Wigner_369j.hpp
	$(COMP)

$(OD)/Module_runModules.o: $(ID)/Module_runModules.cpp \
$(ID)/Module_runModules.hpp $(ID)/DiracOperator.hpp $(ID)/HartreeFockClass.hpp \
$(ID)/DMionisation/Module_atomicKernal.hpp \
$(ID)/Operators.hpp $(ID)/UserInput.hpp $(ID)/Wavefunction.hpp
	$(COMP)

$(OD)/Module_atomicKernal.o: $(ID)/DMionisation/Module_atomicKernal.cpp \
$(ID)/DMionisation/Module_atomicKernal.hpp $(ID)/DMionisation/AKF_akFunctions.hpp \
$(ID)/AtomInfo.hpp $(ID)/ChronoTimer.hpp $(ID)/ContinuumOrbitals.hpp \
$(ID)/Grid.hpp $(ID)/PhysConst_constants.hpp \
$(ID)/Wavefunction.hpp
	$(COMP)

$(OD)/ContinuumOrbitals.o: $(ID)/ContinuumOrbitals.cpp \
$(ID)/ContinuumOrbitals.hpp $(ID)/ADAMS_bound.hpp $(ID)/ADAMS_continuum.hpp \
$(ID)/AtomInfo.hpp $(ID)/Wavefunction.hpp $(ID)/Grid.hpp \
$(ID)/PhysConst_constants.hpp
	$(COMP)

$(OD)/CoulombIntegrals.o: $(ID)/CoulombIntegrals.cpp $(ID)/CoulombIntegrals.hpp\
$(ID)/DiracSpinor.hpp $(ID)/NumCalc_quadIntegrate.hpp
	$(COMP)

$(OD)/dmeXSection.o: $(ID)/DMionisation/dmeXSection.cpp \
$(ID)/DMionisation/AKF_akFunctions.hpp $(ID)/ChronoTimer.hpp \
$(ID)/FileIO_fileReadWrite.hpp $(ID)/Grid.hpp $(ID)/NumCalc_quadIntegrate.hpp \
$(ID)/PhysConst_constants.hpp $(ID)/DMionisation/StandardHaloModel.hpp
	$(COMP)

$(OD)/fitParametric.o: $(ID)/fitParametric.cpp \
$(ID)/AtomInfo.hpp $(ID)/ChronoTimer.hpp $(ID)/FileIO_fileReadWrite.hpp \
$(ID)/HartreeFockClass.hpp $(ID)/NumCalc_quadIntegrate.hpp \
$(ID)/Parametric_potentials.hpp $(ID)/PhysConst_constants.hpp \
$(ID)/Wavefunction.hpp
	$(COMP)

$(OD)/hartreeFock.o: $(ID)/hartreeFock.cpp $(ID)/AtomInfo.hpp \
$(ID)/ChronoTimer.hpp $(ID)/DiracOperator.hpp $(ID)/UserInput.hpp \
$(ID)/HartreeFockClass.hpp $(ID)/Nuclear.hpp \
$(ID)/Operators.hpp \
$(ID)/PhysConst_constants.hpp $(ID)/Wavefunction.hpp
	$(COMP)

$(OD)/HartreeFockClass.o: $(ID)/HartreeFockClass.cpp $(ID)/HartreeFockClass.hpp\
$(ID)/AtomInfo.hpp $(ID)/CoulombIntegrals.hpp \
$(ID)/DiracSpinor.hpp $(ID)/Grid.hpp $(ID)/NumCalc_quadIntegrate.hpp \
$(ID)/Parametric_potentials.hpp $(ID)/Wavefunction.hpp $(ID)/Wigner_369j.hpp
	$(COMP)

$(OD)/Parametric_potentials.o: $(ID)/Parametric_potentials.cpp\
$(ID)/Parametric_potentials.hpp $(ID)/AtomInfo.hpp $(ID)/PhysConst_constants.hpp
	$(COMP)

$(OD)/StandardHaloModel.o: $(ID)/DMionisation/StandardHaloModel.cpp \
$(ID)/DMionisation/StandardHaloModel.hpp $(ID)/Grid.hpp
	$(COMP)

$(OD)/UserInput.o: $(ID)/UserInput.cpp $(ID)/UserInput.hpp \
$(ID)/FileIO_fileReadWrite.hpp
	$(COMP)

$(OD)/Wavefunction.o: $(ID)/Wavefunction.cpp $(ID)/Wavefunction.hpp \
$(ID)/ADAMS_bound.hpp $(ID)/AtomInfo.hpp $(ID)/DiracSpinor.hpp $(ID)/Grid.hpp \
$(ID)/Nuclear.hpp $(ID)/PhysConst_constants.hpp
	$(COMP)

$(OD)/wigner.o: $(ID)/wigner.cpp $(ID)/FileIO_fileReadWrite.hpp \
$(ID)/Wigner_369j.hpp
	$(COMP)

$(OD)/nuclearData.o: $(ID)/nuclearData.cpp $(ID)/Nuclear.hpp \
$(ID)/Nuclear_DataTable.hpp $(ID)/AtomInfo.hpp
	$(COMP)

################################################################################
# Hust to save typing: Many programs depend on these combos:

BASE = $(addprefix $(OD)/, \
 ADAMS_bound.o Wavefunction.o \
)

HF = $(addprefix $(OD)/, \
 HartreeFockClass.o CoulombIntegrals.o Parametric_potentials.o \
)

CNTM = $(addprefix $(OD)/, \
 ADAMS_continuum.o ContinuumOrbitals.o \
)

################################################################################
# Link + build all final programs

$(XD)/h-like: $(BASE) $(OD)/h-like.o
	$(LINK)

$(XD)/fitParametric: $(BASE) $(HF) $(OD)/fitParametric.o \
$(OD)/Parametric_potentials.o
	$(LINK)

$(XD)/hartreeFock: $(BASE) $(HF) $(CNTM) $(OD)/hartreeFock.o $(OD)/UserInput.o \
$(OD)/AKF_akFunctions.o $(OD)/Module_runModules.o $(OD)/Module_atomicKernal.o
	$(LINK)

$(XD)/dmeXSection: $(BASE) $(CNTM) $(HF) $(OD)/dmeXSection.o \
$(OD)/AKF_akFunctions.o $(OD)/StandardHaloModel.o
	$(LINK)

$(XD)/wigner: $(OD)/wigner.o
	$(LINK)

$(XD)/nuclearData: $(OD)/nuclearData.o
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
