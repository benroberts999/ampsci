#Set directories for input files/source code (ID),
# output object files (OD), and executables (XD)
ID =./src
OD =./obj
XD =.

CXX=g++ #clang++ #

OPT=-O3
OMP=-fopenmp #-fopenmp=libiomp5 ##needed for clang++

WARN=-Wpedantic -Wall -Wextra -Wdouble-promotion -Wconversion

CXXFLAGS= -std=c++11 $(OPT) $(OMP) $(WARN) #-I$(ID)
LIBS=-lgsl -lgslcblas

#Command to compile objects and link them
COMP=$(CXX) -c -o $@ $< $(CXXFLAGS)
LINK=$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

################################################################################
#Allow exectuables to be placed in another directory:
ALLEXES = $(addprefix $(XD)/, \
 h-like fitParametric parametricPotential atomicKernal \
 hartreeFock wigner dmeXSection \
)

#Default make rule:
all: checkObj checkXdir $(ALLEXES)

################################################################################
## Dependencies:

# All programs depend on these generic common headers:
COMH = $(addprefix $(ID)/, \
 ATI_atomInfo.h FPC_physicalConstants.h FileIO_fileReadWrite.h \
 INT_quadratureIntegration.h \
)

# Rule for files that have .cpp AND a .h file
# They depend 'only' on their own header, + generic common headers
$(OD)/%.o: $(ID)/%.cpp $(ID)/%.h $(COMH)
	$(COMP)

# Rule for files that _don't_ have a .h header. (mains)
# These also depend on the common headers
$(OD)/%.o: $(ID)/%.cpp $(COMH)
	$(COMP)

# Here: List rules for any other progs that don't fit above rules?
$(OD)/dummy.o: $(ID)/dummy.cpp $(COMH) $(ID)/otherHeader.h
	$(COMP)

################################################################################
# Hust to save typing: Many programs depend on these combos:

BASE = $(addprefix $(OD)/, \
 ADAMS_solveLocalBS.o ElectronOrbitals.o \
 MAlg_matrixAlgebra.o \
)

HF = $(addprefix $(OD)/, \
 HartreeFockClass.o PRM_parametricPotentials.o WIG_369j.o \
)

CNTM = $(addprefix $(OD)/, \
 ADAMS_solveLocalContinuum.o ContinuumOrbitals.o \
)

################################################################################
# Link + build all final programs

$(XD)/h-like: $(BASE) $(OD)/h-like.o $(OD)/ChronoTimer.o
	$(LINK)

$(XD)/fitParametric: $(BASE) $(OD)/fitParametric.o \
$(OD)/PRM_parametricPotentials.o $(OD)/ChronoTimer.o
	$(LINK)

$(XD)/parametricPotential: $(BASE) $(OD)/parametricPotential.o \
$(OD)/PRM_parametricPotentials.o $(OD)/ChronoTimer.o
	$(LINK)

$(XD)/atomicKernal: $(BASE) $(CNTM) $(HF) \
$(OD)/atomicKernal.o $(OD)/AKF_akFunctions.o $(OD)/ChronoTimer.o \
$(OD)/ExponentialGrid.o
	$(LINK)

$(XD)/hartreeFock: $(BASE) $(HF) $(OD)/hartreeFock.o $(OD)/ChronoTimer.o
	$(LINK)

$(XD)/dmeXSection: $(BASE) $(CNTM) $(HF) $(OD)/dmeXSection.o \
$(OD)/AKF_akFunctions.o \
$(OD)/StandardHaloModel.o $(OD)/ChronoTimer.o $(OD)/ExponentialGrid.o
	$(LINK)

$(XD)/wigner: $(OD)/wigner.o $(OD)/WIG_369j.o
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
