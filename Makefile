IDIR =./src

CXX=g++
CXXFLAGS=-I$(IDIR) -std=c++11 -Wall -fopenmp -O #-Wextra -Wpedantic
LIBS=-lgsl -lgslcblas -lm

all: h-like.x fitParametric.x parametricPotential.x atomicKernal.x hartree.x \
 wigner.x

################################################################################
## All programs depend on these header/object files:

DEPS = $(addprefix $(IDIR)/, \
 adamsSolveLocalBS.h ATI_atomInfo.h ElectronOrbitals.h \
 INT_quadratureIntegration.h MAT_matrixAlgebraGSL.h physicalConstants.h \
)

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

# All programs depend on these objects:
OBJ = $(addprefix $(IDIR)/, \
 adamsSolveLocalBS.o ElectronOrbitals.o INT_quadratureIntegration.o \
 MAT_matrixAlgebraGSL.o \
)

################################################################################
# List all other dependencies:

$(IDIR)/PRM_parametricPotentials.o: \
$(IDIR)/PRM_parametricPotentials.cpp $(IDIR)/PRM_parametricPotentials.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(IDIR)/WIG_369j.o: \
$(IDIR)/WIG_369j.cpp $(IDIR)/WIG_369j.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(IDIR)/adamsSolveLocalContinuum.o: \
$(IDIR)/adamsSolveLocalContinuum.cpp $(IDIR)/adamsSolveLocalContinuum.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(IDIR)/HF_hartree.o: \
$(IDIR)/HF_hartree.cpp $(IDIR)/HF_hartree.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(IDIR)/ContinuumOrbitals.o: \
$(IDIR)/ContinuumOrbitals.cpp $(IDIR)/ContinuumOrbitals.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(IDIR)/akFunctions.o: \
$(IDIR)/akFunctions.cpp $(IDIR)/akFunctions.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

CNTM = $(addprefix $(IDIR)/, \
 adamsSolveLocalContinuum.o ContinuumOrbitals.o \
)

################################################################################
# All final programs

COMP = $(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

h-like.x: $(OBJ) $(IDIR)/h-like.o
	$(COMP)

fitParametric.x: $(OBJ) $(IDIR)/fitParametric.o \
$(IDIR)/PRM_parametricPotentials.o
	$(COMP)

parametricPotential.x: $(OBJ) $(IDIR)/parametricPotential.o \
$(IDIR)/PRM_parametricPotentials.o
	$(COMP)

atomicKernal.x: $(OBJ) $(IDIR)/atomicKernal.o $(IDIR)/akFunctions.o $(IDIR)/WIG_369j.o \
$(IDIR)/PRM_parametricPotentials.o $(CNTM) $(IDIR)/HF_hartree.o
	$(COMP)

# dmeXSection.x: $(OBJ) $(IDIR)/dmeXSection.o $(IDIR)/akFunctions.o $(IDIR)/WIG_369j.o
# 	$(COMP)

hartree.x: $(OBJ) $(IDIR)/hartree.o $(IDIR)/PRM_parametricPotentials.o \
$(IDIR)/HF_hartree.o
	$(COMP)

wigner.x: $(IDIR)/wigner.o $(IDIR)/WIG_369j.o
	$(COMP)

.PHONY: clean
clean:
	rm -f *.x *~ $(IDIR)/*.o $(IDIR)/*~
