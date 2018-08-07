IDIR =./src
ODIR =./obj

CXX=g++
CXXFLAGS=-I$(IDIR) -std=c++11 -Wall -fopenmp -O #-Wextra -Wpedantic
LIBS=-lgsl -lgslcblas -lm

all: h-like.x fitParametric.x parametricPotential.x atomicKernal.x hartree.x \
 wigner.x dmeXSection.x

################################################################################
## All programs depend on these header/object files:

DEPS = $(addprefix $(IDIR)/, \
 ADAMS_solveLocalBS.h ATI_atomInfo.h ElectronOrbitals.h \
 INT_quadratureIntegration.h MAT_matrixAlgebraGSL.h FPC_physicalConstants.h \
)

$(ODIR)/%.o: $(IDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

# All programs depend on these objects:
OBJ = $(addprefix $(ODIR)/, \
 ADAMS_solveLocalBS.o ElectronOrbitals.o INT_quadratureIntegration.o \
 MAT_matrixAlgebraGSL.o \
)

################################################################################
# List all other dependencies:

$(ODIR)/PRM_parametricPotentials.o: \
$(IDIR)/PRM_parametricPotentials.cpp $(IDIR)/PRM_parametricPotentials.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/WIG_369j.o: \
$(IDIR)/WIG_369j.cpp $(IDIR)/WIG_369j.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/ADAMS_solveLocalContinuum.o: \
$(IDIR)/ADAMS_solveLocalContinuum.cpp $(IDIR)/ADAMS_solveLocalContinuum.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/HF_hartree.o: \
$(IDIR)/HF_hartree.cpp $(IDIR)/HF_hartree.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/ContinuumOrbitals.o: \
$(IDIR)/ContinuumOrbitals.cpp $(IDIR)/ContinuumOrbitals.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/AKF_akFunctions.o: \
$(IDIR)/AKF_akFunctions.cpp $(IDIR)/AKF_akFunctions.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

CNTM = $(addprefix $(ODIR)/, \
 ADAMS_solveLocalContinuum.o ContinuumOrbitals.o \
)

################################################################################
# All final programs

COMP = $(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

h-like.x: $(OBJ) $(ODIR)/h-like.o
	$(COMP)

fitParametric.x: $(OBJ) $(ODIR)/fitParametric.o \
$(ODIR)/PRM_parametricPotentials.o
	$(COMP)

parametricPotential.x: $(OBJ) $(ODIR)/parametricPotential.o \
$(ODIR)/PRM_parametricPotentials.o
	$(COMP)

atomicKernal.x: $(OBJ) $(ODIR)/atomicKernal.o $(ODIR)/AKF_akFunctions.o \
$(ODIR)/WIG_369j.o $(ODIR)/PRM_parametricPotentials.o $(CNTM) \
$(ODIR)/HF_hartree.o
	$(COMP)

dmeXSection.x: $(OBJ) $(ODIR)/dmeXSection.o $(ODIR)/AKF_akFunctions.o \
$(ODIR)/WIG_369j.o $(CNTM)
	$(COMP)

hartree.x: $(OBJ) $(ODIR)/hartree.o $(ODIR)/PRM_parametricPotentials.o \
$(ODIR)/HF_hartree.o
	$(COMP)

wigner.x: $(ODIR)/wigner.o $(ODIR)/WIG_369j.o
	$(COMP)

.PHONY: clean
clean:
	rm -f *.x *~ $(ODIR)/*.o $(IDIR)/*~
