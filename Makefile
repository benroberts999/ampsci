IDIR =./src
ODIR =./obj

CXX=g++
CXXFLAGS=-I$(IDIR) -std=c++11 -Ofast -fopenmp -Wall -Wextra -Wpedantic
LIBS=-lgsl -lgslcblas -lm

all: checkObj h-like.x fitParametric.x parametricPotential.x atomicKernal.x \
hartreeFock.x wigner.x #dmeXSection.x 


################################################################################
## All programs depend on these header/object files:

DEPS = $(addprefix $(IDIR)/, \
 ADAMS_solveLocalBS.h ATI_atomInfo.h ElectronOrbitals.h \
 INT_quadratureIntegration.h MAT_matrixAlgebraGSL.h FPC_physicalConstants.h \
)

$(ODIR)/%.o: $(IDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

# All programs depend on these objects:
BASE = $(addprefix $(ODIR)/, \
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

$(ODIR)/HF_hartreeFock.o: \
$(IDIR)/HF_hartreeFock.cpp $(IDIR)/HF_hartreeFock.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/ContinuumOrbitals.o: \
$(IDIR)/ContinuumOrbitals.cpp $(IDIR)/ContinuumOrbitals.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/AKF_akFunctions.o: \
$(IDIR)/AKF_akFunctions.cpp $(IDIR)/AKF_akFunctions.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/SHM_standardHaloModel.o: \
$(IDIR)/SHM_standardHaloModel.cpp $(IDIR)/SHM_standardHaloModel.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/SBF_sphericalBesselFunctions.o: \
$(IDIR)/SBF_sphericalBesselFunctions.cpp $(IDIR)/SBF_sphericalBesselFunctions.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

CNTM = $(addprefix $(ODIR)/, \
 ADAMS_solveLocalContinuum.o ContinuumOrbitals.o \
)

################################################################################
# All final programs

COMP = $(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

h-like.x: $(BASE) $(ODIR)/h-like.o
	$(COMP)

fitParametric.x: $(BASE) $(ODIR)/fitParametric.o \
$(ODIR)/PRM_parametricPotentials.o
	$(COMP)

parametricPotential.x: $(BASE) $(ODIR)/parametricPotential.o \
$(ODIR)/PRM_parametricPotentials.o
	$(COMP)

atomicKernal.x: $(BASE) $(ODIR)/atomicKernal.o $(ODIR)/AKF_akFunctions.o \
$(ODIR)/WIG_369j.o $(ODIR)/PRM_parametricPotentials.o $(CNTM) \
$(ODIR)/HF_hartreeFock.o $(ODIR)/SBF_sphericalBesselFunctions.o
	$(COMP)

hartreeFock.x: $(BASE) $(ODIR)/hartreeFock.o \
$(ODIR)/PRM_parametricPotentials.o \
$(ODIR)/WIG_369j.o $(ODIR)/HF_hartreeFock.o
	$(COMP)

dmeXSection.x: $(BASE) $(ODIR)/dmeXSection.o $(ODIR)/AKF_akFunctions.o \
$(ODIR)/SBF_sphericalBesselFunctions.o \
$(ODIR)/WIG_369j.o $(ODIR)/PRM_parametricPotentials.o $(CNTM) \
$(ODIR)/HF_hartreeFock.o $(ODIR)/SHM_standardHaloModel.o
	$(COMP)

wigner.x: $(ODIR)/wigner.o $(ODIR)/WIG_369j.o
	$(COMP)

################################################################################

checkObj:
	@if [ ! -d $(ODIR) ]; then \
	echo '\n ERROR: Directory: '$(ODIR)' doesnt exist - please create it!\n'; \
	false; \
	else \
	echo 'OK'; \
	fi

.PHONY: clean do_the_chicken_dance checkObj
clean:
	rm -f *.x *~ $(ODIR)/*.o $(IDIR)/*~
do_the_chicken_dance:
	@echo 'Why would I do that?'
