IDIR =./src

CXX=g++
CXXFLAGS=-I$(IDIR) -std=c++11 -Wall -fopenmp -O
LIBS=-lgsl -lgslcblas -lm

all: h-like.x fitParametric.x parametricPotential.x #testTF.x

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


################################################################################
# All final programs

COMP = $(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

#testTF.x: $(IDIR)/TMF_thomasFermi.cpp
#	$(COMP)

h-like.x: $(OBJ) $(IDIR)/h-like.o
	$(COMP)

fitParametric.x: $(OBJ) $(IDIR)/fitParametric.o $(IDIR)/PRM_parametricPotentials.o
	$(COMP)

parametricPotential.x: $(OBJ) $(IDIR)/parametricPotential.o $(IDIR)/PRM_parametricPotentials.o
	$(COMP)

.PHONY: clean
clean:
	rm -f *.x *~ $(IDIR)/*.o $(IDIR)/*~
