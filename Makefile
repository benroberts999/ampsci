CXX=g++
CXXFLAGS=-std=c++11 -Wall -fopenmp -O
LIBS=-lgsl -lgslcblas -lm

# All programs depend on these:
DEPS = adamsSolveLocalBS.h atomInfo.h ElectronOrbitals.h \
       INT_quadratureIntegration.h MAT_matrixAlgebraGSL.h physicalConstants.h

%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

OBJ = adamsSolveLocalBS.o ElectronOrbitals.o INT_quadratureIntegration.o \
      MAT_matrixAlgebraGSL.o


COMP = $(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)


all: h-like.x fitParametric.x

h-like.x: h-like.o  $(OBJ)
	$(COMP)

fitParametric.x: fitParametric.o  PRM_parametricPotentials.o \
PRM_parametricPotentials.h $(OBJ)
	$(COMP)


.PHONY: clean
clean:
	rm -f *.x *~ *.o
