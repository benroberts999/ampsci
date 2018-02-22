IDIR =./src

CXX=g++
CXXFLAGS=-I$(IDIR) -std=c++11 -Wall -fopenmp -O
LIBS=-lgsl -lgslcblas -lm


# All programs depend on these:
_DEPS = adamsSolveLocalBS.h atomInfo.h ElectronOrbitals.h \
       INT_quadratureIntegration.h MAT_matrixAlgebraGSL.h physicalConstants.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))


%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)


_OBJ = adamsSolveLocalBS.o ElectronOrbitals.o INT_quadratureIntegration.o \
      MAT_matrixAlgebraGSL.o
OBJ = $(patsubst %,$(IDIR)/%,$(_OBJ))


COMP = $(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)


all: h-like.x fitParametric.x

h-like.x: $(IDIR)/h-like.o  $(OBJ)
	$(COMP)

fitParametric.x: $(IDIR)/fitParametric.o  $(IDIR)/PRM_parametricPotentials.o \
$(IDIR)/PRM_parametricPotentials.h $(OBJ)
	$(COMP)


.PHONY: clean
clean:
	rm -f *.x *~ *.o
