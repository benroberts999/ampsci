# Dependencies for Dirac

$(BD)/ContinuumOrbitals.o: $(SD)/Dirac/ContinuumOrbitals.cpp \
$(SD)/Dirac/ContinuumOrbitals.hpp $(SD)/Adams/DiracODE.hpp \
$(SD)/Dirac/DiracSpinor.hpp $(SD)/Dirac/Wavefunction.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Physics/AtomData.hpp $(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(BD)/DiracOperator.o: $(SD)/Dirac/DiracOperator.cpp \
$(SD)/Dirac/DiracOperator.hpp $(SD)/Dirac/DiracSpinor.cpp \
$(SD)/Dirac/DiracSpinor.hpp $(SD)/Angular/Angular_369j.hpp
	$(COMP)

$(BD)/DiracSpinor.o: $(SD)/Dirac/DiracSpinor.cpp $(SD)/Dirac/DiracSpinor.hpp \
$(SD)/Maths/Grid.hpp $(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Physics/AtomData.hpp
	$(COMP)

$(BD)/BSplineBasis.o: $(SD)/Dirac/BSplineBasis.cpp $(SD)/Dirac/BSplineBasis.hpp\
$(SD)/Maths/LinAlg_MatrixVector.hpp $(SD)/Maths/BSplines.hpp \
$(SD)/Dirac/Hamiltonian.hpp \
$(SD)/Dirac/Wavefunction.hpp $(SD)/Dirac/DiracSpinor.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/Physics/AtomData.hpp
	$(COMP)

$(BD)/Wavefunction.o: $(SD)/Dirac/Wavefunction.cpp $(SD)/Dirac/Wavefunction.hpp\
$(SD)/Adams/DiracODE.hpp $(SD)/Dirac/DiracSpinor.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/Physics/AtomData.hpp \
$(SD)/Dirac/BSplineBasis.hpp \
$(SD)/Physics/NuclearPotentials.hpp $(SD)/Physics/PhysConst_constants.hpp
	$(COMP)
