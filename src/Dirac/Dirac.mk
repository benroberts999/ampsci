# Dependencies for Dirac

$(OD)/ContinuumOrbitals.o: $(SD)/Dirac/ContinuumOrbitals.cpp \
$(SD)/Dirac/ContinuumOrbitals.hpp $(SD)/Adams/DiracODE.hpp \
$(SD)/Dirac/DiracSpinor.hpp $(SD)/Dirac/Wavefunction.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Physics/AtomInfo.hpp $(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(OD)/DiracSpinor.o: $(SD)/Dirac/DiracSpinor.cpp $(SD)/Dirac/DiracSpinor.hpp \
$(SD)/Maths/Grid.hpp $(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Physics/AtomInfo.hpp
	$(COMP)

$(OD)/Wavefunction.o: $(SD)/Dirac/Wavefunction.cpp $(SD)/Dirac/Wavefunction.hpp\
$(SD)/Adams/DiracODE.hpp $(SD)/Dirac/DiracSpinor.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/Physics/AtomInfo.hpp \
$(SD)/Physics/Nuclear.hpp $(SD)/Physics/PhysConst_constants.hpp
	$(COMP)
