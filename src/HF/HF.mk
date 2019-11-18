# Dependencies for HF

$(OD)/CoulombIntegrals.o: $(SD)/HF/CoulombIntegrals.cpp \
$(SD)/HF/CoulombIntegrals.hpp $(SD)/Dirac/DiracSpinor.hpp \
$(SD)/Maths/Grid.hpp $(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Physics/AtomData.hpp $(SD)/Angular/Wigner_369j.hpp
	$(COMP)

$(OD)/HartreeFockClass.o: $(SD)/HF/HartreeFockClass.cpp \
$(SD)/HF/HartreeFockClass.hpp $(SD)/Adams/Adams_Greens.hpp \
$(SD)/Adams/DiracODE.hpp $(SD)/Dirac/DiracSpinor.hpp \
$(SD)/Dirac/Wavefunction.hpp $(SD)/HF/CoulombIntegrals.hpp $(SD)/Maths/Grid.hpp\
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/Physics/Parametric_potentials.hpp \
$(SD)/Angular/Wigner_369j.hpp
	$(COMP)
