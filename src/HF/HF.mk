# Dependencies for HF

$(BD)/CoulombIntegrals.o: \
$(SD)/HF/CoulombIntegrals.cpp $(SD)/HF/CoulombIntegrals.hpp \
$(SD)/Dirac/DiracSpinor.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Physics/AtomData.hpp $(SD)/Angular/Wigner_369j.hpp \
$(SD)/IO/SafeProfiler.hpp
	$(COMP)

$(BD)/ExternalField.o: \
$(SD)/HF/ExternalField.cpp $(SD)/HF/ExternalField.hpp \
$(SD)/HF/CoulombIntegrals.hpp $(SD)/HF/HartreeFockClass.hpp\
$(SD)/Dirac/DiracSpinor.hpp $(SD)/Dirac/DiracOperator.hpp \
$(SD)/Maths/Grid.hpp $(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Angular/Angular.hpp $(SD)/IO/SafeProfiler.hpp
	$(COMP)

$(BD)/HartreeFockClass.o: \
$(SD)/HF/HartreeFockClass.cpp $(SD)/HF/HartreeFockClass.hpp \
$(SD)/Adams/Adams_Greens.hpp $(SD)/Adams/DiracODE.hpp \
$(SD)/Dirac/DiracSpinor.hpp $(SD)/Dirac/Wavefunction.hpp \
$(SD)/HF/CoulombIntegrals.hpp $(SD)/Maths/Grid.hpp\
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/Physics/Parametric_potentials.hpp \
$(SD)/Angular/Wigner_369j.hpp $(SD)/IO/SafeProfiler.hpp
	$(COMP)
