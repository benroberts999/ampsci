# Dependencies for HF

$(BD)/Coulomb.o: \
$(SD)/HF/Coulomb.cpp $(SD)/HF/Coulomb.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/Angular/Angular_369j.hpp \
$(SD)/IO/SafeProfiler.hpp
	$(COMP)

$(BD)/YkTable.o: \
$(SD)/HF/YkTable.cpp $(SD)/HF/YkTable.hpp $(SD)/HF/Coulomb.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Angular/Angular_tables.hpp $(SD)/Angular/Angular_369j.hpp \
$(SD)/IO/SafeProfiler.hpp
	$(COMP)

$(BD)/MixedStates.o: \
$(SD)/HF/MixedStates.cpp $(SD)/HF/MixedStates.hpp \
$(SD)/HF/Coulomb.hpp $(SD)/HF/HartreeFockClass.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/IO/SafeProfiler.hpp
	$(COMP)

$(BD)/ExternalField.o: \
$(SD)/HF/ExternalField.cpp $(SD)/HF/ExternalField.hpp \
$(SD)/HF/Coulomb.hpp $(SD)/HF/HartreeFockClass.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp $(SD)/DiracOperator/DiracOperator.hpp \
$(SD)/Maths/Grid.hpp $(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Angular/Angular_tables.hpp $(SD)/IO/SafeProfiler.hpp
	$(COMP)

$(BD)/HartreeFockClass.o: \
$(SD)/HF/HartreeFockClass.cpp $(SD)/HF/HartreeFockClass.hpp \
$(SD)/DiracODE/Adams_Greens.hpp $(SD)/DiracODE/DiracODE.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp $(SD)/Wavefunction/Wavefunction.hpp \
$(SD)/HF/YkTable.hpp $(SD)/HF/Coulomb.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/Physics/Parametric_potentials.hpp \
$(SD)/Angular/Angular_369j.hpp $(SD)/IO/SafeProfiler.hpp
	$(COMP)
