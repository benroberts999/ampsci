# Dependencies for HF

$(BD)/Coulomb.o: \
$(SD)/Coulomb/Coulomb.cpp $(SD)/Coulomb/Coulomb.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/Angular/Angular_369j.hpp \
$(SD)/IO/SafeProfiler.hpp
	$(COMP)

$(BD)/YkTable.o: \
$(SD)/Coulomb/YkTable.cpp $(SD)/Coulomb/YkTable.hpp $(SD)/Coulomb/Coulomb.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Angular/Angular_tables.hpp $(SD)/Angular/Angular_369j.hpp \
$(SD)/IO/SafeProfiler.hpp
	$(COMP)
