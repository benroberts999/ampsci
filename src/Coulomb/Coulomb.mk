# Dependencies for Coulomb

$(BD)/Coulomb.o: $(SD)/Coulomb/Coulomb.cpp \
$(SD)/Coulomb/Coulomb.hpp \
$(SD)/Angular/Angular_369j.hpp \
$(SD)/Angular/Angular_tables.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp
	$(COMP)

$(BD)/YkTable.o: $(SD)/Coulomb/YkTable.cpp \
$(SD)/Angular/Angular_369j.hpp \
$(SD)/Angular/Angular_tables.hpp \
$(SD)/Coulomb/Coulomb.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp
	$(COMP)
