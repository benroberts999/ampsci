# Dependencies for HF

$(BD)/MixedStates.o: $(SD)/HF/MixedStates.cpp \
$(SD)/HF/MixedStates.hpp \
$(SD)/Angular/Angular_369j.hpp \
$(SD)/DiracODE/Adams_Greens.hpp \
$(SD)/DiracODE/DiracODE.hpp \
$(SD)/Coulomb/Coulomb.hpp \
$(SD)/HF/HartreeFockClass.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp
	$(COMP)

$(BD)/ExternalField.o: $(SD)/HF/ExternalField.cpp \
$(SD)/HF/ExternalField.hpp \
$(SD)/Angular/Angular_369j.hpp \
$(SD)/Angular/Angular_tables.hpp \
$(SD)/Coulomb/Coulomb.hpp \
$(SD)/DiracOperator/DiracOperator.hpp \
$(SD)/HF/HartreeFockClass.hpp \
$(SD)/HF/MixedStates.hpp \
$(SD)/IO/ChronoTimer.hpp \
$(SD)/Wavefunction/BSplineBasis.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp \
$(SD)/Wavefunction/Wavefunction.hpp
	$(COMP)

$(BD)/HartreeFockClass.o: $(SD)/HF/HartreeFockClass.cpp \
$(SD)/HF/HartreeFockClass.hpp \
$(SD)/Angular/Angular_369j.hpp \
$(SD)/Coulomb/Coulomb.hpp \
$(SD)/Coulomb/YkTable.hpp \
$(SD)/DiracODE/Adams_Greens.hpp \
$(SD)/DiracODE/DiracODE.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/MBPT/CorrelationPotential.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Physics/Parametric_potentials.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp \
$(SD)/Wavefunction/Wavefunction.hpp
	$(COMP)
