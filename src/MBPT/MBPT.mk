# Dependencies for MBPT

$(BD)/CorrelationPotential.o: $(SD)/MBPT/CorrelationPotential.cpp \
$(SD)/MBPT/CorrelationPotential.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp $(SD)/Wavefunction/Wavefunction.hpp \
$(SD)/Coulomb/Coulomb.hpp $(SD)/Coulomb/YkTable.hpp \
$(SD)/Angular/Angular_tables.hpp
	$(COMP)
