# Dependencies for Operators

$(BD)/DiracOperator.o: $(SD)/DiracOperator/DiracOperator.cpp \
$(SD)/DiracOperator/DiracOperator.hpp $(SD)/Wavefunction/DiracSpinor.cpp \
$(SD)/Wavefunction/DiracSpinor.hpp $(SD)/Angular/Angular_369j.hpp
	$(COMP)
