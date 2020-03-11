# Dependencies for Operators

$(BD)/DiracOperator.o: $(SD)/Operators/DiracOperator.cpp \
$(SD)/Operators/DiracOperator.hpp $(SD)/Wavefunction/DiracSpinor.cpp \
$(SD)/Wavefunction/DiracSpinor.hpp $(SD)/Angular/Angular_369j.hpp
	$(COMP)
