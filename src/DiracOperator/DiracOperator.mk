# Dependencies for Operators

$(BD)/DiracOperator.o: $(SD)/DiracOperator/DiracOperator.cpp \
$(SD)/DiracOperator/DiracOperator.hpp \
$(SD)/DiracOperator/Operators.hpp \
$(SD)/Angular/Angular_369j.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Maths/SphericalBessel.hpp \
$(SD)/Physics/NuclearPotentials.hpp \
$(SD)/Physics/PhysConst_constants.hpp
	$(COMP)
