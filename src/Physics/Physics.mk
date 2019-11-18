# Dependencies for Physics

$(OD)/AtomData.o: $(SD)/Physics/AtomData.cpp $(SD)/Physics/AtomData.hpp \
$(SD)/Physics/AtomData_PeriodicTable.hpp
	$(COMP)

$(OD)/Angular.o: $(SD)/Physics/Angular.cpp $(SD)/Physics/Angular.hpp \
$(SD)/Physics/Wigner_369j.hpp
	$(COMP)

$(OD)/Nuclear.o: $(SD)/Physics/NuclearPotentials.cpp \
$(SD)/Physics/NuclearPotentials.hpp \
$(SD)/Maths/Grid.hpp $(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Physics/AtomData.hpp $(SD)/Physics/NuclearData.hpp \
$(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(OD)/NuclearData.o: $(SD)/Physics/NuclearData.cpp \
$(SD)/Physics/NuclearData.hpp
	$(COMP)

$(OD)/Parametric_potentials.o: $(SD)/Physics/Parametric_potentials.cpp \
$(SD)/Physics/Parametric_potentials.hpp
	$(COMP)
