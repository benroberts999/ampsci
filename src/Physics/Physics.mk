# Dependencies for Physics

$(OD)/AtomInfo.o: $(SD)/Physics/AtomInfo.cpp $(SD)/Physics/AtomInfo.hpp \
$(SD)/Physics/AtomInfo_PeriodicTable.hpp
	$(COMP)

$(OD)/Angular.o: $(SD)/Physics/Angular.cpp $(SD)/Physics/Angular.hpp \
$(SD)/Physics/Wigner_369j.hpp
	$(COMP)

$(OD)/Nuclear.o: $(SD)/Physics/Nuclear.cpp $(SD)/Physics/Nuclear.hpp \
$(SD)/Maths/Grid.hpp $(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Physics/AtomInfo.hpp $(SD)/Physics/Nuclear_DataTable.hpp \
$(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(OD)/Parametric_potentials.o: $(SD)/Physics/Parametric_potentials.cpp \
$(SD)/Physics/Parametric_potentials.hpp
	$(COMP)
