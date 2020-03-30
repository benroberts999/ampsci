# Dependencies for Physics

$(BD)/AtomData.o: $(SD)/Physics/AtomData.cpp \
$(SD)/Physics/AtomData.hpp \
$(SD)/Physics/AtomData_PeriodicTable.hpp
	$(COMP)

$(BD)/Nuclear.o: $(SD)/Physics/NuclearPotentials.cpp \
$(SD)/Physics/NuclearPotentials.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Physics/AtomData.hpp \
$(SD)/Physics/NuclearData.hpp \
$(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(BD)/NuclearData.o: $(SD)/Physics/NuclearData.cpp \
$(SD)/Physics/NuclearData.hpp
	$(COMP)

$(BD)/Parametric_potentials.o: $(SD)/Physics/Parametric_potentials.cpp \
$(SD)/Physics/Parametric_potentials.hpp
	$(COMP)

$(BD)/RadiativePotential.o: $(SD)/Physics/RadiativePotential.cpp \
$(SD)/Physics/RadiativePotential.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/Maths/Interpolator.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Physics/PhysConst_constants.hpp
	$(COMP)
