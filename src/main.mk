# Dependencies for each 'Main'

$(BD)/diracSCAS.o: $(SD)/diracSCAS.cpp \
$(SD)/IO/ChronoTimer.hpp \
$(SD)/IO/FRW_fileReadWrite.hpp \
$(SD)/IO/UserInput.hpp \
$(SD)/Maths/Interpolator.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Modules/Module_runModules.hpp \
$(SD)/Physics/PhysConst_constants.hpp \
$(SD)/Wavefunction/Wavefunction.hpp \
$(SD)/MBPT/CorrelationPotential.hpp
	$(COMP)

$(BD)/unitTests.o: $(SD)/unitTests.cpp \
$(SD)/qip/Check.hpp $(SD)/qip/Vector.hpp \
$(SD)/HF/HartreeFock_test.hpp $(SD)/HF/MixedStates_test.hpp \
$(SD)/HF/ExternalField_test.hpp \
$(BD)/diracSCAS.o
	$(COMP)

$(BD)/periodicTable.o: $(SD)/periodicTable.cpp \
$(SD)/Physics/NuclearData.hpp \
$(SD)/Physics/AtomData.hpp \
$(SD)/Physics/AtomData_PeriodicTable.hpp \
$(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(BD)/wigner.o: $(SD)/wigner.cpp \
$(SD)/IO/FRW_fileReadWrite.hpp \
$(SD)/Angular/Angular_369j.hpp
	$(COMP)
