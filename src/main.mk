# Dependencies for each 'Main'

$(BD)/hartreeFock.o: $(SD)/hartreeFock.cpp \
$(SD)/IO/ChronoTimer.hpp \
$(SD)/IO/FRW_fileReadWrite.hpp \
$(SD)/IO/UserInput.hpp \
$(SD)/Maths/Interpolator.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Modules/Module_runModules.hpp \
$(SD)/Physics/PhysConst_constants.hpp \
$(SD)/Wavefunction/Wavefunction.hpp
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
