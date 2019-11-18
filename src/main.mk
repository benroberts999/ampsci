# Dependencies for each 'Main'

$(OD)/hartreeFock.o: $(SD)/hartreeFock.cpp $(SD)/Dirac/Wavefunction.hpp \
$(SD)/IO/ChronoTimer.hpp $(SD)/IO/UserInput.hpp \
$(SD)/Modules/Module_runModules.hpp
	$(COMP)

$(OD)/periodicTable.o: $(SD)/periodicTable.cpp \
$(SD)/Physics/NuclearData.hpp $(SD)/Physics/AtomData.hpp \
$(SD)/Physics/AtomData_PeriodicTable.hpp $(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(OD)/wigner.o: $(SD)/wigner.cpp $(SD)/IO/FileIO_fileReadWrite.hpp \
$(SD)/Angular/Wigner_369j.hpp
	$(COMP)
