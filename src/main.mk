# Dependencies for each 'Main'

$(OD)/hartreeFock.o: $(SD)/hartreeFock.cpp $(SD)/Dirac/Wavefunction.hpp \
$(SD)/IO/ChronoTimer.hpp $(SD)/IO/UserInput.hpp \
$(SD)/Modules/Module_runModules.hpp
	$(COMP)

$(OD)/periodicTable.o: $(SD)/periodicTable.cpp $(SD)/Physics/Nuclear.hpp \
$(SD)/Physics/Nuclear_DataTable.hpp $(SD)/Physics/AtomInfo.hpp \
$(SD)/Physics/AtomInfo_PeriodicTable.hpp
	$(COMP)

$(OD)/wigner.o: $(SD)/wigner.cpp $(SD)/IO/FileIO_fileReadWrite.hpp \
$(SD)/Physics/Wigner_369j.hpp
	$(COMP)
