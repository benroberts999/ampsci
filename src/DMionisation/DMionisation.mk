# Dependencies for DMionisation

$(OD)/AKF_akFunctions.o: $(SD)/DMionisation/AKF_akFunctions.cpp \
$(SD)/DMionisation/AKF_akFunctions.hpp \
$(SD)/Dirac/ContinuumOrbitals.hpp $(SD)/Dirac/Wavefunction.hpp \
$(SD)/IO/FileIO_fileReadWrite.hpp $(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Maths/SphericalBessel.hpp $(SD)/Physics/AtomData.hpp \
$(SD)/Physics/PhysConst_constants.hpp $(SD)/Physics/Wigner_369j.hpp
	$(COMP)

$(OD)/dmeXSection.o: $(SD)/DMionisation/dmeXSection.cpp \
$(SD)/DMionisation/AKF_akFunctions.hpp $(SD)/IO/ChronoTimer.hpp \
$(SD)/IO/FileIO_fileReadWrite.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/Physics/PhysConst_constants.hpp \
$(SD)/DMionisation/StandardHaloModel.hpp
	$(COMP)

$(OD)/Module_atomicKernal.o: $(SD)/DMionisation/Module_atomicKernal.cpp \
$(SD)/DMionisation/Module_atomicKernal.hpp \
$(SD)/DMionisation/AKF_akFunctions.hpp $(SD)/Dirac/ContinuumOrbitals.hpp \
$(SD)/Dirac/Wavefunction.hpp $(SD)/HF/HartreeFockClass.hpp \
$(SD)/IO/ChronoTimer.hpp $(SD)/IO/UserInput.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(OD)/StandardHaloModel.o: $(SD)/DMionisation/StandardHaloModel.cpp \
$(SD)/DMionisation/StandardHaloModel.hpp
	$(COMP)
