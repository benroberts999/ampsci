# Dependencies for Modules

MODULELIST = $(addprefix $(BD)/, \
 Module_runModules.o Module_atomicKernal.o AKF_akFunctions.o \
 Module_matrixElements.o Module_fitParametric.o Module_pnc.o Module_tests.o \
)

$(BD)/Module_runModules.o: $(SD)/Modules/Module_runModules.cpp \
$(SD)/Modules/Module_runModules.hpp \
$(SD)/DMionisation/Module_atomicKernal.hpp \
$(SD)/IO/UserInput.hpp \
$(SD)/Modules/Module_fitParametric.hpp \
$(SD)/Modules/Module_matrixElements.hpp \
$(SD)/Modules/Module_pnc.hpp \
$(SD)/Modules/Module_tests.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp \
$(SD)/Wavefunction/Wavefunction.hpp
	$(COMP)

$(BD)/Module_tests.o: $(SD)/Modules/Module_tests.cpp \
$(SD)/Modules/Module_tests.hpp \
$(SD)/DiracOperator/Operators.hpp \
$(SD)/HF/HartreeFock.hpp \
$(SD)/IO/UserInput.hpp \
$(SD)/Physics/PhysConst_constants.hpp \
$(SD)/Wavefunction/Hamiltonian.hpp \
$(SD)/Wavefunction/Wavefunction.hpp
	$(COMP)

$(BD)/Module_fitParametric.o: $(SD)/Modules/Module_fitParametric.cpp \
$(SD)/Modules/Module_fitParametric.hpp \
$(SD)/HF/HartreeFock.hpp \
$(SD)/IO/UserInput.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Physics/AtomData.hpp \
$(SD)/Physics/NuclearPotentials.hpp \
$(SD)/Physics/Parametric_potentials.hpp \
$(SD)/Physics/PhysConst_constants.hpp \
$(SD)/Wavefunction/Wavefunction.hpp
	$(COMP)

$(BD)/Module_matrixElements.o: $(SD)/Modules/Module_matrixElements.cpp \
$(SD)/Modules/Module_matrixElements.hpp \
$(SD)/DiracOperator/DiracOperator.hpp \
$(SD)/DiracOperator/Operators.hpp \
$(SD)/HF/ExternalField.hpp \
$(SD)/IO/ChronoTimer.hpp \
$(SD)/IO/UserInput.hpp \
$(SD)/MBPT/CorrelationPotential.hpp \
$(SD)/MBPT/DiagramRPA.hpp \
$(SD)/Physics/NuclearPotentials.hpp \
$(SD)/Physics/PhysConst_constants.hpp \
$(SD)/Physics/RadiativePotential.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp \
$(SD)/Wavefunction/Wavefunction.hpp
	$(COMP)

$(BD)/Module_pnc.o: $(SD)/Modules/Module_pnc.cpp \
$(SD)/Modules/Module_pnc.hpp \
$(SD)/DiracOperator/DiracOperator.hpp \
$(SD)/DiracOperator/Operators.hpp \
$(SD)/HF/ExternalField.hpp \
$(SD)/HF/MixedStates.hpp \
$(SD)/IO/UserInput.hpp \
$(SD)/Physics/NuclearData.hpp \
$(SD)/Physics/NuclearPotentials.hpp \
$(SD)/Wavefunction/Wavefunction.hpp
	$(COMP)
