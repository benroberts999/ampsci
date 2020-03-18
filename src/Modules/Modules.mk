# Dependencies for Modules

MODULELIST = $(addprefix $(BD)/, \
 Module_runModules.o Module_atomicKernal.o AKF_akFunctions.o \
 Module_matrixElements.o Module_fitParametric.o Module_pnc.o Module_tests.o \
)

$(BD)/Module_runModules.o: $(SD)/Modules/Module_runModules.cpp \
$(SD)/Modules/Module_runModules.hpp $(SD)/DMionisation/Module_atomicKernal.hpp \
$(SD)/DiracOperator/DiracOperator.hpp $(SD)/DiracOperator/Operators.hpp \
$(SD)/Wavefunction/Hamiltonian.hpp \
$(SD)/Wavefunction/Wavefunction.hpp $(SD)/IO/UserInput.hpp \
$(SD)/Modules/Module_fitParametric.hpp $(SD)/Modules/Module_matrixElements.hpp \
$(SD)/Physics/PhysConst_constants.hpp $(SD)/HF/HartreeFockClass.hpp
	$(COMP)

$(BD)/Module_tests.o: \
$(SD)/Modules/Module_tests.cpp $(SD)/Modules/Module_tests.hpp \
$(SD)/DiracOperator/Operators.hpp $(SD)/Wavefunction/Hamiltonian.hpp \
$(SD)/Wavefunction/Wavefunction.hpp $(SD)/IO/UserInput.hpp \
$(SD)/Physics/PhysConst_constants.hpp $(SD)/HF/HartreeFockClass.hpp
	$(COMP)

$(BD)/Module_fitParametric.o: $(SD)/Modules/Module_fitParametric.cpp \
$(SD)/Modules/Module_fitParametric.hpp $(SD)/Wavefunction/Wavefunction.hpp \
$(SD)/HF/HartreeFockClass.hpp $(SD)/IO/UserInput.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Physics/NuclearPotentials.hpp $(SD)/Physics/Parametric_potentials.hpp \
$(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(BD)/Module_matrixElements.o: $(SD)/Modules/Module_matrixElements.cpp \
$(SD)/Modules/Module_matrixElements.hpp $(SD)/DiracOperator/DiracOperator.hpp \
$(SD)/DiracOperator/Operators.hpp $(SD)/Wavefunction/Wavefunction.hpp $(SD)/IO/UserInput.hpp \
$(SD)/Physics/NuclearPotentials.hpp $(SD)/Physics/PhysConst_constants.hpp \
$(SD)/HF/ExternalField.hpp
	$(COMP)

$(BD)/Module_pnc.o: $(SD)/Modules/Module_pnc.cpp $(SD)/Modules/Module_pnc.hpp\
$(SD)/HF/ExternalField.hpp $(SD)/DiracOperator/Operators.hpp \
$(SD)/Wavefunction/Wavefunction.hpp $(SD)/IO/UserInput.hpp
	$(COMP)
