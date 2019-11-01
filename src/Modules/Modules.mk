# Dependencies for Modules

$(OD)/Module_runModules.o: $(SD)/Modules/Module_runModules.cpp \
$(SD)/Modules/Module_runModules.hpp $(SD)/DMionisation/Module_atomicKernal.hpp \
$(SD)/Dirac/DiracOperator.hpp $(SD)/Dirac/Operators.hpp \
$(SD)/Dirac/Wavefunction.hpp $(SD)/IO/UserInput.hpp \
$(SD)/Modules/Module_fitParametric.hpp $(SD)/Modules/Module_matrixElements.hpp \
$(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(OD)/Module_fitParametric.o: $(SD)/Modules/Module_fitParametric.cpp \
$(SD)/Modules/Module_fitParametric.hpp $(SD)/Dirac/Wavefunction.hpp \
$(SD)/HF/HartreeFockClass.hpp $(SD)/IO/UserInput.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Physics/Nuclear.hpp $(SD)/Physics/Parametric_potentials.hpp \
$(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(OD)/Module_matrixElements.o: $(SD)/Modules/Module_matrixElements.cpp \
$(SD)/Modules/Module_matrixElements.hpp $(SD)/Dirac/DiracOperator.hpp \
$(SD)/Dirac/Operators.hpp $(SD)/Dirac/Wavefunction.hpp $(SD)/IO/UserInput.hpp \
$(SD)/Physics/Nuclear.hpp $(SD)/Physics/PhysConst_constants.hpp
	$(COMP)
