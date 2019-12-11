# Dependencies for Modules

$(BD)/Module_runModules.o: $(SD)/Modules/Module_runModules.cpp \
$(SD)/Modules/Module_runModules.hpp $(SD)/DMionisation/Module_atomicKernal.hpp \
$(SD)/Dirac/DiracOperator.hpp $(SD)/Dirac/Operators.hpp \
$(SD)/Dirac/Wavefunction.hpp $(SD)/IO/UserInput.hpp \
$(SD)/Modules/Module_fitParametric.hpp $(SD)/Modules/Module_matrixElements.hpp \
$(SD)/Physics/PhysConst_constants.hpp $(SD)/HF/HartreeFockClass.hpp
	$(COMP)

$(BD)/Module_fitParametric.o: $(SD)/Modules/Module_fitParametric.cpp \
$(SD)/Modules/Module_fitParametric.hpp $(SD)/Dirac/Wavefunction.hpp \
$(SD)/HF/HartreeFockClass.hpp $(SD)/IO/UserInput.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Physics/NuclearPotentials.hpp $(SD)/Physics/Parametric_potentials.hpp \
$(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(BD)/Module_matrixElements.o: $(SD)/Modules/Module_matrixElements.cpp \
$(SD)/Modules/Module_matrixElements.hpp $(SD)/Dirac/DiracOperator.hpp \
$(SD)/Dirac/Operators.hpp $(SD)/Dirac/Wavefunction.hpp $(SD)/IO/UserInput.hpp \
$(SD)/Physics/NuclearPotentials.hpp $(SD)/Physics/PhysConst_constants.hpp
	$(COMP)

$(BD)/Module_pnc.o: $(SD)/Modules/Module_pnc.cpp \
$(SD)/Modules/Module_runModules.hpp  \
$(SD)/Dirac/DiracOperator.hpp $(SD)/Dirac/Operators.hpp \
$(SD)/Dirac/Wavefunction.hpp $(SD)/IO/UserInput.hpp
	$(COMP)
