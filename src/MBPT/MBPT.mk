# Dependencies for MBPT

$(BD)/CorrelationPotential.o: $(SD)/MBPT/CorrelationPotential.cpp \
$(SD)/MBPT/CorrelationPotential.hpp \
$(SD)/Angular/Angular_tables.hpp \
$(SD)/Coulomb/Coulomb.hpp \
$(SD)/Coulomb/YkTable.hpp \
$(SD)/IO/FRW_fileReadWrite.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Maths/Interpolator.hpp \
$(SD)/Maths/LinAlg_MatrixVector.hpp \
$(SD)/Physics/AtomData.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp \
$(SD)/Wavefunction/Wavefunction.hpp
	$(COMP)
