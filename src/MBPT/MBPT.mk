# Dependencies for MBPT

$(BD)/CorrelationPotential.o: $(SD)/MBPT/CorrelationPotential.cpp \
$(SD)/MBPT/CorrelationPotential.hpp \
$(SD)/Angular/Angular_tables.hpp \
$(SD)/Coulomb/Coulomb.hpp \
$(SD)/Coulomb/YkTable.hpp \
$(SD)/HF/HartreeFock.hpp \
$(SD)/IO/FRW_fileReadWrite.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/MBPT/GreenMatrix.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Maths/Interpolator.hpp \
$(SD)/Maths/LinAlg_MatrixVector.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp
	$(COMP)

$(BD)/GoldstoneSigma2.o: $(SD)/MBPT/GoldstoneSigma2.cpp \
$(SD)/MBPT/GoldstoneSigma2.hpp \
$(SD)/Angular/Angular_tables.hpp \
$(SD)/Coulomb/Coulomb.hpp \
$(SD)/Coulomb/YkTable.hpp \
$(SD)/IO/FRW_fileReadWrite.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/MBPT/CorrelationPotential.hpp \
$(SD)/MBPT/GreenMatrix.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Maths/LinAlg_MatrixVector.hpp
	$(COMP)

$(BD)/FeynmanSigma.o: $(SD)/MBPT/FeynmanSigma.cpp \
$(SD)/MBPT/FeynmanSigma.hpp \
$(SD)/Angular/Angular_tables.hpp \
$(SD)/Coulomb/Coulomb.hpp \
$(SD)/Coulomb/YkTable.hpp \
$(SD)/DiracODE/DiracODE.hpp \
$(SD)/HF/HartreeFock.hpp \
$(SD)/IO/FRW_fileReadWrite.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/MBPT/CorrelationPotential.hpp \
$(SD)/MBPT/GreenMatrix.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Maths/LinAlg_MatrixVector.hpp
	$(COMP)

$(BD)/DiagramRPA.o: $(SD)/MBPT/DiagramRPA.cpp \
$(SD)/MBPT/DiagramRPA.hpp \
$(SD)/Angular/Angular_369j.hpp \
$(SD)/Angular/Angular_tables.hpp \
$(SD)/Coulomb/Coulomb.hpp \
$(SD)/Coulomb/YkTable.hpp \
$(SD)/IO/FRW_fileReadWrite.hpp \
$(SD)/DiracOperator/DiracOperator.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp
	$(COMP)
