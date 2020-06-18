# Dependencies for Maths

$(BD)/Grid.o: $(SD)/Maths/Grid.cpp \
$(SD)/Maths/Grid.hpp
	$(COMP)

$(BD)/LinAlg_MatrixVector.o: $(SD)/Maths/LinAlg_MatrixVector.cpp \
$(SD)/Maths/LinAlg_MatrixVector.hpp \
$(SD)/IO/SafeProfiler.hpp
	$(COMP)
