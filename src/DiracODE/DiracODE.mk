# Dependencies for DiracODE

$(BD)/Adams_bound.o: $(SD)/DiracODE/Adams_bound.cpp \
$(SD)/DiracODE/Adams_bound.hpp \
$(SD)/DiracODE/Adams_coefs.hpp \
$(SD)/DiracODE/DiracODE.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Maths/LinAlg_MatrixVector.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp
	$(COMP)

$(BD)/Adams_continuum.o: $(SD)/DiracODE/Adams_continuum.cpp \
$(SD)/DiracODE/Adams_continuum.hpp \
$(SD)/DiracODE/Adams_bound.hpp \
$(SD)/DiracODE/DiracODE.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp \
$(SD)/Maths/Grid.hpp
	$(COMP)

$(BD)/Adams_Greens.o: $(SD)/DiracODE/Adams_Greens.cpp \
$(SD)/DiracODE/Adams_Greens.hpp \
$(SD)/DiracODE/Adams_bound.hpp \
$(SD)/DiracODE/DiracODE.hpp \
$(SD)/IO/SafeProfiler.hpp \
$(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp \
$(SD)/Wavefunction/DiracSpinor.hpp
	$(COMP)
