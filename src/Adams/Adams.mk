# Dependencies for ADAMS

$(BD)/Adams_bound.o: $(SD)/Adams/Adams_bound.cpp $(SD)/Adams/Adams_bound.hpp \
$(SD)/Adams/Adams_coefs.hpp $(SD)/Dirac/DiracSpinor.hpp \
$(SD)/Adams/DiracODE.hpp $(SD)/Maths/Grid.hpp $(SD)/Maths/Matrix_linalg.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/IO/SafeProfiler.hpp
	$(COMP)

$(BD)/Adams_continuum.o: $(SD)/Adams/Adams_continuum.cpp \
$(SD)/Adams/Adams_continuum.hpp $(SD)/Adams/Adams_bound.hpp \
$(SD)/Adams/DiracODE.hpp $(SD)/Maths/Grid.hpp $(SD)/Adams/Adams_coefs.hpp
	$(COMP)

$(BD)/Adams_Greens.o: $(SD)/Adams/Adams_Greens.cpp \
$(SD)/Adams/Adams_Greens.hpp $(SD)/Adams/Adams_bound.hpp \
$(SD)/Dirac/DiracSpinor.hpp $(SD)/Adams/DiracODE.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/Adams/Adams_coefs.hpp \
$(SD)/IO/SafeProfiler.hpp
	$(COMP)
