# Dependencies for ADAMS

$(OD)/Adams_bound.o: $(SD)/Adams/Adams_bound.cpp $(SD)/Adams/Adams_bound.hpp \
$(SD)/Adams/Adams_coefs.hpp $(SD)/Dirac/DiracSpinor.hpp \
$(SD)/Adams/DiracODE.hpp $(SD)/Maths/Grid.hpp $(SD)/Maths/Matrix_linalg.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp
	$(COMP)

$(OD)/Adams_continuum.o: $(SD)/Adams/Adams_continuum.cpp \
$(SD)/Adams/Adams_continuum.hpp $(SD)/Adams/Adams_bound.hpp \
$(SD)/Adams/DiracODE.hpp $(SD)/Maths/Grid.hpp $(SD)/Adams/Adams_coefs.hpp
	$(COMP)

$(OD)/Adams_Greens.o: $(SD)/Adams/Adams_Greens.cpp \
$(SD)/Adams/Adams_Greens.hpp $(SD)/Adams/Adams_bound.hpp \
$(SD)/Dirac/DiracSpinor.hpp $(SD)/Adams/DiracODE.hpp $(SD)/Maths/Grid.hpp \
$(SD)/Maths/NumCalc_quadIntegrate.hpp $(SD)/Adams/Adams_coefs.hpp
	$(COMP)
