# Dependencies for IO

$(BD)/UserInput.o: $(SD)/IO/UserInput.cpp $(SD)/IO/UserInput.hpp \
$(SD)/IO/FRW_fileReadWrite.hpp
	$(COMP)
