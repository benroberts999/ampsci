# Dependencies for IO

$(BD)/UserInput.o: $(SD)/IO/UserInput.cpp $(SD)/IO/UserInput.hpp \
$(SD)/IO/FileIO_fileReadWrite.hpp
	$(COMP)
