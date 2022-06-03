# All makefile targets

################################################################################
#Allow exectuables to be placed in another directory:
ALLEXES = $(addprefix $(XD)/, \
 ampsci unitTests wigner dmeXSection periodicTable \
)

DEFAULTEXES = $(addprefix $(XD)/, \
 ampsci periodicTable unitTests \
)

#Default make rule:
all: GitInfo checkObj checkXdir $(DEFAULTEXES)

################################################################################
# Automatically generate dependency files for each cpp file, + compile:

# All the files that in in src/{subsir}/.cpp don't have a main()
# Compile them into objects; reflect the src/ directory tree.
$(BD)/*/%.o: $(SD)/*/%.cpp
	@mkdir -p $(@D)
	$(COMP)

# All the files that in in src/.cpp *do* have a main(), compile
$(BD)/%.o: $(SD)/%.cpp
	@mkdir -p $(@D)
	$(COMP)

# include the dependency files
-include $(BD)/*.d $(BD)/*/*.d

################################################################################
# List all objects in sub-directories (i.e., that don't conatin a main())
# e.g., for each src/sub_dir/file.cpp -> build/sub_dir/file.o
OBJS = $(subst $(SD),$(BD),$(subst .cpp,.o,$(wildcard $(SD)/*/*.cpp)))

################################################################################
# Link + build all final programs

$(XD)/ampsci: $(BD)/ampsci.o $(OBJS)
	$(LINK)

$(XD)/unitTests: $(BD)/unitTests.o $(OBJS)
	$(LINK)

$(XD)/dmeXSection: $(BD)/dmeXSection.o $(OBJS)
	$(LINK)

$(XD)/wigner: $(BD)/wigner.o
	$(LINK)

$(XD)/periodicTable: $(BD)/periodicTable.o $(BD)/Physics/AtomData.o \
$(BD)/Physics/NuclearData.o
	$(LINK)

# Add git version info to compile flags
NOW:=$(shell date +%Y-%m-%d' '%H:%M' '%Z 2>/dev/null)
CXXFLAGS+=-D GITREVISION="$(shell git rev-parse --short HEAD 2>/dev/null)"
CXXFLAGS+=-D GITBRANCH="$(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)"
CXXFLAGS+=-D GITMODIFIED="$(shell git status -s 2>/dev/null)"
CXXFLAGS+=-D CXXVERSION="$(shell $(CXX) --version 2>/dev/null | sed -n '1p' | sed s/'('/'['/ | sed s/')'/']'/)"
CXXFLAGS+=-D COMPTIME="$(NOW)"

################################################################################
################################################################################

# CPrint some git history to screen
SHELL=/bin/bash
string=$(shell $(CXX) --version 2>/dev/null | sed -n '1p' | sed s/'('/'['/ | sed s/')'/']'/)
GitInfo:
	@echo Git Branch: $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)
	@echo Git Revision: $(shell git rev-parse --short HEAD 2>/dev/null)
	@echo Modified: $(shell git status -s 2>/dev/null)

checkObj:
	@if [ ! -d $(BD) ]; then \
	  echo '\n ERROR: Directory: '$(BD)' doesnt exist - please create it!\n'; \
	  false; \
	fi

checkXdir:
	@if [ ! -d $(XD) ]; then \
		echo '\n ERROR: Directory: '$(XD)' doesnt exist - please create it!\n'; \
		false; \
	fi

.PHONY: clean docs doxy do_the_chicken_dance GitInfo checkObj checkXdir
clean:
	rm -f -v $(ALLEXES)
	rm -rf -v $(BD)/*.o $(BD)/*.d $(BD)/*/
# Make the 'ampsci.pdf' physics documentation
docs:
	( cd ./doc/tex && make )
	cp ./doc/tex/ampsci.pdf ./doc/ampsci.pdf
	( cd ./doc/tex && make clean)
# Make the doxygen code documentation
doxy:
	doxygen ./src/Doxyfile
	( cd ./doc/latex && make )
	cp ./doc/latex/refman.pdf ./doc/documentation.pdf
	make docs
	cp ./doc/ampsci.pdf ./docs/ampsci.pdf 2>/dev/null || :
	( cd ./doc/latex && make clean)
do_the_chicken_dance:
	@echo 'Why would I do that?'
clang_format:
	clang-format -i src/*.cpp src/*/*.cpp src/*/*.hpp
