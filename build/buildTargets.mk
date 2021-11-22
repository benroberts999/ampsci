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
all: $(SD)/git.info checkObj checkXdir $(DEFAULTEXES)

################################################################################
# Automatically generate dependency files for each cpp file, + compile:
# I make all files depend on '$(SD)/git.info'. This achieves two things:
# (1) It forces git.info to be built first (even in a parallel build)
# (2) It forces a clean make when changing branches
# The '|' means 'order only' pre-req

# All the files that in in src/{subsir}/.cpp don't have a main()
# Compile them into objects; reflect the src/ directory tree.
$(BD)/*/%.o: $(SD)/*/%.cpp | $(SD)/git.info
	@mkdir -p $(@D)
	$(COMP)

# All the files that in in src/.cpp *do* have a main(), compile
$(BD)/%.o: $(SD)/%.cpp $(SD)/git.info
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

################################################################################
################################################################################

# Check to see if this is a git repo, so $(SD)/git.info will work even if not
ifneq ("$(wildcard .git/HEAD)","")
  GIT_FILES = .git/HEAD .git/index
endif
# Create the 'git.info' file (c++ header file)
$(SD)/git.info: $(GIT_FILES)
	@echo Git Files: $(GIT_FILES)
	@echo "// git.info: auto-generated file" > $@
	@echo "#pragma once" >> $@
	@echo "namespace GitInfo {" >> $@
	@echo Git version: $(shell git rev-parse --short HEAD 2>/dev/null)
	@echo "const char *gitversion = \"$(shell git rev-parse --short HEAD 2>/dev/null)\";" >> $@
	@echo "const char *gitbranch = \"$(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)\";" >> $@
	@echo Git branch : $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)
	@echo "} // namespace GitInfo" >> $@

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

.PHONY: clean docs doxy do_the_chicken_dance checkObj checkXdir
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
	make docs
	doxygen ./src/Doxyfile
	( cd ./doc/latex && make )
	cp ./doc/latex/refman.pdf ./doc/documentation.pdf
	cp ./doc/ampsci.pdf ./docs/ampsci.pdf 2>/dev/null || :
	( cd ./doc/latex && make clean)
do_the_chicken_dance:
	@echo $(OBJS)
	@echo 'Why would I do that?'
clang_format:
	clang-format -i src/*.cpp src/*/*.cpp src/*/*.hpp
