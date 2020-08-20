# All makefile targets

################################################################################
#Allow exectuables to be placed in another directory:
ALLEXES = $(addprefix $(XD)/, \
 diracSCAS unitTests wigner dmeXSection periodicTable \
)

DEFAULTEXES = $(addprefix $(XD)/, \
 diracSCAS wigner periodicTable \
)

# if dev or debug build, make unitTests too
ifneq ($(Build),release)
  DEFAULTEXES+=$(XD)/unitTests
endif

#Default make rule:
all: $(SD)/git.info checkObj checkXdir $(DEFAULTEXES)

################################################################################
# Automatically generate dependency files for each cpp file, + compile:
# I make all files depend on '$(SD)/git.info'. This achieves two things:
# (1) It forces git.info to be built first (even in a parallel build)
# (2) It forces a clean make when changing branches

$(BD)/%.o: $(SD)/*/%.cpp $(SD)/git.info
	$(COMP)

$(BD)/%.o: $(SD)/%.cpp $(SD)/git.info
	$(COMP)

-include $(BD)/*.d

################################################################################
# Link + build all final programs

# List all objects in sub-directories (i.e., that don't conatin a main())
OBJS = $(addprefix $(BD)/,$(notdir $(subst .cpp,.o,$(wildcard $(SD)/*/*.cpp))))

$(XD)/diracSCAS: $(BD)/diracSCAS.o $(OBJS)
	$(LINK)

$(XD)/unitTests: $(BD)/unitTests.o $(OBJS)
	$(LINK)

$(XD)/dmeXSection: $(BD)/dmeXSection.o $(OBJS)
	$(LINK)

$(XD)/wigner: $(BD)/wigner.o
	$(LINK)

$(XD)/periodicTable: $(BD)/periodicTable.o $(BD)/AtomData.o \
$(BD)/NuclearData.o
	$(LINK)

################################################################################
################################################################################

# Check to see if this is a git repo, so $(SD)/git.info will work even if not
ifneq ("$(wildcard .git/HEAD)","")
  GIT_FILES = .git/HEAD .git/index
endif
# Make the 'git.info' file (c++ header file)
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
	rm -f $(ALLEXES) $(BD)/*.o $(BD)/*.d
# Make the 'diracSCAS.pdf' physics documentation
docs:
	( cd ./doc/tex && make )
	cp ./doc/tex/diracSCAS.pdf ./doc/diracSCAS.pdf
	( cd ./doc/tex && make clean)
# Make the doxygen code documentation
doxy:
	doxygen ./src/Doxyfile
	( cd ./doc/latex && make )
	cp ./doc/latex/refman.pdf ./doc/documentation.pdf
	cp ./doc/diracSCAS.pdf ./doc/html/diracSCAS.pdf 2>/dev/null || :
	( cd ./doc/latex && make clean)
do_the_chicken_dance:
	@echo 'Why would I do that?'
