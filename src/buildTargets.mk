# BUILD TARGETS:

DEFAULTEXES := ampsci
ifneq ($(MODE),release)
  DEFAULTEXES += tests
endif

ALLEXES := ampsci tests

.DEFAULT_GOAL := all

all:  $(DEFAULTEXES)
full: $(ALLEXES)


.PHONY: all full clean remove_junk docs doxy clang_format pre_commit includes force do_the_chicken_dance

clean:
	rm -fv $(ALLEXES)
	rm -rfv $(BUILD)/* $(BUILD)/*/
	rm -fv *deleteme*

remove_junk:
	rm -f -v *deleteme*

# Builds the ampsci.pdf (physics) docs
docs:
	( cd ./doc/tex && make )
	cp ./doc/tex/ampsci.pdf ./doc/ampsci.pdf
	( cd ./doc/tex && make clean)

# Builds code documentation using doxygen, if available
doxy:
	doxygen ./doc/doxygen/Doxyfile
	( cd ./doc/latex && make )
	cp ./doc/latex/refman.pdf ./doc/documentation.pdf
	make docs
	cp ./doc/ampsci.pdf ./docs/ampsci.pdf 2>/dev/null || :
	( cd ./doc/latex && make clean)

do_the_chicken_dance:
	@echo 'Why would I do that?'

# Don't run clang format on external libraries
CLANG_FORMAT_EXCLUDES := $(SRC)/fmt $(SRC)/catch2 $(SRC)/json

clang_format:
	@echo "Running clang format (whole project)"
	find $(SRC)/ \( $(foreach dir,$(CLANG_FORMAT_EXCLUDES),-path "$(dir)" -o ) -false \) -prune -o \( -iname '*.cpp' -o -iname '*.hpp' -o -iname '*.ipp' \) -print | xargs -r clang-format -i -verbose

# Makes the includes
includes:
	./$(SRC)/build_includes.sh
	@echo "Running clang format (on includes)"
	find $(SRC)/ -iname 'include.hpp' | xargs clang-format -i -verbose

license-year:
	@year=$$(date +%Y); \
	sed -i.bak -E "s/(Copyright \(c\) )([0-9]{4})(–[0-9]{4})?/\1\2–$${year}/" LICENSE; \
	rm -f LICENSE.bak; \
	echo "Updated license year to $${year}"


pre_commit:
	make includes
	make clang_format
	make ampsci tests
	make test
	make remove_junk