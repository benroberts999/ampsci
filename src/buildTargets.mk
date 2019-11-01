# All makefile targets

################################################################################
#Allow exectuables to be placed in another directory:
ALLEXES = $(addprefix $(XD)/, \
 hartreeFock wigner dmeXSection periodicTable \
)

DEFAULTEXES = $(addprefix $(XD)/, \
 hartreeFock wigner periodicTable \
)

#Default make rule:
all: checkObj checkXdir $(DEFAULTEXES) doneMessage

################################################################################

# Dependencies for each object target, broken by sub-directories:
include $(SD)/Adams/Adams.mk
include $(SD)/Dirac/Dirac.mk
include $(SD)/DMionisation/DMionisation.mk
include $(SD)/HF/HF.mk
include $(SD)/IO/IO.mk
include $(SD)/Maths/Maths.mk
include $(SD)/Modules/Modules.mk
include $(SD)/Physics/Physics.mk
# Dependencies for each 'main' object target
include $(SD)/main.mk

################################################################################
# Just to save typing: Many programs depend on these combos:

BASE = $(addprefix $(OD)/, \
 Adams_bound.o Wavefunction.o DiracSpinor.o AtomInfo.o Nuclear.o Grid.o \
)

HF = $(addprefix $(OD)/, \
 HartreeFockClass.o CoulombIntegrals.o Parametric_potentials.o \
 Adams_Greens.o \
)

CNTM = $(addprefix $(OD)/, \
 Adams_continuum.o ContinuumOrbitals.o \
)

MODS = $(addprefix $(OD)/, \
 Module_runModules.o Module_atomicKernal.o AKF_akFunctions.o \
 Module_matrixElements.o Module_fitParametric.o \
)

################################################################################
# Link + build all final programs

$(XD)/hartreeFock: $(BASE) $(HF) $(CNTM) $(OD)/hartreeFock.o \
$(OD)/UserInput.o $(MODS)
	$(LINK)

$(XD)/dmeXSection: $(BASE) $(CNTM) $(HF) $(OD)/dmeXSection.o \
$(OD)/AKF_akFunctions.o $(OD)/StandardHaloModel.o
	$(LINK)

$(XD)/wigner: $(OD)/wigner.o
	$(LINK)

$(XD)/periodicTable: $(OD)/periodicTable.o $(OD)/AtomInfo.o $(OD)/Nuclear.o
	$(LINK)

################################################################################
################################################################################

checkObj:
	@if [ ! -d $(OD) ]; then \
	  echo '\n ERROR: Directory: '$(OD)' doesnt exist - please create it!\n'; \
	  false; \
	fi

checkXdir:
	@if [ ! -d $(XD) ]; then \
		echo '\n ERROR: Directory: '$(XD)' doesnt exist - please create it!\n'; \
		false; \
	fi

doneMessage:
		@echo 'done'

.PHONY: clean do_the_chicken_dance checkObj checkXdir
clean:
	rm -f $(ALLEXES) $(OD)/*.o
do_the_chicken_dance:
	@echo 'Why would I do that?'
