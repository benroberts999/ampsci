LATEX_CMD=pdflatex

all: diracSCAS.pdf

pdf: diracSCAS.pdf

diracSCAS.pdf: clean diracSCAS.tex
	$(LATEX_CMD) diracSCAS
	#bibtex diracSCAS #
	$(LATEX_CMD) diracSCAS
	latex_count=8 ; \
	while egrep -s 'Rerun (LaTeX|to get cross-references right)' diracSCAS.log && [ $$latex_count -gt 0 ] ;\
	    do \
	      echo "Rerunning latex...." ;\
	      $(LATEX_CMD) diracSCAS ;\
	      latex_count=`expr $$latex_count - 1` ;\
	    done
	#bibtex diracSCAS #
	$(LATEX_CMD) diracSCAS


clean:
	rm -f *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *.brf *.blg diracSCAS.pdf
