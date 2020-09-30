LATEX_CMD=pdflatex

all: ampsci.pdf

pdf: ampsci.pdf

ampsci.pdf: clean ampsci.tex
	$(LATEX_CMD) ampsci
	#bibtex ampsci #
	$(LATEX_CMD) ampsci
	latex_count=8 ; \
	while egrep -s 'Rerun (LaTeX|to get cross-references right)' ampsci.log && [ $$latex_count -gt 0 ] ;\
	    do \
	      echo "Rerunning latex...." ;\
	      $(LATEX_CMD) ampsci ;\
	      latex_count=`expr $$latex_count - 1` ;\
	    done
	#bibtex ampsci #
	$(LATEX_CMD) ampsci


clean:
	rm -f *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *.brf *.blg ampsci.pdf
