
doc=pdflatex
bib=bibtex
delete=rm -f

thesis.pdf: *.tex thesis.bib VTthesis.cls Makefile
	$(doc) thesis
	$(bib) thesis
	$(doc) thesis
	$(doc) thesis

clean:
	$(delete) thesis.bbl
	$(delete) thesis.blg
	$(delete) thesis.nlo
	$(delete) thesis.lof
	$(delete) thesis.log
	$(delete) thesis.lot
	$(delete) thesis.toc
	$(delete) thesis.out
	$(delete) thesis.aux

