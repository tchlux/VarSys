
doc=pdflatex
bib=bibtex
delete=rm -f

thesis.pdf: *.tex thesis.bib VTthesis.cls Makefile
	$(doc) thesis > /dev/null
	$(bib) thesis > /dev/null
	$(doc) thesis > /dev/null
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

