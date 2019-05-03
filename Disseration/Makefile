
doc=pdflatex
bib=bibtex
delete=rm -f

prelim.pdf: *.tex prelim.bib VTthesis.cls Makefile
	$(doc) prelim
	$(bib) prelim
	$(doc) prelim
	$(doc) prelim

clean:
	$(delete) prelim.bbl
	$(delete) prelim.blg
	$(delete) prelim.nlo
	$(delete) prelim.lof
	$(delete) prelim.log
	$(delete) prelim.lot
	$(delete) prelim.toc
	$(delete) prelim.out
	$(delete) prelim.aux

