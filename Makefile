## R CMD check complains if we use naked "R" or "Rscript", but
## "$(R_HOME)/bin/Rscript" doesn't necessarily get the right path ...
R = R
Rscr = Rscript

vignettes: vignettes/fitode.pdf vignettes/fitode-ms.pdf vignettes/fitode.bib
	mv vignettes/fitode*.pdf inst/doc

%.pdf: %.Rnw
	$(Rscr) -e "knitr::knit2pdf('$*.Rnw')"
	rm $*.tex

%.tex: %.Rnw
	$(Rscr) -e "knitr::knit('$*.Rnw')"

%.R: %.Rnw
	$(Rscr) -e "knitr::purl('$*.Rnw')"

%.tex: %.R
	time $(R) CMD BATCH --vanilla $*.R

%.pdf: %.tex
	texi2pdf $<

clean:
	rm -f unnamed*
	rm -f *.aux *.tikz *.out *.log *.blg *.toc *.bbl *~ \#*

fresh: clean
	rm -f fitode-ms.pdf fitode.pdf
	rm cache/*
