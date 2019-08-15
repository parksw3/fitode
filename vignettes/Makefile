fitode-ms.pdf: fitode.bib

%.tex: %.Rnw
	Rscript -e "library(knitr); knit('$*.Rnw')"

%.R: %.Rnw
	Rscript -e "library(knitr); purl('$*.Rnw')"

%.tex: %.R
	time R CMD BATCH --vanilla $*.R

%.pdf: %.tex
	texi2pdf $<
clean:
	rm -f unnamed*
	rm -f *.{aux,tikz,out,log,blg,toc,bbl}
fresh: clean
	rm -f fitode-ms.pdf fitode.pdf