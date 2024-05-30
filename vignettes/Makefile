RSCRIPT=Rscript --vanilla -e
OLDVERSION=1ed1332e
brauer-ms.pdf: brauer-ms.Rnw brauer-ms.bib
	$(RSCRIPT) "knitr::knit2pdf('brauer-ms.Rnw')"

FIGS=Bombay-figure-1.pdf  Bombay-figure-2.pdf  phila-figure-1.pdf  plot_stochastic_SIR-1.pdf  plot_stochastic_SIR-2.pdf

brauer-figs: brauer-ms.pdf
	mkdir -p ../outputs/brauer-figs
	cd figure; cp $(FIGS) ../../outputs/brauer-figs

post: brauer-ms.pdf
	scp brauer-ms.pdf ms.mcmaster.ca:~/public_html/misc/fitode_EarnParkBolker.pdf

clean:
	rm -f *.log *.aux *.bbl *.blg *.out *~

## make PDF of diff. This will fail if repo is not in a clean state
## (i.e., commit changes first ...)
## Makefile has to be removed when reverting to HEAD, not 100% clear why
##  this is the case
diff: brauer-ms.Rnw brauer-ms.bib
	git checkout $(OLDVERSION)
	$(RSCRIPT) "knitr::knit('brauer-ms.Rnw')"
	mv brauer-ms.tex brauer-ms-old.tex
	git checkout master
	$(RSCRIPT) "knitr::knit('brauer-ms.Rnw')"
	latexdiff brauer-ms-old.tex brauer-ms.tex >diff.tex
	texi2dvi -p diff.tex
