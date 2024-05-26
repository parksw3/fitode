RSCRIPT="Rscript --vanilla -e"
brauer-ms.pdf: brauer-ms.Rnw brauer-ms.bib
	$(RSCRIPT) "knitr::knit2pdf('brauer-ms.Rnw')"

post: brauer-ms.pdf
	scp brauer-ms.pdf ms.mcmaster.ca:~/public_html/misc/fitode_EarnParkBolker.pdf

clean:
	rm -f *.log *.aux *.bbl *.blg *.out *~

diff:
	git checkout 1ed1332e
	$(RSCRIPT) "knitr::knit('brauer-ms.Rnw')"
	mv brauer-ms.tex brauer-ms-old.tex
	git checkout master
	$(RSCRIPT) "knitr::knit('brauer-ms.Rnw')"
	latexdiff brauer-ms-old.tex brauer-ms.tex >diff.tex
	texi2dvi -p diff.tex
