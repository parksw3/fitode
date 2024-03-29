brauer-ms.pdf: brauer-ms.Rnw brauer-ms.bib
	Rscript -e "knitr::knit2pdf('brauer-ms.Rnw')"

post: brauer-ms.pdf
	scp brauer-ms.pdf ms.mcmaster.ca:~/public_html/misc/fitode_EarnParkBolker.pdf
