all: _book/resdoc.pdf bibclean dropbox

_book/resdoc.pdf: index.Rmd 00-introduction.Rmd 01-decision-context.Rmd 02-objectives.Rmd 03-om.Rmd 04-mp.Rmd 05-simulation.Rmd 06-results.Rmd 07-discussion.Rmd 08-app-biological-data.Rmd 09-app-survey-data.Rmd 10-app-catch-data.Rmd 11-app-om.Rmd 12-app-mp.Rmd 98-environment.Rmd 99-references.Rmd mpframework-refs.bib refs.bib 14-app-cosewic.Rmd
	Rscript -e "bookdown::render_book('index.Rmd')"

trash-cache:
	rm -rf knitr-cache-pdf
	rm -rf knitr-figs-pdf

continuous:
	while true; do make _book/resdoc.pdf --silent; sleep 0.2; done

bibclean:
	bibtool -s --print.use.tab=Off --print.line.length=80 -i mpframework-refs.bib -o mpframework-refs.bib
	bibtool -s --print.use.tab=Off --print.line.length=80 -i refs.bib -o refs.bib

dropbox:
	cp -r ../figs ~/Dropbox/temp/ye-figs/
	cp -r ../mse/figures ~/Dropbox/temp/ye-figs/mse/
