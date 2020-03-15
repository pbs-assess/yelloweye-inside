all: _book/resdoc.pdf

_book/resdoc.pdf: index.Rmd 00-introduction.Rmd 01-decision-context.Rmd 02-objectives.Rmd 03-om.Rmd 04-mp.Rmd 05-simulation.Rmd 06-results.Rmd 07-discussion.Rmd 08-app-biological-data.Rmd 09-app-survey-data.Rmd 10-app-catch-data.Rmd 11-app-om.Rmd 12-app-mp.Rmd 98-environment.Rmd 99-references.Rmd bib/mpframework-refs.bib bib/refs.bib
	Rscript -e "bookdown::render_book('index.Rmd')"

continuous:
	while true; do make --silent; sleep 0.2; done

bibclean:
	bibtool -s --print.use.tab=Off --print.line.length=80 -i bib/mpframework-refs.bib -o bib/mpframework-refs.bib
	bibtool -s --print.use.tab=Off --print.line.length=80 -i bib/refs.bib -o bib/refs.bib
