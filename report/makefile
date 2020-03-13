all: _book/resdoc.pdf

_book/resdoc.pdf: index.Rmd 01-introduction.Rmd 02-methods.Rmd 03-mse-results.Rmd 03-mse-results2.Rmd 04-discussion.Rmd 05-appendix-om-description.Rmd 06-appendix-MP-descriptions.Rmd 07-appendix-biological-data.Rmd 08-appendix-catch-data.Rmd 09-appendix-survey-data.Rmd 10-appendix-management-changes.Rmd 11-references.Rmd bib/mpframework-refs.bib bib/refs.bib
	Rscript -e "bookdown::render_book('index.Rmd')"

continuous:
	while true; do make --silent; sleep 0.2; done

bibclean:
	bibtool -s --print.use.tab=Off --print.line.length=80 -i bib/refs.bib -o bib/refs.bib
	bibtool -s --print.use.tab=Off --print.line.length=80 -i bib/mpframework-refs.bib -o refs/mpframework-refs.bib
