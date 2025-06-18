MAKEFLAGS += --no-print-directory

all:
	@make html
	@echo "Finished rendering."
	@echo "To continously rebuild on change, use 'make watch'."

html:
	@Rscript -e "rmarkdown::render('reproduce_results.Rmd', quiet=TRUE)"

watch:
	@watchexec -e Rmd -w . --clear -- make html

