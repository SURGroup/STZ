# Utility to automatically build LaTeX files
compile=latexmk -pdf -silent

# Utility to automatically remove all intermediate files (e.g. ".aux", ".log")
cleanup=latexmk -c

d_pdf=overview.pdf bar_coords.pdf

all: $(d_pdf)

%.pdf: %.tex elas.bib
	$(compile) $<
	$(cleanup) $<

bar_coords.pdf: problem_1.jpg problem_2.jpg

clean:
	rm -rf $(d_pdf)

.PHONY: clean
