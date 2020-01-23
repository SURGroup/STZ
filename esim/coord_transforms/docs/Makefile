# Utility to automatically build LaTeX files
compile=latexmk -pdf -silent

# Utility to automatically remove all intermediate files (e.g. ".aux", ".log")
cleanup=latexmk -c

d_pdf=transformations.pdf

all: $(d_pdf)

%.pdf: %.tex
	$(compile) $<
	$(cleanup) $<

clean:
	rm -rf $(d_pdf)

.PHONY: clean
