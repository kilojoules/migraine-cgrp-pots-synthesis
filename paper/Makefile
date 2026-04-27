PAPER = main

.PHONY: all clean watch

all: $(PAPER).pdf

$(PAPER).pdf: $(PAPER).tex refs.bib
	latexmk -pdf -bibtex -interaction=nonstopmode $(PAPER).tex

watch:
	latexmk -pdf -bibtex -pvc -interaction=nonstopmode $(PAPER).tex

clean:
	latexmk -C
	rm -f *.bbl *.run.xml *.synctex.gz
