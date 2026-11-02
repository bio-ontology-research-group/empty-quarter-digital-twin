# Build the supplement first so that the main text can resolve its references to supplementary items (package xr).
PDFLATEX = pdflatex -interaction=nonstopmode -halt-on-error

all: supplementary.pdf main.pdf

supplementary.pdf: FORCE
	$(PDFLATEX) supplementary.tex
	$(PDFLATEX) supplementary.tex

main.pdf: supplementary.pdf FORCE
	$(PDFLATEX) main.tex
	bibtex main
	$(PDFLATEX) main.tex
	$(PDFLATEX) main.tex

FORCE:
.PHONY: all FORCE
